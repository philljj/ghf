#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "basis.h"
#include "matrix.h"
#include "util.h"

#if !defined M_PI
  #define M_PI  3.14159265359
#endif

#define MAX_GEOM_FILE_SIZE 1024

static double get_r_diff_sq(const r_t a, const r_t b);
static r_t    get_Rp(const size_t mu, const size_t nu);

static bool init_basis_set(const char * file);
static bool init_atom_list(const char * file);
static void init_R_list(void);
static bool load_geom_file(const char * file);
static bool load_shell(const char * p);
static void dump_basis_set(void);
static void dump_atom_list(void);
static const char * skip_line(const char * p);
static const char * skip_space(const char * p);
static const char * skip_float(const char * p);

static atom_t       atom_list[MAX_ATOMS];
static double *     atom_basis = 0;
static size_t       n_atoms = 0;
static size_t       n_basis = 0;
static size_t       n_coeff = 0;
static double *     norm_l = 0;
static r_t *        R_list = 0;
static basis_set_t  basis_set;
static char         geom_file[MAX_GEOM_FILE_SIZE];
static const char * basis_tag = "basis\n";
static const char * geom_tag = "geometry\n";

/*
/
/  A geometry file will have this form:
/
/    basis
/    3.0 2.0 1.0 0.4
/
/    geometry
/    h 0 0 0
/    h 0.75 0 0
/
/  From the basis heading, construct the basis map of exponents
/  per given element number Z.
/
/  From the geometry heading, construct the atom list.
/
*/

bool
init_geom_basis(const char * file)
{
    if (!load_geom_file(file)) { return false; }
    if (!init_basis_set(file)) { return false; }
    if (!init_atom_list(file)) { return false; }

    if (n_coeff == 0) {
        fprintf(stderr, "error: n_coeff == 0, invalid basis\n");
        return false;
    }

    n_basis = n_atoms * n_coeff;
    atom_basis = safer_calloc(n_basis, sizeof(double), 0);

    for (size_t a = 0; a < n_atoms; ++a) {
        for (size_t i = 0; i < n_coeff; ++i) {
            atom_basis[a * n_coeff + i] = basis_set.exp[i];
        }
    }

    norm_l = safer_calloc(n_basis, sizeof(double), 0);

    init_R_list();

    dump_basis_set();
    dump_atom_list();

    return true;
}




static bool
load_geom_file(const char * file)
{
    if (!file || !*file) {
        fprintf(stderr, "error: prepare_basis_set: no file\n");
        return false;
    }

    struct stat st;

    if (stat(file, &st) < 0) {
        int errsv = errno;
        fprintf(stderr, "error: stat %s failed: %s\n", file,
                strerror(errsv));
        return false;
    }

    if (st.st_size == 0) {
        fprintf(stderr, "error: file %s is empty\n", file);
        return false;
    }

    if ((size_t) st.st_size > sizeof(geom_file)) {
        fprintf(stderr, "error: file %s is too large\n", file);
        return false;
    }

    int fd = open(file, O_RDONLY);

    if (fd <= 0) {
        int errsv = errno;
        fprintf(stderr, "error: open %s failed: %s\n", file,
                strerror(errsv));
        return false;
    }

    ssize_t n_r = read(fd, geom_file, st.st_size);

    if (n_r <= 0 || n_r != (ssize_t) st.st_size) {
        int errsv = errno;
        fprintf(stderr, "error: read %s failed: %s\n", file,
                strerror(errsv));
        return false;
    }

    geom_file[st.st_size] = '\0';

    fprintf(stderr, "loaded geom file:\n%s\n", geom_file);

    close(fd);

    return true;
}



static bool
init_basis_set(const char * file)
{
    memset(basis_set.exp, 0, sizeof(basis_set.exp));

    const char * p = geom_file;
    const char * start = geom_file;
    const char * end = strstr(geom_file, geom_tag);

    if (memcmp(geom_file, basis_tag, strlen(basis_tag)) != 0) {
        fprintf(stderr, "error: invalid geometry file %s\n", file);
        return false;
    }

    if (!end) {
        fprintf(stderr, "error: geometry not listed in file %s\n", file);
        return false;
    }

    p += strlen(basis_tag);

    while (*p != '\0') {
        if ((p - start) == (end - start)) {
            // end of basis listing
            break;
        }
        else if (isdigit(*p)) {
            if (load_shell(p)) {
                p = skip_line(p);
            }
            else {
                fprintf(stderr, "error: invalid basis listing: %s\n", p);
                return false;
            }
        }
        else {
            p = skip_line(p);
        }
    }

    return true;
}



static bool
get_next_token(double *       v_p,
               const char * * p_p)
{
    // Begins at first float:
    //   <float> <float> <float>
    const char * p = *p_p;
    const char * q = p;
    bool         done = false;

    q = skip_float(p);

    if (*q == '\n') { ++q; done = true; }

    *v_p = atof(p);
    *p_p = q;

    return done;
}



static bool
init_atom_list(const char * file)
{
    const char * p = strstr(geom_file, geom_tag);
    size_t       num_atoms = 0;

    if (!p) {
        fprintf(stderr, "error: geometry not listed in file %s\n", file);
        return false;
    }

    p += strlen(geom_tag);

    while (*p) {
        if (num_atoms >= MAX_ATOMS) {
            fprintf(stderr, "error: too many atoms in geom file %s\n", file);
            return false;
        }

        atom_t * new_atom = &atom_list[num_atoms];

        ++num_atoms;

        if (memcmp(p, "g ", strlen("g ")) == 0) {
            new_atom->z = 0;
            p += 2;
        } else if (memcmp(p, "h ", strlen("g ")) == 0) {
            new_atom->z = 1;
            p += 2;
        } else if (memcmp(p, "he ", strlen("he ")) == 0) {
            new_atom->z = 2;
            p += 3;
        }
        else {
            fprintf(stderr, "error: unsupported element: %s\n", p);
            return false;
        }

        get_next_token(&new_atom->R.x, &p);
        get_next_token(&new_atom->R.y, &p);
        get_next_token(&new_atom->R.z, &p);
    }

    n_atoms = num_atoms;

    return true;
}



static bool
load_shell(const char * p)
{
    // Line has form
    //   3.0 2.0 0.5 0.2
    //
    // Pointer p begins at initial float

    bool done = false;
    int  j = 0;

    while (!done) {
        double exp = 0;
        done = get_next_token(&exp, &p);

        if (exp <= 0) {
            fprintf(stderr, "error: invalid exponent %s\n", p);
            return false;
        }
        else {
            basis_set.exp[j] = exp;
            ++n_coeff;
            ++j;

            if (j == MAX_BASIS_PER_ATOM) {
                fprintf(stderr, "error: too many exponents in basis\n");
                fprintf(stderr, "       Max basis allowed is: %d\n", MAX_BASIS_PER_ATOM);
                fprintf(stderr, "       You have:             %d\n", j);
                return false;
            }
        }
    }

    return true;
}



static void
dump_basis_set(void)
{
    printf("dump_basis_set\n");

    for (size_t k = 0; k < MAX_BASIS_PER_ATOM; ++k) {
        if (basis_set.exp[k] == 0) {
            break;
        }

        printf("%f ", basis_set.exp[k]);
    }

    printf("\n");

    return;
}



static void
dump_atom_list()
{
    for (size_t i = 0; i < n_atoms; ++i) {
        atom_t * atom = &atom_list[i];
        printf("atom %zu Z=%zu\n", i, atom->z);
        printf("  x=%f\n", atom->R.x);
        printf("  y=%f\n", atom->R.y);
        printf("  z=%f\n", atom->R.z);
    }

    return;
}



const char *
skip_line(const char * p)
{
    while (*p != '\0') {
        if (*p++ == '\n') { break; }
    }

    return p;
}



const char *
skip_space(const char * p)
{
    while (*p == ' ') {
        ++p;
    }

    return p;
}



const char *
skip_float(const char * p)
{
    while (isdigit(*p) || *p == '.' || *p == '-') {
        ++p;
    }

    return skip_space(p);
}



static void
init_R_list(void)
{

    /*
    *   R_list is nbasis in length. E.g. if 8 in length, if
    *   R_list was static it would look like:
    *
    *     static r_t  R_list[8] = {{0, 0, 0},
    *                              {0, 0, 0},
    *                              {0, 0, 0},
    *                              {0, 0, 0},
    *                              {0.7344, 0, 0},
    *                              {0.7344, 0, 0},
    *                              {0.7344, 0, 0},
    *                              {0.7344, 0, 0}};
    */

    R_list = safer_calloc(n_basis, sizeof(r_t), 0);

    for (size_t i = 0; i < n_atoms; ++i) {
        for (size_t j = 0; j < n_coeff; ++j) {
            R_list[i * j] = atom_list[i].R;
        }
    }

    return;
}



void
build_overlap(double * S)
{
    double b_sum;
    double b_pro;
    double pref;
    double r_sq;
    double d_exp;

    //
    // Construct (not-normalized) overlap integrals:
    //
    //   (A|B) = [Pi/(a + b)]^(3/2) * exp[-ab * (|Ra - Rb|^2) / (a + b)]
    //

    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            b_sum = atom_basis[i] + atom_basis[j];
            b_pro = atom_basis[i] * atom_basis[j];

            pref = pow((M_PI / b_sum), 1.5);

            r_sq = get_r_diff_sq(R_list[i], R_list[j]);

            d_exp = exp(-b_pro * r_sq / b_sum);

            S[i * n_basis + j] = pref * d_exp;
        }
    }

    //
    // Normalize them:
    //
    //   (a A | a A) = 1
    //   a = sqrt(1 / (A|A))
    //

    for (size_t i = 0; i < n_basis; ++i) {
        norm_l[i] = 1.0 / sqrt(S[i * n_basis + i]);
    }

    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            S[i * n_basis + j] *= (norm_l[i] * norm_l[j]);
        }
    }

    return;
}



void
build_core_hamiltonian(double * H)
{
    // Build core Hamiltonian.
    double * T = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * Z = safer_calloc(n_basis * n_basis, sizeof(double), 0);

    double b_sum;
    double b_prod;
    double r_sq;
    double t;
    double t_1;
    double t_2;
    double t_3;
    double t_4;
    double Zc  = 1;
    double fo;

    //
    // Build normalized kinetic energy integrals:
    //
    //   (A|t|B) = ab/(a + b)[3 - r_sq * 2ab/(a + b)] *
    //             [Pi/(a + b)]^(3/2) *
    //             exp[-r_sq * ab/(a + b)]
    //

    for (size_t mu = 0; mu < n_basis; ++mu) {
        for (size_t nu = 0; nu < n_basis; ++nu) {
            b_sum = atom_basis[mu] + atom_basis[nu];
            b_prod = atom_basis[mu] * atom_basis[nu];

            r_sq = get_r_diff_sq(R_list[mu], R_list[nu]);

            t_1 = b_prod / b_sum;
            t_2 = 3 - r_sq * (2 * b_prod / b_sum);
            t_3 = pow((M_PI / b_sum), 1.5);
            t_4 = exp(-b_prod * r_sq / b_sum);

            T[mu * n_basis + nu] = norm_l[mu] * norm_l[nu] * t_1
                                   * t_2 * t_3 * t_4;
        }
    }

    //
    // Build normalized nuclear energy integrals:
    //
    //   (A|Vz|B) = -2Pi * Zc /(a + b) *
    //              exp[-r_sq_ab * ab/(a + b)] *
    //              F[r_sq_pc * (a + b)]
    //
    // where
    //
    //   F[t] = (1/2) * sqrt(Pi/t) * erf[sqrt(t)]
    //

    for (size_t mu = 0; mu < n_basis; ++mu) {
        for (size_t nu = 0; nu < mu + 1; ++nu) {
            for (size_t c = 0; c < n_atoms; ++c) {
                b_sum = atom_basis[mu] + atom_basis[nu];
                b_prod = atom_basis[mu] * atom_basis[nu];

                r_sq = get_r_diff_sq(R_list[mu], R_list[nu]);

                t_1 = - 2 * M_PI * Zc / b_sum;

                t_2 = exp(-r_sq * b_prod / b_sum);

                if ((mu / n_coeff == nu / n_coeff) && (mu / n_coeff == c)) {
                    // mu, nu, c all same atom center. r diff will be zero.
                    fo = 1;
                }
                else {
                    r_t Rp = get_Rp(mu, nu);

                    t = b_sum * get_r_diff_sq(Rp, R_list[n_coeff * c]);
                    fo = 0.5 * sqrt(M_PI / t) * erf(sqrt(t));
                }

                Z[mu * n_basis + nu] += (norm_l[mu] * norm_l[nu] * t_1
                                        * t_2 * fo);
            }

            Z[nu * n_basis + mu] = Z[mu * n_basis + nu];
        }
    }

    sum_mat(H, T, Z, n_basis);

    free(T);
    free(Z);

    return;
}



double
two_elec_int(const size_t a,
             const size_t b,
             const size_t c,
             const size_t d)
{
    //
    // Build normalized two electron integral:
    //
    //   (AB|CD) = 2 Pi^(5/2) *
    //             1 / [(a + b)(c + d)sqrt(a + b + c + d)] *
    //             exp[-r_sq_ab * ab / (a + b) - r_sq_cd * cd / (c + d)] *
    //             F[r_sq_pq * (a + b)(c + d) / (a + b + c + d)]
    //
    // where
    //
    //   F[t] = (1/2) * sqrt(Pi/t) * erf[sqrt(t)]
    //

    double result = 0;

    double ab_sum = atom_basis[a] + atom_basis[b];
    double cd_sum = atom_basis[c] + atom_basis[d];
    double abcd_sum = ab_sum + cd_sum;
    double ab_prod = atom_basis[a] * atom_basis[b];
    double cd_prod = atom_basis[c] * atom_basis[d];

    double r_sq_ab = get_r_diff_sq(R_list[a], R_list[b]);
    double r_sq_cd = get_r_diff_sq(R_list[c], R_list[d]);

    r_t    R_p = get_Rp(a, b);
    r_t    R_q = get_Rp(c, d);

    double r_sq_pq = get_r_diff_sq(R_p, R_q);

    double t_1 = 2 * pow(M_PI, 2.5);
    double t_2 = ab_sum * cd_sum * sqrt(abcd_sum);
    double t_3 = -r_sq_ab * ab_prod / ab_sum -r_sq_cd * cd_prod / cd_sum;

    double t = r_sq_pq * ab_sum * cd_sum / abcd_sum;

    double fo = 1;

    if (t > 0.0001) {
        fo = 0.5 * sqrt(M_PI / t) * erf(sqrt(t));
    }

    result = norm_l[a] * norm_l[b] * norm_l[c] * norm_l[d] *
             t_1 * exp(t_3) * fo / t_2;

    return result;
}



static double
get_r_diff_sq(const r_t a,
              const r_t b)
{
    double x = a.x - b.x;
    double y = a.y - b.y;
    double z = a.z - b.z;

    return ((x * x) + (y * y) + (z * z));
}



static r_t
get_Rp(const size_t mu,
       const size_t nu)
{
    //
    // Rp = (a * R_a + b * R_b)/(a + b)
    //

    r_t    Rp;

    if (mu == nu) {
        Rp.x = R_list[mu].x;
        Rp.y = R_list[mu].y;
        Rp.z = R_list[mu].z;

        return Rp;
    }

    double b_mu = atom_basis[mu];
    double b_nu = atom_basis[nu];
    double b_sum = b_mu + b_nu;

    Rp.x = ((b_mu * R_list[mu].x) + (b_nu * R_list[nu].x)) / b_sum;
    Rp.y = ((b_mu * R_list[mu].y) + (b_nu * R_list[nu].y)) / b_sum;
    Rp.z = ((b_mu * R_list[mu].z) + (b_nu * R_list[nu].z)) / b_sum;

    return Rp;
}



size_t
get_n_basis(void)
{
    return n_basis;
}
