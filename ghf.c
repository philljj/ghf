#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#define E_CONV   0.000000001
#define MAX_ITER 20

typedef struct {
    double x;
    double y;
    double z;
} r_t;

typedef struct {
    size_t z;
    r_t    R;
} atom_t;

static void   print_usage_and_die(void) __attribute__((__noreturn__));
static void   init_basis(void);
static double get_r_diff_sq(const r_t a, const r_t b);
static r_t    get_Rp(const size_t mu, const size_t nu);
static void   build_overlap_matrices(double * S, double * X, double * Y);
static void   build_overlap(double * S);
static void   build_density_matrix(double * P, const double * C);
static void   scf_procedure(const double * S, const double * X,
                            const double * Y, const double * H);
static void   population_analysis(const double * P, const double * S,
                                  const double * Y);
static void   build_core_hamiltonian(double * H);
static double two_elec_int(const size_t a, const size_t b,
                                      const size_t c, const size_t d);
static void   build_fock(double * F, const double * H, const double * P);
static void   diag_fock(double * C, const double * F, const double * X);
static void   call_dsyev(double * C, double * c);
static void   sum_mat(double * c, const double * a, const double * b,
                      const size_t len);
static void   mult_mat(double * C, const double * A, const double * B,
                       const size_t len);
static void   build_spectral_mat(double * C, const double * U,
                                 const double * s, const double exp);
extern void   dsyev(char* jobz, char* uplo, int* n, double* a, int* lda,
                    double* w, double* work, int* lwork, int* info );
static void   print_matrix(const double * M, const size_t len,
                           const char * what);

static double hartree_fock_energy(const double * P, const double * H,
                                  const double * F);

static size_t   n_atoms = 2;
static size_t   n_ele = 2;
static double   b_coeff[] = {3.0, 2.0, 1.0, 0.4};
static size_t   n_coeff;
static double * basis;
static size_t   debug = 0;

static r_t      R_list[8] = {{0, 0, 0},
                             {0, 0, 0},
                             {0, 0, 0},
                             {0, 0, 0},
                             {0.7344, 0, 0},
                             {0.7344, 0, 0},
                             {0.7344, 0, 0},
                             {0.7344, 0, 0}};

static size_t   n_basis = 0;
static double * norm_l = 0;



int
main(int    argc,
     char * argv[])
{
    if (argc < 2) {
        print_usage_and_die();
    }

    int    opt = 0;
    size_t threads = 0;

    while ((opt = getopt(argc, argv, "di:t:?")) != -1) {
        switch (opt) {
        case 'd':
            debug = 1;
            break;
        case 'i':
            printf("using input file: %s\n", optarg);
            break;

        case 't':
            threads = strtoul(optarg, 0, 10);
            printf("using threads: %zu\n", threads);
            break;

        case '?':
        default:
            print_usage_and_die();
        }
    }

    init_basis();

    double * S = calloc(n_basis * n_basis, sizeof(double));
    double * X = calloc(n_basis * n_basis, sizeof(double));
    double * Y = calloc(n_basis * n_basis, sizeof(double));
    double * H = calloc(n_basis * n_basis, sizeof(double));

    norm_l = calloc(n_basis, sizeof(double));

    build_overlap_matrices(S, X, Y);

    build_core_hamiltonian(H);
    if (debug) { print_matrix(H, n_basis, "Core Hamiltonian"); }

    scf_procedure(S, X, Y, H);

    return EXIT_SUCCESS;
}



static void
scf_procedure(const double * S,
              const double * X,
              const double * Y,
              const double * H)
{
    double   old_energy = 0;
    double   new_energy = 0;
    double * F = calloc(n_basis * n_basis, sizeof(double));
    double * C = calloc(n_basis * n_basis, sizeof(double));
    double * P = calloc(n_basis * n_basis, sizeof(double));

    for (size_t z = 0; z < MAX_ITER; ++z) {
        build_fock(F, H, P);
        diag_fock(C, F, X);
        build_density_matrix(P, C);
        population_analysis(P, S, Y);

        old_energy = new_energy;
        new_energy = hartree_fock_energy(P, H, F);

        double diff_energy = fabs(new_energy - old_energy);

        if (diff_energy < E_CONV) {
            printf("SCF converged at %zu iterations.\n", z);
            printf("Final GHF energy: %.9f\n", new_energy);
            printf("  dE:          %.9f\n", diff_energy);
            break;
        }
        else {
            printf("Iteration: %zu.\n", z);
            printf("  GHF energy: %.9f\n", new_energy);
            printf("  dE:          %.9f\n", diff_energy);
            sleep(1);
        }
    }

    return;
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

    double b_mu = basis[mu];
    double b_nu = basis[nu];
    double b_sum = b_mu + b_nu;

    Rp.x = ((b_mu * R_list[mu].x) + (b_nu * R_list[nu].x)) / b_sum;
    Rp.y = ((b_mu * R_list[mu].y) + (b_nu * R_list[nu].y)) / b_sum;
    Rp.z = ((b_mu * R_list[mu].z) + (b_nu * R_list[nu].z)) / b_sum;

    return Rp;
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
            b_sum = basis[i] + basis[j];
            b_pro = basis[i] * basis[j];

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



static void
build_overlap_matrices(double * S,
                       double * X,
                       double * Y)
{
    //
    // Build overlap matrix S, and then:
    //
    //   X = S^(-1/2)
    //   Y = S^(1/2)
    //

    build_overlap(S);

    if (debug) { print_matrix(S, n_basis, "S"); }

    double * U = calloc(n_basis * n_basis, sizeof(double));
    double * eig_val = calloc(n_basis, sizeof(double));

    memcpy(U, S, n_basis * n_basis * sizeof(double));

    call_dsyev(U, eig_val);

    build_spectral_mat(X, U, eig_val, -0.5);

    if (debug) { print_matrix(X, n_basis, "X"); }

    build_spectral_mat(Y, U, eig_val, 0.5);

    if (debug) { print_matrix(Y, n_basis, "Y"); }

    free(U);
    free(eig_val);

    return;
}



static void
build_density_matrix(double *       P,
                     const double * C)
{
    //
    // Given coefficient matrix C, build density matrix
    //
    //   P_ij = 2 Sum[ C_ia * C_ja, a, 0, N_elec / 2]
    //

    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            P[i * n_basis + j] = 0;

            for (size_t a = 0; a <  n_ele / 2; ++a) {
                P[i * n_basis + j] += 2 * C[i * n_basis + a] *
                                          C[j * n_basis + a];
            }
        }
    }

    if (debug) { print_matrix(P, n_basis, "P"); }

    return;
}



static void
population_analysis(const double * P,
                    const double * S,
                    const double * Y)
{
    double * T = calloc(n_basis * n_basis, sizeof(double));
    double * V = calloc(n_basis * n_basis, sizeof(double));

    mult_mat(T, P, S, n_basis);

    double n = 0;

    for (size_t i = 0; i < n_basis; ++i) {
        n += T[i * n_basis + i];
    }

    printf("n=%f\n", n);

    mult_mat(T, Y, P, n_basis);
    mult_mat(V, T, Y, n_basis);

    if (debug) { print_matrix(V, n_basis, "P in Lowdin basis"); }

    return;
}



void
build_core_hamiltonian(double * H)
{
    // Build core Hamiltonian.
    double * T = calloc(n_basis * n_basis, sizeof(double));
    double * Z = calloc(n_basis * n_basis, sizeof(double));

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
            b_sum = basis[mu] + basis[nu];
            b_prod = basis[mu] * basis[nu];

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
                b_sum = basis[mu] + basis[nu];
                b_prod = basis[mu] * basis[nu];

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



static double
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

    double ab_sum = basis[a] + basis[b];
    double cd_sum = basis[c] + basis[d];
    double abcd_sum = ab_sum + cd_sum;
    double ab_prod = basis[a] * basis[b];
    double cd_prod = basis[c] * basis[d];

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



static void
call_dsyev(double * C,
           double * c)
{
    //
    // Pass in matrix C to be diagonalized, and optionally
    // vector c to receive the eigenvalues. If c is null,
    // scratch space for eigenvalues is allocated, and eigvals
    // are not stored.
    //
    // Eigenvectors are stored in C.
    //

    double *   eig_val = 0;
    double *   work;
    double     wkopt;
    int        info;
    int        lda = n_basis;
    int        lwork = -1;
    int        n = n_basis;

    if (c) {
        eig_val = c;
    }
    else {
        eig_val = calloc(n_basis, sizeof(double));
    }

    dsyev("Vectors", "Upper", &n, C, &lda, eig_val, &wkopt, &lwork, &info);

    lwork = (int) wkopt;
    work = calloc(lwork, sizeof(double));

    dsyev("Vectors", "Upper", &n, C, &lda, eig_val, work, &lwork, &info);

    free(work);

    if (info) {
        fprintf(stderr, "error: dsyev failed to compute eigenvalues\n");
        exit(EXIT_FAILURE);
    }

    if (!c) {
        printf("eigenvalues\n");
        for (size_t i = 0; i < n_basis; ++i) {
            printf( " %6.16f", eig_val[i]);
        }
        printf("\n");

        free(eig_val);
    }

    return;
}



static void
build_fock(double *       F,
           const double * H,
           const double * P)
{

    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            F[i * n_basis + j] = H[i * n_basis + j];

            for (size_t k = 0; k < n_basis; ++k) {
                for (size_t l = 0; l < n_basis; ++l ) {
                    F[i * n_basis + j] += P[k * n_basis + l] *
                                          two_elec_int(i, j, k, l);

                    F[i * n_basis + j] -= 0.5 *
                                          P[k * n_basis + l] *
                                          two_elec_int(i, k, j, l);

                }
            }
        }
    }

    if (debug) { print_matrix(F, n_basis, "Fock matrix"); }

    return;
}



static void
diag_fock(double *       C,
          const double * F,
          const double * X)
{
    //
    // Given Fock matrix F, build
    //
    //   F` = X^t F X
    //
    // then solve
    //
    //   F`C` = C`e
    //
    //   C = X C`
    //

    double * Fx = calloc(n_basis * n_basis, sizeof(double));
    double * T  = calloc(n_basis * n_basis, sizeof(double));

    mult_mat(T, X, F, n_basis);
    mult_mat(Fx, T, X, n_basis);

    call_dsyev(Fx, 0);

    mult_mat(C, X, Fx, n_basis);

    if (debug) { print_matrix(C, n_basis, "C"); }


    free(Fx);
    free(T);

    return;
}



static void
sum_mat(double *       C,
        const double * A,
        const double * B,
        const size_t   len)
{
    //
    // C = A + B
    //

    for (size_t i = 0; i < len; ++i) {
        for (size_t j = 0; j < len; ++j) {
            C[i * len + j] = A[i * len + j] +
                             B[i * len + j];
        }
    }

    return;
}



static void
mult_mat(double *       C,
         const double * A,
         const double * B,
         const size_t   len)
{
    //
    // C = Transpose[A] * B
    //

    for (size_t i = 0; i < len; ++i) {
        for (size_t j = 0; j < len; ++j) {
            C[i * len + j] = 0;

            for (size_t k = 0; k < len; ++k) {
                C[i * len + j] += A[i * len + k] *
                                  B[j * len + k];
            }
        }
    }

    return;
}



static void
build_spectral_mat(double *       C,
                   const double * U,
                   const double * s,
                   const double   exp)
{
    //
    // Given a vector of eigenvalues s_i, and unitary transformation
    // matrix U, such that:
    //
    //   s = U^t * S * U
    //
    // Build new spectral matrix C, given by:
    //
    //   C = U * s^(exp) * U^t
    //
    //   C_ij = U_ik * [s^(exp)]_k * U_kj
    //

    double * s_exp = calloc(n_basis, sizeof(double));

    for (size_t k = 0; k < n_basis; ++k) {
        s_exp[k] = pow(s[k], exp);
    }

    for (size_t i = 0; i < n_basis; i++) {
        for (size_t j = 0; j < n_basis; ++j) {
            C[i * n_basis + j] = 0;

            for (size_t k = 0; k < n_basis; ++k) {
                C[i * n_basis + j] += U[k * n_basis + i] *
                                      s_exp[k] *
                                      U[k * n_basis + j];
            }
        }
    }

    free(s_exp);
    return;
}



static void
init_basis(void)
{
    n_coeff = sizeof(b_coeff) / sizeof(b_coeff[0]);

    n_basis = n_atoms * n_coeff;
    basis = calloc(n_basis, sizeof(double));

    for (size_t a = 0; a < n_atoms; ++a) {
        for (size_t i = 0; i < n_coeff; ++i) {
            basis[a * n_coeff + i] = b_coeff[i];
        }
    }

    return;
}



static void
print_usage_and_die(void)
{
    fprintf(stderr, "usage\n");
    exit(EXIT_FAILURE);
}



void
print_matrix(const double * M,
             const size_t   len,
             const char *   what)
{
    printf("%s\n", what);

    for (size_t mu = 0; mu < len; ++mu) {
        for (size_t nu = 0; nu < len; ++nu) {
            printf(" %6.16f", M[mu * len + nu]);
        }

        printf("\n");
    }

    printf("\n");

    return;
}



static double
hartree_fock_energy(const double * P,
                    const double * H,
                    const double * F)
{
    double energy = 0;

    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            energy += 0.5 * 0.5 * P[i * n_basis + j] *
                      (H[i * n_basis + j] + F[i * n_basis + j]);
        }
    }

    printf("scf energy = %6.16f\n", energy);

    return energy;
}
