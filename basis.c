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
#include "log.h"

#define MAX_GEOM_FILE_SIZE 1024

static bool init_basis_map(const char * file);
static bool init_atom_list(const char * file);
static bool load_geom_file(const char * file);
static bool load_shell(size_t Z, const char * * geom_line, int * left);
static void dump_basis_map(void);
static void dump_atom_list(void);

static atom_t       atom_list[MAX_ATOMS];
static R_t *        R_list = 0;
static atom_basis_t basis_map[MAX_Z];
static char         geom_file[MAX_GEOM_FILE_SIZE];
static const char * basis_tag = "basis\n";
static const char * geom_tag = "geometry\n";

/*
/
/  A geometry file will have this form:
/
/    basis
/    h s 3.0 2.0 1.0 0.4
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
    if (!init_basis_map(file)) { return false; }
    if (!init_atom_list(file)) { return false; }

    dump_basis_map();
    dump_atom_list();

    return true;
}



static bool
init_basis_map(const char * file)
{
    memset(basis_map, 0, sizeof(basis_map));

    const char * p = geom_file;
    int          left = strlen(p);

    if (memcmp(geom_file, basis_tag, strlen(basis_tag)) != 0) {
        fprintf(stderr, "error: invalid geometry file %s\n", file);
        return false;
    }

    p += strlen(basis_tag);
    left -= strlen(basis_tag);

    if (left <= 0) {
        fprintf(stderr, "error: invalid geometry file %s\n", file);
        return false;
    }

    switch (*p) {
    case 'h':
    case 'H':
        if (!load_shell(1, &p, &left)) {
            fprintf(stderr, "error: invalid basis listing: %s\n", p);
            return false;
        }
        else {
            break;
        }

    default:
        fprintf(stderr, "error: invalid basis listing: %s\n", p);
        return false;
    }

    return true;
}



static bool
get_next_token(double *       v_p,
               const char * * p_p,
               int *          left_p)
{
    const char * p = *p_p;
    const char * q = p;
    bool         done = false;
    size_t       len = 0;

    while (*q != '\0') {
        if (*q == ' ') { ++q; break; }
        if (*q == '\n') { ++q; done = true; break; }
        ++q;
    }

    len = (size_t) (p - q);
    *v_p = atof(p);

    *p_p = q;
    *left_p -= len;

    return done;
}



static bool
init_atom_list(const char * file)
{
    const char * p = strstr(geom_file, geom_tag);
    int          left = strlen(p);
    size_t       num_atoms = 0;

    if (!p || !left) {
        log_error("geometry not listed in file");
        return false;
    }

    p += strlen(geom_tag);
    left -= strlen(geom_tag);

    while (left && *p) {
        if (num_atoms >= MAX_ATOMS) {
            log_error("too many atoms in geom file");
            return false;
        }

        atom_t * new_atom = &atom_list[num_atoms];

        ++num_atoms;

        if (memcmp(p, "h ", 2) == 0) {
            new_atom->z = 1;
            p += 2;
            left -= 2;
        } else if (memcmp(p, "he ", 3) == 0) {
            new_atom->z = 2;
            p += 3;
            left -= 3;
        }
        else {
            fprintf(stderr, "error: unsupported element: %s\n", p);
            return false;
        }

        get_next_token(&new_atom->R.x, &p, &left);
        get_next_token(&new_atom->R.y, &p, &left);
        get_next_token(&new_atom->R.z, &p, &left);
    }

    return true;
}



static bool
load_geom_file(const char * file)
{
    if (!file || !*file) {
        fprintf(stderr, "error: prepare_basis_map: no file\n");
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
        fprintf(stderr, "error: fle %s is empty\n", file);
        return false;
    }

    if ((size_t) st.st_size > sizeof(geom_file)) {
        fprintf(stderr, "error: fle %s is too large\n", file);
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
load_shell(size_t         Z,
           const char * * geom_line,
           int *          left_p)
{
    const char * p = *geom_line;
    int          left = *left_p;
    char         type = '\0';
    size_t       i = 0;

    // Line has form
    //   h s 3.0 2.0 0.5 0.2
    //
    // Pointer p begins at initial "h s ..."

    if (left < 2) { goto error; }
    p += 2;
    left -= 2;

    type = *p;

    if (type != 's' && type != 'p') {
        fprintf(stderr, "error: invalid shell type: %c\n", type);
        return false;
    }

    if (left < 2) { goto error; }
    p += 2;
    left -= 2;

    while (basis_map[Z].shells[i].exp[0] > 0) { 
        ++i;
        if (i == MAX_SHELLS) { goto error; }
    }

    struct atom_shell_t * new_shell = &basis_map[Z].shells[i];

    new_shell->type = type;

    bool done = false;
    int  j = 0;

    while (!done) {
        double exp = 0;
        done = get_next_token(&exp, &p, &left);

        if (exp == 0) {
            fprintf(stderr, "error: invalid exponent %s\n", p);
            return false;
        }

        new_shell->exp[j] = exp;

        ++j;
        if (j == MAX_BASIS_PER_ATOM) { goto error; }
    }

    *geom_line = p;
    *left_p = left;

    return true;

error:
    fprintf(stderr, "error: invalid basis listing: %s\n", p);
    return false;
}



static void
dump_basis_map(void)
{
    printf("dump_basis_map\n");

    atom_basis_t * ab;

    for (size_t i = 0; i < MAX_Z; ++i) {
        ab = &basis_map[i];

        for (size_t j = 0; j < MAX_SHELLS; ++j) {
            printf("basis %c\n", ab->shells[j].type);

            for (size_t k = 0; k < MAX_BASIS_PER_ATOM; ++k) {
                if (ab->shells[j].exp[k] == 0) {
                    break;
                }

                printf("%f ", ab->shells[j].exp[k]);
            }

            printf("\n");
        }
    }

    return;
}



static void
dump_atom_list()
{

    return;
}
