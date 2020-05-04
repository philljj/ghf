#if !defined(BASIS_H)
#define BASIS_H

#include <stdbool.h>

#define MAX_ATOMS          1024
#define MAX_BASIS_PER_ATOM 16
#define MAX_Z              2
#define MAX_SHELLS         4

struct atom_shell_t {
    /*
    * 
    */
    char   type;                    // s, p, d, ...
    double exp[MAX_BASIS_PER_ATOM]; // List of gaussian exponents.
};

struct atom_basis_t {
    /*
    * An atom is a bundle of MAX_SHELLS.
    */
    struct atom_shell_t shells[MAX_SHELLS];
};

typedef struct atom_basis_t atom_basis_t;

typedef struct {
    double x;
    double y;
    double z;
} R_t;

typedef struct {
    size_t z;
    R_t    R;
} atom_t;

bool init_geom_basis(const char * file);

#endif /* if !defined(BASIS_H) */
