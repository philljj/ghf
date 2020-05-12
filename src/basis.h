#if !defined(BASIS_H)
#define BASIS_H

#include <stdbool.h>

#define MAX_ATOMS          1024
#define MAX_BASIS_PER_ATOM 16
#define MAX_Z              2
#define MAX_SHELLS         4

struct basis_set_t {
    double exp[MAX_BASIS_PER_ATOM]; // List of gaussian exponents.
};

typedef struct basis_set_t basis_set_t;

typedef struct {
    double x;
    double y;
    double z;
} r_t;

typedef struct {
    size_t Z;
    r_t    R;
} atom_t;

bool init_geom_basis(const char * file, bool debug);
void build_overlap(double * S);
void build_core_hamiltonian(double * H);
double two_elec_int(const size_t a, const size_t b,
                    const size_t c, const size_t d);
void   build_density_matrix(double * P, const double * C);
double nuclear_rep_energy(void);
size_t get_n_basis(void);
bool   is_debug(void);
#endif /* if !defined(BASIS_H) */
