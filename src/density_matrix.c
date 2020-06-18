#include "spin_types.h"
#include "options.h"
#include "basis.h"
#include "matrix.h"
#include "util.h"

static void build_dens_mat_block(double * P, const double * Cx,
                                 const double * Cz, const size_t n_ele,
                                 const size_t occupancy);
static void build_dens_mat_rhf(spin_matrix_t * P, const spinor_t * C);
static void build_dens_mat_uhf(spin_matrix_t * P, const spinor_t * C);
static void build_dens_mat_ghf(spin_matrix_t * P, const spinor_t * C);
static void build_dens_mat_block(double * P, const double * Cx,
                                 const double * Cy, const size_t   n_ele,
                                 const size_t   occupancy);



void
build_density_matrix(spin_matrix_t *  P,
                     const spinor_t * C)
{
    switch (get_hf_type()) {
    case RHF:
        return build_dens_mat_rhf(P, C);
    case UHF:
        return build_dens_mat_uhf(P, C);
    case GHF:
        return build_dens_mat_ghf(P, C);
    }
}



static void
build_dens_mat_rhf(spin_matrix_t *  P,
                   const spinor_t * C)
{
    size_t occupancy = 2;
    size_t n_ele = P->n_ele;

    build_dens_mat_block(P->uu, C->u, C->u, n_ele, occupancy);

    return;
}



static void
build_dens_mat_uhf(spin_matrix_t *  P,
                   const spinor_t * C)
{
    size_t occupancy = 1;
    size_t n_ele_u = P->n_ele_u;
    size_t n_ele_d = P->n_ele_d;

    build_dens_mat_block(P->uu, C->u, C->u, n_ele_u, occupancy);
    build_dens_mat_block(P->dd, C->d, C->d, n_ele_d, occupancy);

    return;
}



static void
build_dens_mat_ghf(spin_matrix_t *  P,
                   const spinor_t * C)
{
    size_t occupancy = 1;
    size_t n_ele = P->n_ele;

    build_dens_mat_block(P->uu, C->u, C->u, n_ele, occupancy);
    build_dens_mat_block(P->dd, C->d, C->d, n_ele, occupancy);
    build_dens_mat_block(P->ud, C->u, C->d, n_ele, occupancy);
    build_dens_mat_block(P->du, C->d, C->u, n_ele, occupancy);

    return;
}



static void
build_dens_mat_block(double *       P,
                     const double * Cx,
                     const double * Cy,
                     const size_t   n_ele,
                     const size_t   occupancy)
{
    //
    // Given coefficient matrix C, build density matrix
    //
    //   P_ij = 2 Sum[ C_ia * C_ja, a, 0, N_elec / occupancy]
    //

    size_t n_basis = get_n_basis();
    size_t n_mo = n_ele / occupancy;

    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            P[i * n_basis + j] = 0;

            for (size_t a = 0; a <  n_mo ; ++a) {
                P[i * n_basis + j] += occupancy * Cx[i * n_basis + a]
                                                * Cy[j * n_basis + a];
            }
        }
    }

    return;
}
