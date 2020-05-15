#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "basis.h"
#include "matrix.h"
#include "g_matrix.h"
#include "util.h"

static double two_elec_int_in_mem(const size_t a, const size_t b,
                                  const size_t c, const size_t d);

static size_t   n_basis = 0;
static size_t   n_basis_sq = 0;
static size_t   n_basis_cub = 0;
static size_t   n_basis_quad = 0;
static double * v_ao = 0;
static size_t   threads = 1;



void
precalculate_two_elec_int(size_t opt_threads)
{
    if (v_ao) {
        fprintf(stderr, "error: double initialization of v_ao\n");
        exit(EXIT_FAILURE);
    }

    threads = opt_threads;

    n_basis = get_n_basis();

    n_basis_sq   = n_basis * n_basis;
    n_basis_cub  = n_basis * n_basis * n_basis;
    n_basis_quad = n_basis * n_basis * n_basis * n_basis;

    v_ao = safer_calloc(n_basis_quad, sizeof(double), "init_R_list");

    //#pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            for (size_t k = 0; k < n_basis; ++k) {
                for (size_t l = 0; l < n_basis; ++l ) {
                    v_ao[i * n_basis_cub + j * n_basis_sq + k * n_basis + l] = two_elec_int(i, j, k, l);
                }
            }
        }
    }

    return;
}



void
free_two_elec_int(void)
{
    if (v_ao) {
        free(v_ao);
        v_ao = 0;
    }
}



double
two_elec_int_in_mem(const size_t a,
                    const size_t b,
                    const size_t c,
                    const size_t d)
{
    // Might be wasteful to return this through function rather
    // than access directly. Will worry about that when it is a bottleneck.
    return v_ao[a * n_basis_cub + b * n_basis_sq + c * n_basis + d];
}



void
build_G_matrix_in_memory(double *       G,
                         const double * P)
{
    //
    // Build the 2-electron part of the Fock matrix:
    //   Fij = Hij + Gij
    //

    if (n_basis == 0) { n_basis = get_n_basis(); }

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            for (size_t k = 0; k < n_basis; ++k) {
                for (size_t l = 0; l < n_basis; ++l ) {
                    G[i * n_basis + j] += P[k * n_basis + l] *
                                          two_elec_int_in_mem(i, j, k, l);

                    G[i * n_basis + j] -= 0.5 *
                                          P[k * n_basis + l] *
                                          two_elec_int_in_mem(i, k, j, l);

                }
            }
        }
    }

    return;
}



void
build_G_matrix_on_the_fly(double *       G,
                          const double * P)
{
    //
    // Build the 2-electron part of the Fock matrix:
    //   Fij = Hij + Gij
    //

    if (n_basis == 0) { n_basis = get_n_basis(); }

    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            for (size_t k = 0; k < n_basis; ++k) {
                for (size_t l = 0; l < n_basis; ++l ) {
                    G[i * n_basis + j] += P[k * n_basis + l] *
                                          two_elec_int(i, j, k, l);

                    G[i * n_basis + j] -= 0.5 *
                                          P[k * n_basis + l] *
                                          two_elec_int(i, k, j, l);

                }
            }
        }
    }

    return;
}
