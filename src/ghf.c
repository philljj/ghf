#include <ctype.h>
#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "basis.h"
#include "matrix.h"
#include "util.h"

#define E_CONV   0.000000001
#define MAX_ITER 20

struct fock_matrix_t {
    /*
    * Fock matrix is nbasis * nbasis square matrix.
    *
    * F is the actual Fock matrix array.
    *
    * in_memory indicates if 2-electron integrals are stored
    * or re-calculated on the fly.
    *
    * The make_fock callback is loaded based on the type of
    * SCF calculation to be performed.
    */

    size_t     nbasis;
    double *   F;
    bool       in_memory;
    void     (*make_fock)(struct fock_matrix_t *f, const double *P);
};

static void   print_usage_and_die(void) __attribute__((__noreturn__));
static void   build_overlap_matrices(double * S, double * X, double * Y);
static void   scf_procedure(const double * S, const double * X,
                            const double * Y, const double * H);
static void   population_analysis(const double * P, const double * S,
                                  const double * Y);
static void   build_fock(double * F, const double * H, const double * P);
static void   diag_fock(double * C, const double * F, const double * X);
static double hartree_fock_energy(const double * P, const double * H,
                                  const double * F);

static size_t debug = 0;
static size_t n_basis = 0;
static size_t sleep_time = 0;



int
main(int    argc,
     char * argv[])
{
    if (argc < 2) {
        print_usage_and_die();
    }

    int          opt = 0;
    size_t       threads = 0;
    const char * input_file = 0;

    while ((opt = getopt(argc, argv, "dsf:t:?")) != -1) {
        switch (opt) {
        case 'd':
            debug = 1;
            break;

        case 'f':
            printf("using input file: %s\n", optarg);
            input_file = optarg;
            break;

        case 't':
            threads = strtoul(optarg, 0, 10);
            printf("using threads: %zu\n", threads);
            break;

        case 's':
            sleep_time = 1;
            printf("using sleep: %zu\n", sleep_time);
            break;

        case '?':
        default:
            print_usage_and_die();
        }
    }

    if (!init_geom_basis(input_file)) {
        return EXIT_FAILURE;
    }

    n_basis = get_n_basis();

    if (n_basis == 0) {
        fprintf(stderr, "error: n_basis == 0\n");
        return EXIT_FAILURE;
    }

    double * S = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * X = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * Y = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * H = safer_calloc(n_basis * n_basis, sizeof(double), 0);

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
    double * F = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * C = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * P = safer_calloc(n_basis * n_basis, sizeof(double), 0);

    for (size_t z = 0; z < MAX_ITER; ++z) {
        build_fock(F, H, P);
        diag_fock(C, F, X);
        build_density_matrix(P, C);
        population_analysis(P, S, Y);

        if (debug) { print_matrix(P, n_basis, "P"); }

        old_energy = new_energy;
        new_energy = hartree_fock_energy(P, H, F);

        double diff_energy = fabs(new_energy - old_energy);

        if (diff_energy < E_CONV) {
            printf("SCF converged at %zu iterations.\n", z);
            printf("Final GHF energy: %.9f\n", new_energy);
            printf("  dE:          %.9f\n", diff_energy);
            return;
        }
        else {
            printf("Iteration: %zu.\n", z);
            printf("  GHF energy: %.9f\n", new_energy);
            printf("  dE:          %.9f\n", diff_energy);
        }

        if (sleep_time) { sleep(sleep_time); }
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

    double * U = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * eig_val = safer_calloc(n_basis, sizeof(double), 0);

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
population_analysis(const double * P,
                    const double * S,
                    const double * Y)
{
    double * T = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * V = safer_calloc(n_basis * n_basis, sizeof(double), 0);

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

    double * Fx = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * T  = safer_calloc(n_basis * n_basis, sizeof(double), 0);

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
print_usage_and_die(void)
{
    fprintf(stderr, "usage:\n");
    fprintf(stderr, "  ghf -f <path to geom file> [-ds]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "options:\n");
    fprintf(stderr, "  -d    enables debug\n");
    fprintf(stderr, "  -s    sleep for 1s between iterations\n");
    exit(EXIT_FAILURE);
}



static double
hartree_fock_energy(const double * P,
                    const double * H,
                    const double * F)
{
    double energy = 0;

    for (size_t i = 0; i < n_basis; ++i) {
        for (size_t j = 0; j < n_basis; ++j) {
            energy +=  0.5 * P[i * n_basis + j] *
                      (H[i * n_basis + j] + F[i * n_basis + j]);
        }
    }

    printf("scf energy = %6.16f\n", energy);

    return energy;
}

