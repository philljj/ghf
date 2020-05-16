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
#include "g_matrix.h"
#include "options.h"
#include "util.h"

#define E_CONV          0.000000001
#define MAX_ITER        20

static void   print_usage_and_die(void) __attribute__((__noreturn__));
static void   build_overlap_matrices(double * S, double * X, double * Y);
static void   scf_procedure(const double * S, const double * X,
                            const double * Y, const double * H);
static void   population_analysis(const double * P, const double * S,
                                  const double * Y);
static void   build_fock(double * F, const double * H, const double * P);
static void (*build_G_matrix)(double *G, const double *P);
static void   diag_fock(double * C, const double * F, const double * X);
static double hartree_fock_energy(const double * P, const double * H,
                                  const double * F);

static size_t n_basis = 0;



int
main(int    argc,
     char * argv[])
{
    if (argc < 2) {
        print_usage_and_die();
    }

    if (!validate_options(argc, argv)) {
        print_usage_and_die();
    }

    if (!init_geom_basis(get_input_file(), is_debug())) {
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

    if (in_memory()) {
        printf("info: storing two-electron integrals in memory\n");
        build_G_matrix = build_G_matrix_in_memory;
        precalculate_two_elec_int();
    }
    else {
        build_G_matrix = build_G_matrix_on_the_fly;
    }

    scf_procedure(S, X, Y, H);

    free(H);
    free(Y);
    free(X);
    free(S);

    if (in_memory()) {
        free_two_elec_int();
    }

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
    double   nuc_energy = nuclear_rep_energy();

    double * F = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * C = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * P = safer_calloc(n_basis * n_basis, sizeof(double), 0);

    for (size_t z = 0; z < MAX_ITER; ++z) {
        build_fock(F, H, P);
        diag_fock(C, F, X);
        build_density_matrix(P, C);
        population_analysis(P, S, Y);

        if (is_debug()) { print_matrix(P, n_basis, "P"); }

        old_energy = new_energy;
        new_energy = hartree_fock_energy(P, H, F);

        double diff_energy = fabs(new_energy - old_energy);

        if (diff_energy < E_CONV) {
            printf("SCF converged at %zu iterations.\n", z);
            printf("Final GHF energy: %.9f\n", new_energy);
            printf(" energy per atom: %.9f\n", new_energy / get_n_atoms());
            printf("  nuc rep energy: %.9f\n", nuc_energy);
            printf("              dE: %.9f\n", diff_energy);
            printf("\n");
            printf("Total GHF energy: %.9f\n", new_energy + nuc_energy);
            return;
        }
        else if (!is_quiet()) {
            printf("Iteration: %zu.\n", z);
            printf("  electronic energy: %.9f\n", new_energy);
            printf("     nuc rep energy: %.9f\n", nuc_energy);
            printf("                 dE: %.9f\n", diff_energy);
        }

        if (get_sleep_time()) { sleep(get_sleep_time()); }
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

    if (is_debug()) { print_matrix(S, n_basis, "S"); }

    double * U = safer_calloc(n_basis * n_basis, sizeof(double), 0);
    double * eig_val = safer_calloc(n_basis, sizeof(double), 0);

    memcpy(U, S, n_basis * n_basis * sizeof(double));

    call_dsyev(U, eig_val);

    build_spectral_mat(X, U, eig_val, -0.5);

    if (is_debug()) { print_matrix(X, n_basis, "X"); }

    build_spectral_mat(Y, U, eig_val, 0.5);

    if (is_debug()) { print_matrix(Y, n_basis, "Y"); }

    free(U);
    free(eig_val);

    return;
}




static void
population_analysis(const double * P,
                    const double * S,
                    const double * Y)
{
    if (is_quiet()) { return; }

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

    if (is_debug()) { print_matrix(V, n_basis, "P in Lowdin basis"); }

    return;
}



static void
build_fock(double *       F,
           const double * H,
           const double * P)
{
    double * G = safer_calloc(n_basis * n_basis, sizeof(double), "G matrix");

    build_G_matrix(G, P);
    sum_mat(F, H, G, n_basis);

    if (is_debug()) { print_matrix(H, n_basis, "H"); }
    if (is_debug()) { print_matrix(G, n_basis, "G"); }
    if (is_debug()) { print_matrix(F, n_basis, "Fock matrix"); }

    free(G);

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

    if (is_debug()) { print_matrix(C, n_basis, "C"); }

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

    return energy;
}
