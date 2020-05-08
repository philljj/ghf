#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "basis.h"
#include "matrix.h"
#include "util.h"

extern void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
                   double* w, double* work, int* lwork, int* info );



void
sum_mat(double *       C,
        const double * A,
        const double * B,
        const size_t   len)
{
    //
    // C = A + B
    //

    for (size_t i = 0; i < len * len; ++i) {
        C[i] = A[i] + B[i];
    }

    return;
}



void
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



void
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
    //  TODO: refactor this so it's n_basis * n_basis in scaling.

    size_t   n_basis = get_n_basis();
    double * s_exp = safer_calloc(n_basis, sizeof(double), 0);

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



void
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

    size_t   n_basis = get_n_basis();
    double * eig_val = 0;
    double * work;
    double   wkopt;
    int      info;
    int      lda = n_basis;
    int      lwork = -1;
    int      n = n_basis;

    if (c) {
        eig_val = c;
    }
    else {
        eig_val = safer_calloc(n_basis, sizeof(double), 0);
    }

   /*
   *  dsyev.f comments:
   *    If LWORK = -1, then a workspace query is assumed; the routine
   *    only calculates the optimal size of the WORK array, returns
   *    this value as the first entry of the WORK array, and no error
   *    message related to LWORK is issued by XERBLA.
   */

    dsyev_("V", "U", &n, C, &lda, eig_val, &wkopt, &lwork, &info);

    lwork = (int) wkopt;

    if (lwork <= 0) {
        fprintf(stderr, "error: dsyev: lwork = %d\n", lwork);
        exit(EXIT_FAILURE);
    }

    work = safer_calloc(lwork, sizeof(double), 0);

    dsyev_("V", "U", &n, C, &lda, eig_val, work, &lwork, &info);

    free(work);

    if (info == 0) {
        // Success. Debug message?
    }
    else if (info < 0) {
        fprintf(stderr, "error: dsyev: illegal value: info = %d\n", info);
        exit(EXIT_FAILURE);
    }
    else if (info > 0) {
        fprintf(stderr, "error: dsyev: failed to converge: info = %d\n", info);
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
