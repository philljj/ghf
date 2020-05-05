#if !defined (MATRIX_H)
#define MATRIX_H
void sum_mat(double * c, const double * a, const double * b,
             const size_t len);
void mult_mat(double * C, const double * A, const double * B,
              const size_t len);
void build_spectral_mat(double * C, const double * U,
                        const double * s, const double exp);
void call_dsyev(double * C, double * c);
#endif /* !defined (MATRIX_H) */
