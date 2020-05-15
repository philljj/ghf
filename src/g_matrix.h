#if !defined(G_MATRIX_H)
#define G_MATRIX_H
void   free_two_elec_int(void);
void   precalculate_two_elec_int(size_t threads);
void   build_G_matrix_in_memory(double *G, const double *P);
void   build_G_matrix_on_the_fly(double *G, const double *P);
#endif
