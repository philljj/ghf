#if !defined (FOCK_H)
#define FOCK_H
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
    bool       diis;	
    bool       damping;	
    void     (*make_fock)(struct fock_matrix_t *f, const double *P);	
};
#endif
