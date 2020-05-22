#if !defined(SPIN_TYPES_H)
#define SPIN_TYPES_H

typedef spinor_t {
    // Spin up and down coefficient matrix.
    double * u;
    double * d;
    size_t   n_basis;
}

typedef spin_matrix_t {
    // 4 component generalized matrix.
    double * uu;
    double * ud;
    double * du;
    double * dd;
    size_t   n_basis;
}

#endif /* if !defined (SPIN_TYPES_H)*/
