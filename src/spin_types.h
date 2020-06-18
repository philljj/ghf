#if !defined(SPIN_TYPES_H)
#define SPIN_TYPES_H

struct spinor_t {
    // Spin up and down coefficient matrix.
    double * u;
    double * d;
    size_t   n_basis;
}

typedef struct spinor_t spinor_t;

struct spin_matrix_t {
    // 4 component generalized matrix.
    double *  uu;
    double *  ud;
    double *  du;
    double *  dd;
    size_t    n_basis;
    size_t    n_ele_u; /* only used in UHF */
    size_t    n_ele_d; /* only used in UHF */
    size_t    n_ele;   /* only used in RHF or GHF */
}

typedef struct spin_matrix_t spin_matrix_t;

#endif /* if !defined (SPIN_TYPES_H)*/
