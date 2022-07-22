from c_mpfr cimport *
#cimport bootstrap.prec_float.prec_float as PF

cdef class ndarray:
    cdef mpfr_t* data
    cdef int* shape
    cdef int* strides
    cdef public int ndim
    cdef public int size   
    cdef int prec
    # track our parent class
    cdef ndarray parent
    # the number of parent arrays above this
    cdef int level
    cdef list buf
    cdef int bufferflag
    
    cdef ndarray dot(ndarray self, ndarray b)
    cdef ndarray _vv_dot(ndarray self, ndarray b)
    cdef _mv_dot(ndarray self, ndarray b)
    cdef _vm_dot(ndarray self, ndarray b)
    cdef _vm_dot_assign(ndarray self, ndarray b, ndarray r)
    cdef flatset1(self, int i, mpfr_t value)
    cdef print_array(self, int* arr, int ln)
    cdef array2set(self, int* arr, int ln)
    #cpdef PF.prec_float to_pf(self)
