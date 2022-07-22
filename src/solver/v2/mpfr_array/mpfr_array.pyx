# cython: profile=True

# S. El-Showk, S. Rychkov May-June 2013
# contains mpfr array class (called ndarray for analogy with numpy array)


from libc.stdlib cimport malloc, free 

from c_mpfr cimport *    
from c_ludcmp_mpfr cimport * 

import numbers 
import numpy as np 
import solver.v2.prec_float.prec_float as PF
cimport solver.v2.prec_float.prec_float as PF

import utils.stats as stats


cdef extern from "mpfr_array.h": 
    cdef int* one_arr
     
def getshape(plist):
    """recursively return the shape of a python list"""
    if not isinstance(plist, list):
        return (),[plist]
    shape=(len(plist),)
    subshape,flatlist=getshape(plist[0])
    oldsubshape=subshape
    flatlist=[]
    # only grow the shape if all sublists have the same shape 
    for sl in plist:
        subshape,fl1=getshape(sl)
        if subshape != oldsubshape:
            return (),[]
        flatlist += fl1
    return shape + subshape, flatlist

cpdef ndarray array(plist, prec=-1):
    """plist is a python list of doubles, strings, or prec_float objects"""
    cdef int i
    if isinstance(plist, list):
        shape,flatlist=getshape(plist)
        if shape==():
            raise ValueError("array initialized with badly shaped list")
        if prec > 0:
            ret=ndarray(shape, prec=prec)
        else:
            ret=ndarray(shape)
        for i in range(ret.size):
            ret.flatset(i, flatlist[i])
        return ret
    elif isinstance(plist, np.ndarray):
        ret=ndarray(shape=plist.shape)
        flist = plist.flatten()
        for i in range(ret.size):
            mpfr_set_d(ret.data[i], flist[i], MPFR_RNDD)
        return ret
    else:
        raise ValueError("array initialized with non-list")

cpdef ndarray empty(shape, prec=-1):
    if prec > 0:
            ret=ndarray(shape, prec=prec)
    else:
            ret=ndarray(shape)
    return ret

cpdef ndarray zeros(shape, prec=-1):
    if prec > 0:
            ret=ndarray(shape, prec=prec)
    else:
            ret=ndarray(shape)
    for i in range(ret.size):
            mpfr_set_zero(ret.data[i],1)
    return ret

cpdef ndarray ones(shape, prec=-1):
    if prec > 0:
            ret=ndarray(shape, prec=prec)
    else:
            ret=ndarray(shape)
    for i in range(ret.size):
            mpfr_set_si(ret.data[i],<long int>1,MPFR_RNDD)
    return ret 

cpdef ndarray empty_like(ndarray a):
    return ndarray(prec=a.prec, shape=a.getshape())

cpdef ndarray dot(ndarray array1, ndarray array2):
    return array1.dot(array2)
 
cpdef to_double (ndarray array):
    """returns array of rounded doubles, or a single double if array has length 1"""
    if array.size == 1: # single value
        return mpfr_get_d(array.data[0],MPFR_RNDD)
    else:
        temp = np.empty(array.getshape())
        for i in range(array.size):
            temp.flat[i] = mpfr_get_d(array.data[i],MPFR_RNDD)
        return temp
    
cpdef PF.prec_float to_pf(ndarray array):
        cdef PF.prec_float ret
        if (array.size != 1 ):
            print"size=", array.size
            raise TypeError("only for 1-element arrays can transform into prec_float")
        ret = PF.prec_float(0, array.prec)
        mpfr_set(ret.data, array.data[0],MPFR_RNDD)
        return ret
    
cpdef list to_pf_array(ndarray array):
        cdef list ret
        cdef int i
        ret = [PF.prec_float(0, array.prec) for i in range(array.size)]
        for i in range(array.size):
            mpfr_set((<PF.prec_float>ret[i]).data, array.data[i],MPFR_RNDD)
        return ret

# some memory efficient operators

# add B to C and assign to A -- all must be the same size
cpdef add_assign(ndarray A, ndarray B, ndarray C):
    for i in range(A.size):
        mpfr_add(A.data[i], B.data[i],  C.data[i], MPFR_RNDD)
        stats.ndarray_adds+=1

# element-wise multiply B with C and assign to A -- all must be the same size
cpdef mult_assign(ndarray A, ndarray B, ndarray C):
    for i in range(A.size):
        mpfr_mul(A.data[i], B.data[i],  C.data[i], MPFR_RNDD)
        stats.ndarray_adds+=1


# we need to deal with the fact that inverse prec is now hard coded!!
#cpdef ndarray inv(ndarray mat):
#    cdef ndarray ret
#    cdef int status
#    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
#        raise ValueError("Argument of inverse must be a square matrix.")
#    ret=ndarray(mat.getshape(), prec=mat.prec)
#    status = inverse(mat.data, mat.shape[0], ret.data)
#    if status == 0:
#        raise ValueError("Matrix Inversion failed")
#    return ret

#cpdef int inv_assign(ndarray mat, ndarray ret):
#    cdef int status
#    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
#        raise ValueError("Argument of inverse must be a square matrix.")
#    if ret.ndim != 2 or ret.shape[0] != ret.shape[1] or ret.shape[0] != mat.shape[0]:
#        raise ValueError("Return argument of inverse must be a square matrix of size\
#        equal to the initial matrix.")
#    status = inverse(mat.data, mat.shape[0], ret.data)
#    if status == 0:
#        raise ValueError("Matrix Inversion failed")
#    return 1

cpdef int inv_assign(ndarray mat, ndarray ret, ndarray scratch):
    cdef int status
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError("Argument of inverse must be a square matrix.")
    if ret.ndim != 2 or ret.shape[0] != ret.shape[1] or ret.shape[0] != mat.shape[0]:
        raise ValueError("Return argument of inverse must be a square matrix of size\
        equal to the initial matrix.")
    status = inverse(mat.data, mat.shape[0], ret.data, scratch.data)
    if status == 0:
        raise ValueError("Matrix Inversion failed")
    return 1

# LU decomposition block 
# we have to create a class to hold LU decomposition data to allow for common interface with double
# precision version

cdef class LUdcmp_data:
    cdef ndarray lu
    cdef int* indx
    
    def __cinit__(self, ndarray LU):
        self.lu = LU
        self.indx =  <int*> malloc(LU.shape[0] * sizeof(int))
        
    def __dealloc__(self):
        free(self.indx)
        #self.lu.__dealloc__()
        
        
cpdef LUdcmp_data lu_factor(ndarray mat):
    cdef ndarray lu
    cdef LUdcmp_data lu_data
    cdef int* indx
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError("Argument of lu_factor must be a square matrix.")
    lu = ndarray(mat.getshape(), prec = mat.prec)
    indx =  <int*> malloc(mat.shape[0] * sizeof(int))
    status = LUdcmp(mat.data, mat.shape[0], lu.data, indx)
    if status == 0:
        raise ValueError("LU decomposition failed")
    lu_data = LUdcmp_data(lu)
    lu_data.indx = indx
    return lu_data

cpdef int lu_factor_assign(ndarray mat, LUdcmp_data lu_data):
    cdef ndarray lu
    status = LUdcmp(mat.data, mat.shape[0], (<ndarray>lu_data.lu).data, lu_data.indx)
    if status == 0:
        raise ValueError("LU decomposition failed")
    return 1

cpdef ndarray lu_solve(LUdcmp_data lu_data, ndarray b, int trans):
    cdef int i
    cdef ndarray x
    if b.ndim != 1 or b.shape[0] != lu_data.lu.shape[0]:
        raise ValueError("Wrong size array(s) in lu_solve")
    x = ndarray(b.getshape(), prec = b.prec);
    for i in range(0,b.shape[0]):
        x[i] = b[i] #x.flatset1(i, b.data[i])
    if trans == 0:
        LUsolve_in_place(lu_data.lu.data, lu_data.lu.shape[0], lu_data.indx, x.data)
    elif trans == 1:
        ULsolve_in_place(lu_data.lu.data, lu_data.lu.shape[0], lu_data.indx, x.data)
    else:
        raise ValueError("trans must be 0 or 1")
    return x
        
# end of LU decomposition block 

cpdef ndarray copy(ndarray b): 
    cdef int i
    cdef ndarray x
    x = ndarray(b.getshape(), prec = b.prec);
    for i in range(0,b.size):
        x.flatset1(i, b.data[i])
    return x
    
class ndarray_iter:
    
    def __init__(self, ndarray a):
        self.array = a
        self.size=a.size
        self._cntr=0

    # def next(self):
    def __next__(self):   # Carlos edit
        if self._cntr < self.size:
            self._cntr += 1
            return self.array[self._cntr-1]
        else:
            raise StopIteration

cpdef ndarray concatenate(arraylist):
    cdef int i
    cdef ndarray x
    a = <ndarray>arraylist[0]
    b = <ndarray>arraylist[1]
    if a.ndim !=1 or b.ndim !=1:
        raise ValueError("can concatenate only one-dim arrays")
    x = ndarray((a.size+b.size,), prec = b.prec);
    for i in range(0,a.size):
        x.flatset1(i, a.data[i])
    for i in range(0,b.size):
        x.flatset1(a.size+i, b.data[i])
    return x

cdef char* printbuf
cdef int printbuflength = 10000
printbuf = <char*> malloc( 212*4 * printbuflength * sizeof(char)) 



####################  BEGIN --- CLASS ------NDARRAY  #############################

cdef class ndarray:

#    # vars same as numpy.ndarray
#    cdef mpfr_t* data
#    cdef int* shape
#    cdef int* strides
#    cdef public int ndim
#    cdef public int size
#    
#    cdef int prec
#    # track our parent class
#    cdef ndarray parent
#    # the number of parent arrays above this
#    cdef int level
 
    def __cinit__(self, shape=(1,), ndarray parent=None, int p_ind=-1, int prec=212):
        """initialize an mpfr array of given dims with precision prec"""
        self.data=NULL
        self.shape=NULL
        self.strides=NULL
        self.parent=parent
        #self.bufferflag = 0 # if buffer for tolist has been allocated
        stats.ndarray_init += 1 
        if parent is None:
            self.level=0
            self.parent=self
            self.size=1
            ndim=len(shape)
            self.ndim=ndim
            self.shape=<int*> malloc(ndim * sizeof(int))
            self.strides=<int*> malloc(ndim * sizeof(int))
            for i in range(ndim):
                self.shape[i]= shape[i]
            self.prec=prec
            for i in range(ndim):
                self.size = self.size * self.shape[i]
                self.strides[i] = 1
            #self.blocksize = [1. for x in range(ndim)]
            # for array indexing compute the block size for each index
            # e.g. for [10,5,30] we should get [150, 30, 1]
            for i in range(1,ndim):
                self.strides[ndim-i-1]=self.shape[ndim-i]*self.strides[ndim-i]
            data=<mpfr_t*> malloc(self.size*sizeof(mpfr_t))
            self.data=data
            for i in range(self.size):
                mpfr_init2(self.data[i], self.prec)
            stats.ndarray_mem += (2 * self.ndim * sizeof(int) +
                        self.size*sizeof(mpfr_t))
        # this is a sub-array of a larger array so just use its data
        else:
            self.level=parent.level+1
            # if our parent is a vector then we're a scalar
            if parent.ndim == 1:
                self.shape = one_arr
                self.ndim = 1
                self.strides = one_arr
            else:
                self.shape = &parent.shape[1]
                self.ndim = parent.ndim - 1
                self.strides = &parent.strides[1]
            #self.print_array(self.dims, self.ld)
            #self.print_array(parent.dims, parent.ld)
            self.prec=parent.prec
            self.size=1
            for i in range(self.ndim):
                self.size = self.size * self.shape[i]
            self.data = &parent.data[p_ind * parent.strides[0]]
        
    # sub-arrays reference their parents so the parent should never be garbage
    # collected before the sub-array.  Thus we only do garbage collection if
    # we're a top-level parent (self.parent == None).
    #
    # NOTE: ignoring boundary case of slices which need to de-alloc their shape
    # and strides
    #
    def __dealloc__(ndarray self):
        #self.top_parent().refs -= 1
        stats.ndarray_del += 1
        if self.level == 0:
            #print 'deallocating'
            if self.data != NULL:
                for i in range(self.size): 
                    mpfr_clear(self.data[i])
                free(self.data)
            if self.shape != NULL and self.shape != one_arr:
                free(self.shape)
            if self.strides != NULL and self.shape != one_arr:
                free(self.strides)
            stats.ndarray_free += (2 * self.ndim * sizeof(int) +
                        self.size*sizeof(mpfr_t))

    ##  Overloaded functions


    def __setitem__(ndarray self, int i, value):
        """Set array element from string or number or set subarray from another
        (mpfr) ndarray. """
        cdef int j,k, sd
        # set from another ndarray
        if isinstance(value, ndarray):
            sd=self.ndim
            ndval = <ndarray> value
            if self.ndim == 1 and ndval.size == 1: # exceptional case of RHS being a single element
                mpfr_set(self.data[i], ndval.data[0], MPFR_RNDD)
            elif self.ndim-1 != ndval.ndim or self.strides[0] != ndval.size:
                raise TypeError(
                    "can't set from incorrectly shaped ndarray (err 1)")
            else:
                for k in range(ndval.ndim):
                    if self.shape[k+1] != ndval.shape[k]:
                        raise TypeError(
                            "can't set from incorrectly shaped ndarray (err 2)")
                for j in range(self.strides[0]):
                    mpfr_set(self.data[self.strides[0]*i + j],
                        ndval.data[j], MPFR_RNDD)
        # if we're a scalar set from string or number
        elif self.ndim == 1:
            if isinstance(value, numbers.Real):
                mpfr_set_d(self.data[i], value, MPFR_RNDD)
            elif isinstance(value, str):
               # mpfr_set_str(self.data[i], value, 10, MPFR_RNDD)
                bytes_str= bytes(value,'utf-8')  # Carlos edit
                mpfr_set_str(self.data[i], bytes_str, 10, MPFR_RNDD)   # Carlos
            else:
                raise TypeError(
                    "cannot initialize len > 1 array from string or number")
        else:
            raise TypeError(
                "can only initialize mpfr from number, string, or (mpfr) array")


    def __getitem__(self, int i):
        # should throw exception
        if i > self.shape[0]:
            raise IndexError
        return ndarray(parent=self, p_ind=i)

    def __neg__(ndarray self):
        """Sign flip"""
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            mpfr_neg(ret.data[i], self.data[i], MPFR_RNDD)
        return ret    
    
    def __add__(ndarray self, ndarray other):
        """Element-wise addition."""
        if other.size != self.size:
            raise ValueError("adding different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            mpfr_add(ret.data[i], self.data[i],  other.data[i], MPFR_RNDD)
            stats.ndarray_adds+=1
        return ret    

    def __sub__(ndarray self, ndarray other):
        """Element-wise subtraction."""
        if other.size != self.size:
            raise ValueError("subtracting different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            mpfr_sub(ret.data[i], self.data[i],  other.data[i], MPFR_RNDD)
        return ret    

    def __mul__(ndarray self, ndarray other):
        """Element-wise multiplication."""
        if other.size != self.size:
            raise ValueError("multiplying different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            mpfr_mul(ret.data[i], self.data[i],  other.data[i], MPFR_RNDD)
            stats.ndarray_mults+=1
        return ret    

     #def __div__(ndarray self, ndarray other):
    def __truediv__(ndarray self, ndarray other): # Carlos edit
        """Element-wise division."""
        if other.size != self.size:
            raise ValueError("dividing different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            mpfr_div(ret.data[i], self.data[i],  other.data[i], MPFR_RNDD)
        return ret    

    def __richcmp__(ndarray self, ndarray value, int cmptype):
        """Element-wise comparison
        op cmptype
        <	0   
        ==	2   
        >	4   
        <=	1   
        !=	3   
        >=	5"""
        if value is None:
            return False
        cdef int cmp
        if self.size != 1 or value.size !=1:
            raise ValueError("Comparison defined only for ndarrays of size 1 (item 1)")
       
        cmp = mpfr_cmp(self.data[0], value.data[0])
        
        if cmptype == 0:
            if cmp < 0: return True
            else: return False
        if cmptype == 2:
            if cmp == 0: return True
            else: return False
        if cmptype == 4:
            if cmp > 0: return True
            else: return False
        if cmptype == 1:
            if cmp <= 0: return True
            else: return False
        if cmptype == 3:
            if cmp != 0: return True
            else: return False
        if cmptype == 5:
            if cmp >= 0: return True
            else: return False
    
    def __iter__(ndarray self):
        return ndarray_iter(self)
    
#     def __reduce__(self): # Carlos Edit
#         return self.mat_str()


    ## math functions

    cdef ndarray dot(ndarray self, ndarray b):
        if self.ndim ==1 and b.ndim==1:
            return self._vv_dot(b)
        if self.ndim ==2 and b.ndim==1:
            return self._mv_dot(b)
        if self.ndim ==1 and b.ndim==2:
            return self._vm_dot(b)
        if self.ndim ==2 and b.ndim==2:
            return self._mm_dot(b)

    cdef ndarray _vv_dot(ndarray self, ndarray b):
        """compute c=a.b"""
        cdef int i, bsize
        cdef mpfr_t t
        cdef ndarray r 
        bsize=b.size
        
        if self.size != bsize:
            raise ValueError("cannot dot different shape arrays")
        mpfr_init2(t, self.prec)
        mpfr_set_d(t, 0, MPFR_RNDD)
        r=ndarray(prec=self.prec)
        r[0]=0
        for i in range(bsize):
            mpfr_mul(t, self.data[i],  b.data[i], MPFR_RNDD)
            mpfr_add(r.data[0], r.data[0], t, MPFR_RNDD)
            stats.ndarray_mults+=1
            stats.ndarray_adds+=1
        mpfr_clear(t)
        return r
    
    cdef _mv_dot(ndarray self, ndarray b):
        """compute c=a.b with mat a and vec b."""
        cdef mpfr_t t
        cdef int i, j, nr, nc, inc
        cdef ndarray r
        
        nr=self.shape[0]
        nc=self.shape[1]
        #print "self: ",(self.ndim, self.shape[0], self.size)
        #print "b: ",(nr, nc, b.size)
        if nc != b.shape[0]:
            raise ValueError("cannot dot different shaped arrays")
        # a temporary mpfr
        mpfr_init2(t, self.prec)
        mpfr_set_d(t, 0, MPFR_RNDD)
        # result
        r=ndarray(prec=self.prec, shape=(nr,))
        for i in range(nr):
            r[i]=0
            inc=i*nc
            for j in range(nc):
                #print '(',i,',',j,')'
                #print self.array[j]
                #print b.array[i*nr + j], 
                mpfr_mul(t, self.data[inc + j],  b.data[j], MPFR_RNDD)
                mpfr_add(r.data[i], r.data[i], t, MPFR_RNDD)
                stats.ndarray_mults+=1
                stats.ndarray_adds+=1
        mpfr_clear(t)
        return r
    
    cdef _vm_dot(ndarray self, ndarray b):
        """compute c=a.b with vec a and matrix b."""
        cdef mpfr_t t
        cdef ndarray r
        cdef int i, j, nr, nc, inc
        nr=b.shape[0]
        nc=b.shape[1]
        #print "self: ",(self.ndim, self.shape[0], self.size)
        #print "b: ",(nr, nc, b.size)
        if nr != self.shape[0]:
            raise ValueError("cannot dot different shaped arrays")
        # a temporary mpfr
        mpfr_init2(t, self.prec)
        mpfr_set_d(t, 0, MPFR_RNDD)
        # result
        r=ndarray(prec=self.prec, shape=(nc,))
        for i in range(nc):
            r[i]=0 
            for j in range(nr):
                mpfr_mul(t, self.data[j],  b.data[j*nc + i], MPFR_RNDD)
                mpfr_add(r.data[i], r.data[i], t, MPFR_RNDD)
                stats.ndarray_mults+=1
                stats.ndarray_adds+=1
        mpfr_clear(t)
        return r
 
    def _mm_dot(self, ndarray b):
        """compute c=a.b with a, b, c matrices."""
        cdef mpfr_t t
        cdef int i, j, nr, nc, inc
        nr=self.shape[0]
        nc=self.shape[1]
        br=b.shape[0]
        bc=b.shape[1]
        #print "self: ",(self.ndim, self.shape[0], self.size)
        #print "b: ",(nr, nc, b.size)
        if nc != br:
            raise ValueError("array indices don't line up")
        # a temporary mpfr
        mpfr_init2(t, self.prec)
        mpfr_set_d(t, 0, MPFR_RNDD)
        # result
        r=ndarray(prec=self.prec, shape=(nr,bc))
        for i in range(nr):
            for j in range(bc):
                r[i][j]=0
                inc=i*nc
                for k in range(nc):
                    #print '(',i,',',j,')'
                    #print self.array[j]
                    #print b.array[i*nr + j], 
                    mpfr_mul(t, self.data[inc + k],  b.data[k*bc + j], MPFR_RNDD)
                    mpfr_add(r.data[i*bc + j], r.data[i*bc + j], t, MPFR_RNDD)
                    stats.ndarray_mults+=1
                    stats.ndarray_adds+=1
        mpfr_clear(t)
        return r
    
    #------------------- dot_assign functions
    
    cdef _vm_dot_assign(ndarray self, ndarray b, ndarray r):
        """compute c=a.b with vec a and matrix b and assign result to r"""
        cdef mpfr_t t
        cdef int i, j, nr, nc, inc
        nr=b.shape[0]
        nc=b.shape[1]
        #print "self: ",(self.ndim, self.shape[0], self.size)
        #print "b: ",(nr, nc, b.size)
        if nr != self.shape[0]:
            raise ValueError("cannot dot different shaped arrays")
        if nc != r.shape[0]:
            print self.shape[0], nr,nc,r.shape[0]
            raise ValueError("cannot store in a different shaped arrays")
        # a temporary mpfr
        mpfr_init2(t, self.prec)
        mpfr_set_d(t, 0, MPFR_RNDD)
        # result
        for i in range(nc):
            r[i]=0
            for j in range(nr):
                mpfr_mul(t, self.data[j],  b.data[j*nc + i], MPFR_RNDD)
                mpfr_add(r.data[i], r.data[i], t, MPFR_RNDD)
                stats.ndarray_mults+=1
                stats.ndarray_adds+=1
        mpfr_clear(t)
        return r


    ## Printing functions -- make recursive??

    def __str__(self):
        if self.ndim == 1:
            return self.vec_str()
        if self.ndim == 2:
            return self.mat_str()
    
    def __repr__(self): 
        return self.__str__()
    
    
    def tolist1(self):
        global printbuf
        global printbuflength 
    
        cdef int i
        if self.ndim == 1:
            if printbuflength < self.size:
                raise ValueError("Insufficient print buffer")
            else:
                for i in range(self.size):
                    mpfr_sprintf(printbuf+i*212*4,"%.Re",self.data[i])
                return [(printbuf+i*212*4).decode('ascii') for i in range(self.size)]
        else:
            raise TypeError("tolist function defined only for 1dim ndarrays")
            
            
    def tolist(self):
        if self.ndim == 1:
            return [self[i].vec_str() for i in range(self.size)]
        else:
            raise TypeError("tolist function defined only for 1dim ndarrays")
    #    
    #def tolist(self, buf = ""):
    #    if self.ndim == 1:
    #        if buf == "":
    #            buf = [<char*> malloc(2*self.prec * sizeof(char)) for i in range(self.size)]
    #        for i in range(self.size):
    #            mpfr_sprintf(<char*>(buf[i]),"%.Re",self.data[i])
    #        return [buf[i] for i in range(self.size)]
    #    else:
    #       raise TypeError("tolist function defined only for 1dim ndarrays")
    
    def vec_str(self):
        if self.size ==0:
            return ''
        cdef int bs, n, cn, i
        n=0
        cn=1
        # '[]' (entries +', ') + a little padding
        bs= 2+ (self.prec+2) * (self.size + 4)
        cdef char* buf = <char*> malloc(bs * sizeof(char))
        if self.size == 1:
            mpfr_sprintf(buf,"%.Re",self.data[0])
        else:
            buf[0]='['
            for i in range(self.size-1):
                n=mpfr_sprintf(buf + cn, "%.Re, ", self.data[i]) 
                cn +=n
            n=mpfr_sprintf(buf + cn, "%.Re]", self.data[self.size-1])
        #return buf
        return buf.decode('ascii')  # Carlos edit


    def mat_str(self):
        if self.size ==0:
            return ''
        cdef int bs, n, cn, i
        n=0
        cn=1
        # '[]' + '[]' for each row + (entries +', ') + a little padding
        bs= 2+ 2* self.shape[0] + (self.prec+3) * (self.size + 4)
        cdef char* buf = <char*> malloc(bs * sizeof(char))
        buf[0]='['
        for i in range(self.shape[0]):
            buf[cn]='['
            cn+=1
            for j in range(self.shape[1]):
                n=mpfr_sprintf(buf + cn, "%.Re", 
                        self.data[i * self.strides[0] + j]) 
                cn +=n
                # end of matrix
                if i == self.shape[0]-1 and j == self.shape[1]-1:
                    buf[cn]=']'
                    cn+=1
                # end of row
                elif j == self.shape[1]-1:
                    buf[cn]=']'
                    buf[cn+1]=','
                    buf[cn+2]='\n'
                    cn+=3
                else:
                    buf[cn]=','
                    buf[cn+1]=' '
                    cn+=2

        # do the null termination correctly
        sprintf(buf +cn, ']')
        #return buf 
        return buf.decode('ascii')   # Carlos edit
    
    # utility functions
   
    def flatset(self, int i, value):
        """Set element of data array indexed _directly_ by i (i.e. forget shape)"""
        if isinstance(value, numbers.Real):
            mpfr_set_d(self.data[i], value, MPFR_RNDD)
        elif isinstance(value, str):
            #mpfr_set_str(self.data[i], value, 10, MPFR_RNDD)
            bytes_str= bytes(value,'utf-8')  # Carlos edit
            mpfr_set_str(self.data[i], bytes_str, 10, MPFR_RNDD)
        elif isinstance(value, PF.prec_float):
            mpfr_set(self.data[i], (<PF.prec_float>value).data, MPFR_RNDD) 
        else:
            raise TypeError("can only initialize mpfr from number or string or prec_float")
            
    cdef flatset1(self, int i, mpfr_t value):
        mpfr_set(self.data[i], value, MPFR_RNDD)
     
    def getshape(self):
        return self.array2set(self.shape, self.ndim)

    # General tensor functions
    
    
    def set(self, inds, value):
        """Setter for general shape array.  Set the inds entry to value. 
        inds should be [i,j,k, ...]"""
        i=0
    

    def bs(self):
        for i in range(self.ndim):
            print self.strides[i]

    cdef print_array(self, int* arr, int ln):
        for i in range(ln):
            print arr[i]

    cdef array2set(self, int* arr, int ln):
        r=(arr[0],)
        for i in range(1,ln):
            r+=(arr[i],)
        return r
        #return (arr[i] for i in range(ln))


####################  END --- CLASS ------NDARRAY  #############################



