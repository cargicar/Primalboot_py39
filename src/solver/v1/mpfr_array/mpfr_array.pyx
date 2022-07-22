from libc.stdlib cimport malloc, free
from libc.stdio cimport stdout, sprintf

cimport mpfr as m

import numbers
import numpy as np

import traceback

cdef extern from "mpfr_array.h":
    cdef int* one_arr

def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]
    
def getshape(plist):
    """recusrively return the shape of a python list"""
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
    """plist is a python list of doubles, strings, or mpfr_t objects"""
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
            m.mpfr_set_d(ret.data[i], flist[i], m.MPFR_RNDD)
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
            m.mpfr_set_zero(ret.data[i],1)
    return ret

cpdef ndarray ones(shape, prec=-1):
    if prec > 0:
            ret=ndarray(shape, prec=prec)
    else:
            ret=ndarray(shape)
    for i in range(ret.size):
            m.mpfr_set_si(ret.data[i],<long int>1,m.MPFR_RNDD)
    return ret 

#cdef array2set(int* arr, int ln):
#        r=(arr[0],)
#        for i in range(1,ln):
#            r+=(arr[i],)
#        return r
  

cpdef ndarray empty_like(ndarray a):
    return ndarray(prec=a.prec, shape=a.getshape())

cpdef ndarray dot(ndarray array1, ndarray array2):
    return array1.dot(array2)
 
cpdef to_double (ndarray array):
    """returns array of rounded doubles, or a single double if array has length 1"""
    if array.size == 1: # single value
        return m.mpfr_get_d(array.data[0],m.MPFR_RNDD)
    else:
        temp = np.empty(array.getshape())
        for i in range(array.size):
            temp.flat[i] = m.mpfr_get_d(array.data[i],m.MPFR_RNDD)
        return temp
    
# we need to deal with the fact that inverse prec is now hard coded!!
cpdef ndarray inv(ndarray mat):
    cdef ndarray ret
    cdef int status
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError("Argument of inverse must be a square matrix.")
    ret=ndarray(mat.getshape(), prec=mat.prec)
    status = m.inverse(mat.data, mat.shape[0], ret.data)
    if status == 0:
        raise ValueError("Matrix Inversion failed")
    return ret

cpdef int inv_assign(ndarray mat, ndarray ret):
    cdef int status
    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError("Argument of inverse must be a square matrix.")
    if ret.ndim != 2 or ret.shape[0] != ret.shape[1] or ret.shape[0] != mat.shape[0]:
        raise ValueError("Return argument of inverse must be a square matrix of size\
        equal to the initial matrix.")
    status = m.inverse(mat.data, mat.shape[0], ret.data)
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
    status = m.LUdcmp(mat.data, mat.shape[0], lu.data, indx)
    if status == 0:
        raise ValueError("LU decomposition failed")
    lu_data = LUdcmp_data(lu)
    lu_data.indx = indx
    return lu_data
    
cpdef ndarray lu_solve(LUdcmp_data lu_data, ndarray b, int trans):
    cdef int i
    cdef ndarray x
    if b.ndim != 1 or b.shape[0] != lu_data.lu.shape[0]:
        raise ValueError("Wrong size array(s) in lu_solve")
    x = ndarray(b.getshape(), prec = b.prec);
    for i in range(0,b.shape[0]):
        x[i] = b[i] #x.flatset1(i, b.data[i])
    if trans == 0:
        m.LUsolve_in_place(lu_data.lu.data, lu_data.lu.shape[0], lu_data.indx, x.data)
    elif trans == 1:
        m.ULsolve_in_place(lu_data.lu.data, lu_data.lu.shape[0], lu_data.indx, x.data)
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

    def next(self):
        if self._cntr < self.size:
            self._cntr += 1
            return self.array[self._cntr-1]
        else:
            raise StopIteration



        

####################  BEGIN --- CLASS ------NDARRAY  #############################

cdef class ndarray:

    # vars same as numpy.ndarray
    cdef m.mpfr_t* data
    cdef int* shape
    cdef int* strides
    cdef public int ndim
    cdef public int size
    
    cdef int prec
    # track our parent class
    cdef ndarray parent
    # the number of parent arrays above this
    cdef int level
 
    def __cinit__(self, shape=(1,), ndarray parent=None, int p_ind=-1, int prec=212):
        """initialize an mpfr array of given dims with precision prec"""
        self.data=NULL
        self.shape=NULL
        self.strides=NULL
        self.parent=parent
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
            data=<m.mpfr_t*> malloc(self.size*sizeof(m.mpfr_t))
            self.data=data
            for i in range(self.size):
                m.mpfr_init2(self.data[i], self.prec)
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
        if self.level == 0:
            #print 'deallocating'
            if self.data != NULL:
                for i in range(self.size): 
                    m.mpfr_clear(self.data[i])
                free(self.data)
            if self.shape != NULL and self.shape != one_arr:
                free(self.shape)
            if self.strides != NULL and self.shape != one_arr:
                free(self.strides)

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
                m.mpfr_set(self.data[i], ndval.data[0], m.MPFR_RNDD)
            elif self.ndim-1 != ndval.ndim or self.strides[0] != ndval.size:
                raise TypeError(
                    "can't set from incorrectly shaped ndarray (err 1)")
            else:
                for k in range(ndval.ndim):
                    if self.shape[k+1] != ndval.shape[k]:
                        raise TypeError(
                            "can't set from incorrectly shaped ndarray (err 2)")
                for j in range(self.strides[0]):
                    m.mpfr_set(self.data[self.strides[0]*i + j],
                        ndval.data[j], m.MPFR_RNDD)
        # if we're a scalar set from string or number
        elif self.ndim == 1:
            if isinstance(value, numbers.Real):
                m.mpfr_set_d(self.data[i], value, m.MPFR_RNDD)
            elif isinstance(value, str):
                m.mpfr_set_str(self.data[i], value, 10, m.MPFR_RNDD)
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
            m.mpfr_neg(ret.data[i], self.data[i], m.MPFR_RNDD)
        return ret    
    
    def __add__(ndarray self, ndarray other):
        """Element-wise addition."""
        if other.size != self.size:
            raise ValueError("adding different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            m.mpfr_add(ret.data[i], self.data[i],  other.data[i], m.MPFR_RNDD)
        return ret    

    def __sub__(ndarray self, ndarray other):
        """Element-wise subtraction."""
        if other.size != self.size:
            raise ValueError("subtracting different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            m.mpfr_sub(ret.data[i], self.data[i],  other.data[i], m.MPFR_RNDD)
        return ret    

    def __mul__(ndarray self, ndarray other):
        """Element-wise multiplication."""
        if other.size != self.size:
            raise ValueError("multiplying different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            m.mpfr_mul(ret.data[i], self.data[i],  other.data[i], m.MPFR_RNDD)
        return ret    

    def __div__(ndarray self, ndarray other):
        """Element-wise division."""
        if other.size != self.size:
            raise ValueError("dividing different sized ndarray")
        # would it be faster to deep-copy?
        ret=ndarray(prec=self.prec, shape=self.array2set(self.shape, self.ndim))
        for i in range(self.size):
            m.mpfr_div(ret.data[i], self.data[i],  other.data[i], m.MPFR_RNDD)
        return ret    

#    def __richcmp__(ndarray self, value, int cmptype):
#        """Element-wise comparison
#        op cmptype
#        <	0   
#        ==	2   
#        >	4   
#        <=	1   
#        !=	3   
#        >=	5"""
#        if self.size != 1:
#            raise ValueError("Comparison defined only for ndarrays of size 1 (item 1)")
#       
#        if isinstance(value,int):
#            cmp = m.mpfr_cmp_si(self.data[0], <long int> value)
#        elif isinstance(value,ndarray):
#            ndval = <ndarray> value
#            if ndval.size != 1:
#                raise ValueError("Comparison defined only for ndarrays of size 1 (item 2)")
#            cmp = m.mpfr_cmp(self.data[0], ndval.data[0])
#        else:
#            raise ValueError("Comparison defined only for item 2 integer or ndarray")
#        if cmptype == 0:
#            if cmp < 0: return True
#            else: return False
#        if cmptype == 2:
#            if cmp == 0: return True
#            else: return False
#        if cmptype == 4:
#            if cmp > 0: return True
#            else: return False
#        if cmptype == 1:
#            if cmp <= 0: return True
#            else: return False
#        if cmptype == 3:
#            if cmp != 0: return True
#            else: return False
#        if cmptype == 5:
#            if cmp >= 0: return True
#            else: return False

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
       
        cmp = m.mpfr_cmp(self.data[0], value.data[0])
        
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
        cdef m.mpfr_t t
        cdef ndarray r
        bsize=b.size
        
        if self.size != bsize:
            raise ValueError("cannot dot different shape arrays")
        m.mpfr_init2(t, self.prec)
        m.mpfr_set_d(t, 0, m.MPFR_RNDD)
        r=ndarray(prec=self.prec)
        r[0]=0
        for i in range(bsize):
            m.mpfr_mul(t, self.data[i],  b.data[i], m.MPFR_RNDD)
            m.mpfr_add(r.data[0], r.data[0], t, m.MPFR_RNDD)
        return r
    
    cdef _mv_dot(ndarray self, ndarray b):
        """compute c=a.b with mat a and vec b."""
        cdef m.mpfr_t t
        cdef int i, j, nr, nc, inc
        cdef ndarray r
        
        nr=self.shape[0]
        nc=self.shape[1]
        #print "self: ",(self.ndim, self.shape[0], self.size)
        #print "b: ",(nr, nc, b.size)
        if nc != b.shape[0]:
            raise ValueError("cannot dot different shaped arrays")
        # a temporary mpfr
        m.mpfr_init2(t, self.prec)
        m.mpfr_set_d(t, 0, m.MPFR_RNDD)
        # result
        r=ndarray(prec=self.prec, shape=(nr,))
        for i in range(nr):
            r[i]=0
            inc=i*nc
            for j in range(nc):
                #print '(',i,',',j,')'
                #print self.array[j]
                #print b.array[i*nr + j], 
                m.mpfr_mul(t, self.data[inc + j],  b.data[j], m.MPFR_RNDD)
                m.mpfr_add(r.data[i], r.data[i], t, m.MPFR_RNDD)
        return r
    
    cdef _vm_dot(ndarray self, ndarray b):
        """compute c=a.b with vec a and matrix b."""
        cdef m.mpfr_t t
        cdef ndarray r
        cdef int i, j, nr, nc, inc
        nr=b.shape[0]
        nc=b.shape[1]
        #print "self: ",(self.ndim, self.shape[0], self.size)
        #print "b: ",(nr, nc, b.size)
        if nr != self.shape[0]:
            raise ValueError("cannot dot different shaped arrays")
        # a temporary mpfr
        m.mpfr_init2(t, self.prec)
        m.mpfr_set_d(t, 0, m.MPFR_RNDD)
        # result
        r=ndarray(prec=self.prec, shape=(nc,))
        for i in range(nc):
            r[i]=0
            for j in range(nr):
                m.mpfr_mul(t, self.data[j],  b.data[j*nc + i], m.MPFR_RNDD)
                m.mpfr_add(r.data[i], r.data[i], t, m.MPFR_RNDD)
        return r

    def _mm_dot(self, ndarray b):
        """compute c=a.b with a, b, c matrices."""
        cdef m.mpfr_t t
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
        m.mpfr_init2(t, self.prec)
        m.mpfr_set_d(t, 0, m.MPFR_RNDD)
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
                    m.mpfr_mul(t, self.data[inc + k],  b.data[k*bc + j], m.MPFR_RNDD)
                    m.mpfr_add(r.data[i*bc + j], r.data[i*bc + j], t, m.MPFR_RNDD)
        return r


    def _dot(self, ndarray b):
        """compute c=a.b"""
        #error checking!!
        # if self.ld > 1 or b.ld > 1: throw an error
        if self.shape[self.ndim-1] != b.shape[0]:
            raise ValueError("cannot dot different shape arrays")
        cdef m.mpfr_t r,t
        cdef int i
        m.mpfr_init2(r, self.prec)
        m.mpfr_init2(t, self.prec)
        m.mpfr_set_d(r, 0, m.MPFR_RNDD)
        m.mpfr_set_d(t, 0, m.MPFR_RNDD)
        for i in range(b.size):
            m.mpfr_mul(t, self.data[i],  b.data[i], m.MPFR_RNDD)
            m.mpfr_add(r, r, t, m.MPFR_RNDD)
        m.mpfr_out_str(stdout, 10, 0, r, m.MPFR_RNDD)
        print


    ## Printing functions -- make recursive??

    def __str__(self):
        if self.ndim == 1:
            return self.vec_str()
        if self.ndim == 2:
            return self.mat_str()
    
    def __repr__(self): 
        if self.ndim == 1:
            return ("%i" % self.ndim) + self.vec_str()[:80]
        if self.ndim == 2:
            return ("%i" % self.ndim) + self[0].vec_str()[:80]
    
    def tolist(self):
        if self.ndim == 1:
            return [self[i].vec_str() for i in range(self.size)]
        else:
           raise TypeError("tolist function defined only for 1dim ndarrays")
    
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
            m.mpfr_sprintf(buf,"%.Re",self.data[0])
        else:
            buf[0]='['
            for i in range(self.size-1):
                n=m.mpfr_sprintf(buf + cn, "%.Re, ", self.data[i]) 
                cn +=n
            n=m.mpfr_sprintf(buf + cn, "%.Re]", self.data[self.size-1])
        return buf

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
                n=m.mpfr_sprintf(buf + cn, "%.Re", 
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
        return buf
    
    
    # utility functions

    def flatset(self, int i, value):
        """Set element of data array indexed _directly_ by i (i.e. forget shape)"""
        if isinstance(value, numbers.Real):
            m.mpfr_set_d(self.data[i], value, m.MPFR_RNDD)
        elif isinstance(value, str):
            m.mpfr_set_str(self.data[i], value, 10, m.MPFR_RNDD)
        else:
            raise TypeError("can only initialize mpfr from number or string")
            
    cdef flatset1(self, int i, m.mpfr_t value):
        m.mpfr_set(self.data[i], value, m.MPFR_RNDD)
     
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

    ####  -------------- OLD METHODS ------ DELETE ------

    def setd(self, int i, double value):
        """Setter for vector like array"""
        m.mpfr_set_d(self.data[i], value, m.MPFR_RNDD)

    #def __setitem__(self, int i, char* value):
    #    """String setter for vector like array"""
    #    m.mpfr_set_str(self.array[i], value, 10, m.MPFR_RNDD)

    #def __getitem__(self, int i):
    #    if i > self.size:
    #        return ''
    #    cdef char* buf = <char*> malloc(self.prec * sizeof(char))
    #    m.mpfr_sprintf(buf, "%.Re", self.array[i])
    #    return buf
    
    
    # Matrix functions
    
    def setd(self, int i, int j, double value):
        """Setter for vector like array"""
        m.mpfr_set_d(self.data[i*self.strides[0] + j], value, m.MPFR_RNDD)

    def sets(self, int i, int j, char* value):
        """String setter for vector like array"""
        m.mpfr_set_str(self.data[i*self.strides[0] + j], value, 10, m.MPFR_RNDD)

    def gets(self, int i, int j):
        if i * self.strides[0] + j > self.size:
            return ''
        cdef char* buf = <char*> malloc(self.prec * sizeof(char))
        m.mpfr_sprintf(buf, "%.Re", self.data[i * self.strides[0] + j])
        return buf


####################  END --- CLASS ------NDARRAY  #############################


