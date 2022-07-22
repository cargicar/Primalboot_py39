# cython: profile=True

##
# S. Rychkov, S. El-Showk, May-June 2013
# File contains prec_float and prec_interval classes : cython wrappers for mpfr_t and mpfi_t
# with overloaded arithmetic operations
##

# Edit by Carlos Cardona 2022

from libc.stdlib cimport malloc, free

import utils.stats as stats

def isinteger(a):
    return (isinstance(a,int) or isinstance(a,long))
    
def zeros(size,prec):
    return [ prec_float(0,prec=prec) for i in range(size)]

cpdef prec_float sqrt(prec_float x):
     cdef prec_float ret
     ret = prec_float (0, prec = <int>mpfr_get_prec (x.data)) 
     mpfr_sqrt(ret.data, x.data, MPFR_RNDD)
     return ret
     
cpdef double to_double(prec_float x):   
     return mpfr_get_d(x.data,MPFR_RNDD)


# used to unpickle a prec float
def make_prec_float(val, prec):
    return prec_float(val, prec=prec)

####################  BEGIN --- CLASS ------prec_float  #############################
cdef class prec_float:
    
    def __cinit__(self, value, int prec = 212):
        """initialize an mpfr float with precision prec, from integer, double or """
        # not working for some reason
        #if prec < 0:
        #    raise ValueError("precision must be specified")
        
        self.prec=prec
        mpfr_init2(self.data, self.prec)
        
        if isinstance (value, int):
            mpfr_set_si(self.data, value,MPFR_RNDD)
        elif isinstance(value, str):
           # mpfr_set_str(self.data, value, 10, MPFR_RNDD)
            bytes_str=value.encode('ascii') # Carlos edit 
            mpfr_set_str(self.data, bytes_str, 10, MPFR_RNDD)
        elif isinstance(value, float):
            mpfr_set_d(self.data, value,MPFR_RNDD)
        elif isinstance(value, prec_float):
            typed_value = <prec_float>value
            mpfr_set(self.data, typed_value.data, MPFR_RNDD)
        else:
            raise TypeError("Prec_float can be initialized from int, str, double, or prec_float")
        stats.pf_init +=1
        stats.pf_mem += sizeof(mpfr_t)
        
    def __dealloc__(self):
        mpfr_clear(self.data)
        # for some reason we can't access this this way
        stats.pf_del +=1
        stats.pf_free += sizeof(mpfr_t)
        
    ##  Overloaded functions

    def __neg__(prec_float self):
        """Sign flip"""
        ret=prec_float(0, prec = self.prec)
        mpfr_neg(ret.data, self.data, MPFR_RNDD)
        return ret    
    
    def __add__(prec_float self, prec_float other):
        cdef prec_float ret
        ret=prec_float(0, prec=self.prec)
        mpfr_add(ret.data, self.data,  other.data, MPFR_RNDD)
        stats.pf_adds += 1
        return ret   
 
    def __sub__(prec_float self, prec_float other):
        cdef prec_float ret
        ret=prec_float(0, prec=self.prec)
        mpfr_sub(ret.data, self.data,  other.data, MPFR_RNDD)
        return ret   

    def __mul__(prec_float self, prec_float other):
        cdef prec_float ret
        ret=prec_float(0, prec=self.prec)
        mpfr_mul(ret.data, self.data,  other.data, MPFR_RNDD)
        stats.pf_mults += 1
        return ret   
    


    # def __div__(prec_float self, prec_float other):
    def __truediv__(prec_float self, prec_float other):  # Carlos edit
        cdef prec_float ret
        ret=prec_float(0, prec=self.prec)
        mpfr_div(ret.data, self.data,  other.data, MPFR_RNDD)
        return ret   
        
    def __pow__(prec_float self, prec_float other, z):
        cdef prec_float ret
        ret=prec_float(0, prec=self.prec)
        mpfr_pow(ret.data, self.data, other.data, MPFR_RNDD)
        return ret
    
    def log(prec_float self):
        cdef prec_float ret
        ret=prec_float(0, prec=self.prec)
        mpfr_log(ret.data, self.data, MPFR_RNDD)
        return ret
    
    def exp(prec_float self):
        cdef prec_float ret
        ret=prec_float(0, prec=self.prec)
        mpfr_exp(ret.data, self.data, MPFR_RNDD)
        return ret
    
    def sqrt(prec_float self):
        cdef prec_float ret
        ret = prec_float(0, prec = self.prec) 
        mpfr_sqrt(ret.data, self.data, MPFR_RNDD)
        return ret
    
    def abs(prec_float self):
        """abs value"""
        ret=prec_float(0, prec = self.prec)
        mpfr_abs(ret.data, self.data, MPFR_RNDD)
        return ret

#    DON'T USE THESE - BUGGY
#    def __iadd__(prec_float self, prec_float other):
#        mpfr_add(self.data, self.data,  other.data, MPFR_RNDD)  

#    def __isub__(prec_float self, prec_float other):
#        mpfr_sub(self.data, self.data,  other.data, MPFR_RNDD)  

#    def __imul__(prec_float self, prec_float other):
#        mpfr_mul(self.data, self.data,  other.data, MPFR_RNDD)  

#    def __idiv__(prec_float self, prec_float other):
#        mpfr_div(self.data, self.data,  other.data, MPFR_RNDD)  

    def __richcmp__(prec_float self, prec_float value, int cmptype):
        """comparison
        op cmptype
        <	0   
        ==	2   
        >	4   
        <=	1   
        !=	3   
        >=	5"""
        cdef int cmp
       
        cmp = mpfr_cmp(self.data, value.data)
        
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

    ## Printing functions 
    def __str__(self):
        bs= 2*self.prec #  grossly overkill
        cdef char* buf = <char*> malloc(bs * sizeof(char))
        mpfr_sprintf(buf, "%.Re", self.data) 
        #return str(buf)                # Carlos edit
        return buf.decode('ascii')
    
    def __bytes__(self):
        bs= 2*self.prec #  grossly overkill
        cdef char* buf = <char*> malloc(bs * sizeof(char))
        mpfr_sprintf(buf, "%.Re", self.data) 
        return buf               # Carlos edit
    
    def __repr__(self):
        return self.__str__()
    
    def __reduce__(self):
        pstr=self.__str__()
        #ret=(prec_float.__init__, pstr, self.prec)
        return (make_prec_float, (pstr, self.prec))

    
    def bufrepr(self, buf):
        mpfr_sprintf(<char*>buf, "%.Re", self.data) 
        return buf
    
    
    def assign_prec_float(self, value):
        """assigns the already initialized precision float with data from another precision float;
        only works if the two precision floats have the same precision"""
        
        if not isinstance(value, prec_float):
            raise TypeError("assign_prec_float: attempted to assign value which is not prec_float")

        typed_value = <prec_float>value # cast static typed copy to be able to access private fields data and prec

        if not self.prec == typed_value.prec:
            raise TypeError("assign_prec_float: attempted to assign value which has different precision")
        
        mpfr_set(self.data, typed_value.data, MPFR_RNDD)
    
          

####################  BEGIN --- CLASS ------prec_interval  #############################

cdef class prec_interval:

   
    def __cinit__(self, v1, v2, int prec):
        """initialize an mpfi float with precision prec, from integer, double,
        string or prec_float.  If v1 is a string it should be either a number or
        in the form '[a,b]' and then v2 is ignored.  Otherwise v1 and v2 assumed to be of the same type."""
        self.prec=prec
        mpfi_init2(self.data, self.prec)
        
        if isinstance (v1, int):
            mpfi_interv_si(self.data,v1,v2)
        elif isinstance(v1, str):
            mpfi_set_str(self.data, v1, 10)
        elif isinstance(v1, float):
            mpfi_interv_d(self.data, v1, v2)
        elif isinstance(v1, prec_float):
            typed_value = <prec_float>v1
            typed_value2 = <prec_float>v2
            mpfi_interv_fr(self.data, typed_value.data, typed_value2.data)
        else:
            raise TypeError("Prec_interval can be initialized from int, str, double, or prec_float")
   
    def __dealloc__(self):
        mpfi_clear(self.data)
     
    ##  Overloaded functions

    def __neg__(prec_interval self):
        """Sign flip"""
        ret=prec_interval(0, 0, prec = self.prec)
        mpfi_neg(ret.data, self.data)
        return ret     
    
    def __add__(prec_interval self, prec_interval other):
        cdef prec_interval ret
        ret=prec_interval(0, 0, prec=self.prec)
        mpfi_add(ret.data, self.data,  other.data)
        return ret   

    def __sub__(prec_interval self, prec_interval other):
        cdef prec_interval ret
        ret=prec_interval(0, 0, prec=self.prec)
        mpfi_sub(ret.data, self.data,  other.data)
        return ret   

    def __mul__(prec_interval self, prec_interval other):
        cdef prec_interval ret
        ret=prec_interval(0, 0, prec=self.prec)
        mpfi_mul(ret.data, self.data,  other.data)
        return ret   

    def __div__(prec_interval self, prec_interval other):
        cdef prec_interval ret
        ret=prec_interval(0, 0, prec=self.prec)
        mpfi_div(ret.data, self.data,  other.data)
        return ret   
 
    def __pow__(prec_interval self, prec_interval other, z):
        cdef prec_interval ret
        ret=prec_interval(0, 0, prec=self.prec)
        mpfi_log(ret.data, self.data)
        mpfi_mul(ret.data, ret.data, other.data)
        mpfi_exp(ret.data, ret.data)
        return ret
    
    ## unpacks into a pair of prec_float
    def unpack(self):
        x0 = prec_float(0,prec=self.prec)
        x1 = prec_float(0,prec=self.prec)
        mpfi_get_left(x0.data, self.data)
        mpfi_get_right(x1.data, self.data)
        return (x0,x1)
    
    def length(self):
        x0 = prec_float(0,prec=self.prec)
        x1 = prec_float(0,prec=self.prec)
        mpfi_get_left(x0.data, self.data)
        mpfi_get_right(x1.data, self.data)
        return x1 - x0
    
    ## Printing functions 

    def __str__(self):
        return self.unpack().__str__()
    
    def __repr__(self):
        return self.__str__()
          
