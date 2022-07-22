# cython: profile=False

from libc.stdlib cimport malloc, free
from libc.stdio cimport stdout, sprintf, fprintf

cimport mpfr as m

import numbers
import numpy as np
import pickle

#cpdef prec_float prec_float_const(value, int prec = 212):
#     cdef prec_float ret
#     ret = prec_float(prec = prec)   
#     if isinstance (value, int):
#        m.mpfr_set_si(ret.data, value,m.MPFR_RNDD)
#     elif isinstance(value, str):
#        m.mpfr_set_str(ret.data, value, 10, m.MPFR_RNDD)
#     elif isinstance(value, float):
#        m.mpfr_set_d(ret.data, value,m.MPFR_RNDD)
#     else:
#        raise TypeError("Mpfr constants can be set from int, str or double only")   
#     return ret   

def isinteger(a):
    return isinstance(a,int) or isinstance(a,long)
    
def zeros(size,prec):
    return [ prec_float(0,prec=prec) for i in range(size)]

cpdef prec_float sqrt(prec_float x):
     cdef prec_float ret
     ret = prec_float (0, prec = <int>m.mpfr_get_prec (x.data))
     m.mpfr_sqrt(ret.data, x.data, m.MPFR_RNDD)
     return ret
    
cpdef double to_double(prec_float x):
     return m.mpfr_get_d(x.data,m.MPFR_RNDD)
        
####################  BEGIN --- CLASS ------prec_float  #############################

cdef class prec_float:

    # vars same as numpy.prec_float
    cdef m.mpfr_t data
    
    cdef public int prec
   
    def __cinit__(self, value, int prec):
        """initialize an mpfr float with precision prec, from integer, double or """
        self.prec=prec
        m.mpfr_init2(self.data, self.prec)
        
        if isinstance (value, int):
            m.mpfr_set_si(self.data, value,m.MPFR_RNDD)
        elif isinstance(value, str):
            m.mpfr_set_str(self.data, value, 10, m.MPFR_RNDD)
        elif isinstance(value, float):
            m.mpfr_set_d(self.data, value,m.MPFR_RNDD)
        elif isinstance(value, prec_float):
            typed_value = <prec_float>value
            m.mpfr_set(self.data, typed_value.data, m.MPFR_RNDD)
        else:
            raise TypeError("Prec_float can be initialized from int, str, double, or prec_float")
   
    def __dealloc__(self):
        m.mpfr_clear(self.data)
        
    ##  Overloaded functions

    def __neg__(prec_float self):
        """Sign flip"""
        ret=prec_float(0, prec = self.prec)
        m.mpfr_neg(ret.data, self.data, m.MPFR_RNDD)
        return ret    
    
#    def __add__(prec_float self, prec_float other):
#        ret=prec_float(prec=self.prec)
#        m.mpfr_add(ret.data, self.data,  other.data, m.MPFR_RNDD)
#        return ret   
 
    def __add__(self, other):
        # we are in cython so first we need to check who is true self
        if isinstance (self, prec_float):
            typed_self = <prec_float> self
            ret=prec_float(0, prec=typed_self.prec)
            if isinstance (other, prec_float):
                typed_other = <prec_float>other
                m.mpfr_add(ret.data, typed_self.data,typed_other.data, m.MPFR_RNDD)
            elif isinstance (other, int):
                m.mpfr_add_si(ret.data, typed_self.data, other, m.MPFR_RNDD)
            else: 
                raise TypeError("unsupported type in __add__"+str(type(self))+","+str(type(other))) 
        elif isinstance (self, int): # isinstance (other, prec_float) implicit
            typed_other = <prec_float>other
            ret=prec_float(0, prec=typed_other.prec)
            m.mpfr_add_si(ret.data, typed_other.data, self, m.MPFR_RNDD)
        else: 
                raise TypeError("unsupported type in __add__"+str(type(self))+","+str(type(other)))
        return ret   
 
#    def __sub__(prec_float self, prec_float other):
#        ret=prec_float(prec=self.prec)
#        m.mpfr_sub(ret.data, self.data,  other.data, m.MPFR_RNDD)
#        return ret   

    def __sub__(self, other):
        # we are in cython so first we need to check who is true self
        if isinstance (self, prec_float):
            typed_self = <prec_float> self
            ret=prec_float(0, prec=typed_self.prec)
            if isinstance (other, prec_float):
                typed_other = <prec_float>other
                m.mpfr_sub(ret.data, typed_self.data,typed_other.data, m.MPFR_RNDD)
            elif isinstance (other, int):
                m.mpfr_sub_si(ret.data, typed_self.data, other, m.MPFR_RNDD)
            else: 
                raise TypeError("unsupported type in __sub__"+str(type(self))+","+str(type(other))) 
        elif isinstance (self, int): # isinstance (other, prec_float) implicit
            typed_other = <prec_float>other
            ret=prec_float(0, prec=typed_other.prec)
            m.mpfr_si_sub(ret.data, self, typed_other.data, m.MPFR_RNDD)
        else: 
                raise TypeError("unsupported type in __sub__"+str(type(self))+","+str(type(other)))
        return ret  

#    def __mul__(prec_float self, prec_float other):
#        ret=prec_float(prec=self.prec)
#        m.mpfr_mul(ret.data, self.data,  other.data, m.MPFR_RNDD)
#        return ret   

    def __mul__(self, other):
        # we are in cython so first we need to check who is true self
        if isinstance (self, prec_float):
            typed_self = <prec_float> self
            ret=prec_float(0, prec=typed_self.prec)
            if isinstance (other, prec_float):
                typed_other = <prec_float>other
                m.mpfr_mul(ret.data, typed_self.data,typed_other.data, m.MPFR_RNDD)
            elif isinstance (other, int):
                m.mpfr_mul_si(ret.data, typed_self.data, other, m.MPFR_RNDD)
            else: 
                raise TypeError("unsupported type in __mul__:"+str(type(self))+","+str(type(other))) 
        elif isinstance (self, int): # isinstance (other, prec_float) implicit
            typed_other = <prec_float>other
            ret=prec_float(0, prec=typed_other.prec)
            m.mpfr_mul_si(ret.data, typed_other.data, self, m.MPFR_RNDD)
        else: 
                raise TypeError("unsupported type in __mul__:"+str(type(self))+","+str(type(other)))
        return ret  

    def __div__(self, other):
        # we are in cython so first we need to check who is true self
        if isinstance (self, prec_float):
            typed_self = <prec_float> self
            ret=prec_float(0, prec=typed_self.prec)
            if isinstance (other, prec_float):
                typed_other = <prec_float>other
                m.mpfr_div(ret.data, typed_self.data,typed_other.data, m.MPFR_RNDD)
            elif isinteger (other):
                m.mpfr_div_si(ret.data, typed_self.data, other, m.MPFR_RNDD)
            else: 
                raise TypeError("unsupported type in __div__"+str(type(self))+","+str(type(other))) 
        elif isinteger (self): # isinstance (other, prec_float) implicit
            typed_other = <prec_float>other
            ret=prec_float(0, prec=typed_other.prec)
            m.mpfr_si_div(ret.data, self, typed_other.data, m.MPFR_RNDD)
        else: 
                raise TypeError("unsupported type in __div__"+str(type(self))+","+str(type(other)))
        return ret  

#    def __div__(prec_float self, prec_float other):
#        ret=prec_float(prec=self.prec)
#        m.mpfr_div(ret.data, self.data,  other.data, m.MPFR_RNDD)
#        return ret   
        
    def __pow__(prec_float self, other, z):
        ret=prec_float(0, prec=self.prec)
        if isinstance(other,int):
            m.mpfr_pow_si(ret.data, self.data, other, m.MPFR_RNDD)
        elif isinstance(other,prec_float):
            typed_other = <prec_float>other
            m.mpfr_pow(ret.data, self.data, typed_other.data, m.MPFR_RNDD)
        return ret   

#    DON'T USE THESE - BUGGY
#    def __iadd__(prec_float self, prec_float other):
#        m.mpfr_add(self.data, self.data,  other.data, m.MPFR_RNDD)  

#    def __isub__(prec_float self, prec_float other):
#        m.mpfr_sub(self.data, self.data,  other.data, m.MPFR_RNDD)  

#    def __imul__(prec_float self, prec_float other):
#        m.mpfr_mul(self.data, self.data,  other.data, m.MPFR_RNDD)  

#    def __idiv__(prec_float self, prec_float other):
#        m.mpfr_div(self.data, self.data,  other.data, m.MPFR_RNDD)  

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
       
        cmp = m.mpfr_cmp(self.data, value.data)
        
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
        m.mpfr_sprintf(buf, "%.Re", self.data) 
        return buf
    
    def __repr__(self):
        return self.__str__()
          
 
####################  END --- CLASS ------prec_float  #############################


