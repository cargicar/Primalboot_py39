cdef extern from "mpfr_funcs.h":
    double* printarray(int len, double* arr)
    void testmpfr()

cimport mpfr as m

from libc.stdlib cimport malloc
from libc.stdio cimport stdout, FILE, printf

import numpy as np
cimport numpy as np

import mpfr_array as mpa

import random as r
import datetime

#
# See this url:
#
# http://stackoverflow.com/questions/3046305/simple-wrapping-of-c-code-with-cython
#



def fpy(N, np.ndarray[np.double_t, ndim=1] A):
    printarray(len(A), <double*> A.data)

def cythontest():
    a=np.arange(3, dtype=np.float64)
    fpy(len(a), a)
    print a

def testmpfr():
    # test mpfr
    testmpfr()
    cdef m.mpfr_t* s 
    s=<m.mpfr_t*>malloc(4*sizeof(m.mpfr_t))
    m.mpfr_init2(s[0], 212)
    m.mpfr_init2(s[1], 212)
    m.mpfr_init2(s[2], 212)
    m.mpfr_init2(s[3], 212)
    m.mpfr_set_d(s[0], 1.0, m.MPFR_RNDD)
    m.mpfr_set_d(s[1], 3.0, m.MPFR_RNDD)
    m.mpfr_div(s[2], s[0], s[1], m.MPFR_RNDD)
    m.mpfr_out_str(stdout, 10, 0, s[2], m.MPFR_RNDD);
    print
    p="1.1e-3"
    cdef char* sp=p
    printf("string is %s\n", sp)
    m.mpfr_set_str(s[0], sp, 10, m.MPFR_RNDD)
    m.mpfr_out_str(stdout, 10, 0, s[0], m.MPFR_RNDD);
    print

def test_mpfr_array():
    ma=mpa.ndarray(prec=40, shape=(10, 5, 30))
    #print ma.block_size
    ma.bs()
    mb=mpa.ndarray(prec=40, shape=(4,))
    mb[0]="1.2323"
    mb[1]="89.567"
    mb.setd(2, 9284934.23423)
    print mb
    mc=mpa.ndarray(prec=200, shape=(4,))
    mc[0]="2.234241142524353534543543454e10"
    mc[2]="334.928434"
    print mc
    print mc[0]
    md=mpa.ndarray(prec=20, shape=(4,4))
    for i in range(4):
        for j in range(4):
            md.sets(i,j, str(r.random()))
    print "md:"
    print md

def test_vec_vec_dot(n, pr=False):
    ma=mpa.ndarray(prec=200, shape=(n,))
    mb=mpa.ndarray(prec=200, shape=(n,))
    na=np.ndarray(shape=(n,))
    nb=np.ndarray(shape=(n,))
    for i in range(n):
        na[i]=r.random()
        nb[i]=r.random()
        ma[i]=na[i]
        mb[i]=nb[i]
    if pr:
        print ma
        print mb
    t0 = datetime.datetime.now()
    mr=mpa.dot(ma, mb)
    t1 = datetime.datetime.now()
    mtm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    t0 = datetime.datetime.now()
    nr=np.dot(na, nb)
    t1 = datetime.datetime.now()
    ntm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    if pr:
        print 'mpfr result: ', mr
        print 'np result: ', nr
    print '(numvec, mpfrtime, nptime, ratio): ', (n,mtm,ntm,mtm/ntm)


def test_vec_mat_dot(m,n):
    mv=mpa.ndarray(prec=200, shape=(n,))
    mm=mpa.ndarray(prec=200, shape=(m,n))
    nv=np.ndarray(shape=(n,))
    nm=np.ndarray(shape=(m,n))
    for i in range(n):
        nv[i]=r.random()
        mv[i]=nv[i]
    for i in range(m):
        for j in range(n):
            nm[i][j]=r.random()
            mm[i][j]=nm[i][j]
    t0 = datetime.datetime.now()
    mr=mpa.dot(mm, mv)
    t1 = datetime.datetime.now()
    mtm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    t0 = datetime.datetime.now()
    nr=np.dot(nm, nv)
    t1 = datetime.datetime.now()
    ntm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    print 'mpfr result: ', mr[0], mr[1]
    print 'np result: ', nr[0:2]
    print '(numvec, mpfrtime, nptime, ratio): ', (n,mtm,ntm,mtm/ntm)


def test_mpfr_subarray():
    md=mpa.ndarray(prec=20, shape=(6,6))
    for i in range(6):
        for j in range(6):
            md.sets(i,j, str(r.random()))
    print "parent:"
    print md
    ms=mpa.ndarray(parent=md, p_ind=2)
    print "child:"
    print ms
    ms1=mpa.ndarray(parent=ms, p_ind=3)
    print "subchild:"
    print ms1
    print "accessing via m[i][j]"
    print md[2]
    print md[2][3]
    print "should raise index error"
    try:
        print md[10]
    except IndexError:
        print "caught index error"
    print "setting via m[1][1] = 1.3"
    md[1][1]=1.3
    print md
    print "setting via m[1][1] = '1.23456789'"
    md[1][1]='1.23456789'
    print md
    print "setting via m[0][0] = 19"
    md[0][0]=19
    print md
    print "should raise type exception"
    try:
        md[1]=3
    except TypeError as e:
        print "caught type error:", e


def test_mpfr_arithmatic():
    v1=mpa.ndarray(prec=20, shape=(6,))
    for i in range(6):
        v1[i]=r.random()
    m1=mpa.ndarray(prec=20, shape=(6,6))
    for i in range(6):
        for j in range(6):
            m1[i][j]=r.random()
    print "v1:"
    print v1
    print "v1+v1"
    print v1 + v1
    print "v1*v1"
    print v1 * v1
    print "v1*v1 -v1"
    print v1 * v1 -v1
    print "(v1*v1 -v1)/v1"
    print (v1 * v1 -v1)/v1
    print "m1:"
    print m1
    print "m1+m1"
    print m1 + m1
    print "should raise value exception"
    try:
        m1/v1
    except ValueError as e:
        print "caught value error:", e
    print "m1[2]/v1"
    print m1[2]/v1

def test_varying_dot(int rng, int step, int m, int n):
    cdef i, j, ni
    for i in range(1, rng, step):
        ni = i * n
        mv=mpa.ndarray(prec=200, shape=(ni,))
        mm=mpa.ndarray(prec=200, shape=(m,ni))
        nv=np.ndarray(shape=(ni,))
        nm=np.ndarray(shape=(m,ni))
        for i in range(ni):
            nv[i]=r.random()
            mv[i]=nv[i]
        for i in range(m):
            for j in range(ni):
                nm[i][j]=r.random()
                mm[i][j]=nm[i][j]
        t0 = datetime.datetime.now()
        mr=mpa.dot(mm, mv)
        t1 = datetime.datetime.now()
        mtm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
        t0 = datetime.datetime.now()
        nr=np.dot(nm, nv)
        t1 = datetime.datetime.now()
        ntm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
        #print '(m, n, total ops): ', (m, ni, m*ni*ni)
        #print '(mpfrtime, nptime, ratio): ', (mtm,ntm,mtm/ntm)
        print '(mpfrtime/i^2, nptime/i^2, ratio): ',\
            (mtm/(i*i),ntm/(i*i),mtm/ntm)

def test_iter(n):
    ma=mpa.ndarray(prec=200, shape=(n,))
    for i in range(n):
        ma[i]=r.random()
    for x in ma:
        print x

def test_array_set():
    #print mpa.getshape([5,2,3])
    #print mpa.getshape([1,[2,2],3])
    #print mpa.getshape([[1,4],[2,[2,2]],[3,0]])
    #print mpa.getshape([[1,4],[2,5],[3,0]])
    marr=mpa.array([[1,4],[2,5],[3,0]])
    print marr
    marr=mpa.array([[1.234,4.23],[2.23423,5.354645672],[3.2435,0.85e10]])
    print marr
    marr=mpa.array([
        ["23459.3453453465352556456735454535234243e23423",
         "65456.2563465436445645645645e-3434"],
        ["905436845645646543654564.5643656346456456454556476345e-2434534",
            "8945674.456437436456463737329290534434e43535352"]])
    print marr
    print marr * marr
    print marr / marr

def test_set_to_array():
    marr=mpa.array([[1.234,4.23, 8.2],[2.23423,5.354645672,
        2323],[3.2435,0.85e10, -1]], prec=40)
    print marr
    print 'setting m[1][2]=3'
    marr[1][2]=3
    print marr
    vec=mpa.array([9234.2342, 0.4353434, 834534.454], prec=40)
    print 'setting m[1] to:', vec
    marr[1]=vec
    print marr
    marr=mpa.array([[[1,2],[3,4]],[[5,6],[7,8]],[[9,10],[11,12]]])
    print 'Three tensor (no print method yet so printing row by row)'
    print marr[0]
    print marr[1]
    print marr[2]
    vec = mpa.array([[25,26],[27,28]])
    print 'replacing m[1] with (matrix): '
    print vec
    marr[1]=vec
    print 'new three tensor'
    print marr[0]
    print marr[1]
    print marr[2]
    vec = mpa.array([39,310])
    print 'replacing m[2][0] with: ', vec
    marr[2][0]=vec
    print 'new three tensor'
    print marr[0]
    print marr[1]
    print marr[2]


def test_inv_perf(n):
    ones=np.diag(np.ones(n))
    npM=np.random.rand(n*n).reshape([n,n])
    mM=mpa.array(npM)
    mones=mpa.array(ones)
    print npM.shape
    t0 = datetime.datetime.now()
    ninv=np.linalg.inv(npM)
    t1 = datetime.datetime.now()
    nptm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    print 'numpy inverse time for',n,'x',n,' matrix:',nptm
    t0 = datetime.datetime.now()
    minv=mpa.inv(mM)
    t1 = datetime.datetime.now()
    mtm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    print 'mpfr inverse time for',n,'x',n,' matrix:',mtm,'(', mtm/nptm,\
            ' times slower)'
    print "here are some entries from m*minv - 1"
    nz=np.dot(npM, ninv) - ones
    mz=mpa.dot(mM, minv) - mones
    for l in range(10):
        i=r.randint(0,n-1)
        j=r.randint(0,n-1)
        print 'np:',nz[i][j],' mpfr:',mz[i][j]

def test_mat_inv():
    print "here's a matrix"
    ones=mpa.array([[1, 0, 0], [0,1,0], [0,0,1]])
    marr=mpa.array([[1.234,4.23, 8.2],[2.23423,5.354645672,
        2323],[3.2435,0.85e10, -1]])
    print marr
    minv=mpa.inv(marr)
    print "here's its inverse"
    print minv
    print "here m*minv - 1"
    print mpa.dot(marr, minv) - ones
    print 'Now some performance testing'
    test_inv_perf(66)
    test_inv_perf(120)





print "test mpfr_array class"
test_mpfr_array()
print "dot product of length 4 vecs"
test_vec_vec_dot(4, True)
print
print "timing for dot product length 100k"
test_vec_vec_dot(100000)
print
print "testing subarrays:"
test_mpfr_subarray()
print
print "testing array arithmatic:"
test_mpfr_arithmatic()
print
print "dot product of 10x4 mat with 4 vec"
test_vec_mat_dot(10,4)
print
print "dot product of 200x100 mat with 100 vec"
test_vec_mat_dot(200, 100)
print
print "dot product of 200x100 mat with 100 vec"
#test_vec_mat_dot(4000, 66)
print
print "varying dot product size"
#test_varying_dot(20, 2, 2000, 50)
print
print "varying dot product size"
#test_varying_dot(20, 2, 4000, 20)
print
print "testing iterator"
test_iter(5)
print
print "testing array(list) function"
test_array_set()
print
print "testing set m[i]=mpfr_array"
test_set_to_array()
print
print "testing matrix inverstion"
test_mat_inv()
print
