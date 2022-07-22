import setpath

import numpy as np
import os
import time

import config


from version import  Spectrum, lp_problem, NP_PREC, mpfr_array, prec_float

description='a test run'
dim=3
opedim=dim
opespin=2
#nmax=9
# nmax = 6 # Carlos edit # Move nmax and mmax to config.py, because i will be changing those two more often.
# mmax=1
# d=2
#s0=0.1235
#s1=0.1257
s0=0.5179
s1=0.5187
sdelta=0.00005
threshold=1e-12

###################### block below commented out Carlos edit
# if True:
#     s0=0.115
#     s1=0.135
#     sdelta=0.002
#     s0=0.5179
#     #s1=0.518501
#     s1=0.517901
#     sdelta=0.00015
#     #s1=0.518901
#     #s1=0.5186
#     #sdelta=0.0001
#     #sdelta=0.00005
#     #siglist=np.arange(s0,s1,sdelta)
#     nmax=6
#     nmax=3
#     threshold=1e-8
#     threshold=1e-5

# a basic health check (20 points)
# if False:
#     #s0=0.115
#     #s1=0.135
#     #sdelta=0.002
#     s0=0.5160
#     s1=0.5200
#     sdelta=0.0002
#     nmax=9
#     threshold=1e-8

###################### block below commented out Carlos edit
#epsmin=0.81 #d=2
#epsmin=0.91 #d=3
epsmin=1 #d=3
# don't give an nice bound:
#epsmin=0.51 #d=3
#epsmin=0.61 #d=3


# -----------------------
# Cluster params
# -----------------------

queue='8nm'
queue='1nh'
queue='8nh'
bsuboutput='/output.bsub'
clusterhome='/afs/cern.ch/user/s/sheer/'
clusterhome='~/'
clusterpybin=clusterhome+'/Utils/local/bin/python2.7'
clusterCWD=clusterhome+'Dev/python/cern_new/py'
# disables loading ope_bound_analysis which requires libraries not on the
# cluster
clusterMode=True

# -------------------------
# Batch dir/file structure
# -------------------------

# work dir where all batches will be kept
workdir='../workdir/'

# files in jobdir
sigfile='/sig.p'
resfile='/res.p'
runfile='/running'
output='/output'

# for running the reaper
pybin='python2.7'
reaper='reaper.py'



