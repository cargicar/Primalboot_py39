from __future__ import print_function
from __future__ import with_statement

import solver.v1.cb_eval.cb_eval_prec as CBEval
import solver.v1.prec_float.prec_float as PF
from solver.v1.prec_float.prec_float import prec_float as pf
import numpy as np

import config


def save_float ( x, FILE):
    print("%.18e" % x, file = FILE)
def save_int ( x, FILE):
    print("%i" % x, file = FILE)
def save_str ( x, FILE):
    print(x, file = FILE)
def save_pf ( x, FILE):
    print(x.__str__() +" %i" % x.prec, file = FILE)

def read_float( FILE):
    return float(FILE.readline())
def read_int(FILE):
    return int(FILE.readline())
def read_str(FILE):
    return FILE.readline()
def read_pf(FILE):
    line = FILE.readline()
    string,prec=line.split(" ")
    return pf(string,prec = int(prec))

############################################################################
#  
# CBTABLE
# -------
#
# The CB_Table class represents a conformal block table 
#
#
#############################################################################

def make_dlist(eps, lmax, dim_all_max, step, dim_scalar_min ):
    """ uniform step """
    dlist = [[] for l in range(0,lmax+1)]
    for l in range(0,lmax+1,2):
        if l==0:
            dmin = dim_scalar_min
        else:
            dmin = l + 2*eps
        dlist[l] = np.arange(dmin,dim_all_max+step,step).tolist()
    return dlist

def make_dlist1(eps, lmax, dim_all_max, steps, dim_scalar_min ):
    """ steps is a dict ==> {l: [[step1, length1],[ step2, length2],...}"""
    dlist = [[] for l in range(0,lmax+1)]

    for l in range(0,lmax+1,2):
        if l==0:
            dmin = dim_scalar_min + config.fudge
        else:
            dmin = l + 2*eps + config.fudge
            #dmin = l + 2*eps
        dlist[l] = [dmin]
        for step in steps[l]:
            dimlast = dlist[l][-1]
            dlist[l] += np.arange(dimlast+step[0],min(dimlast+step[1],dim_all_max)+step[0],step[0]).tolist()
            if dlist[l][-1]>=dim_all_max:
                break
    return dlist    
    
class CB_Table:
    """ help by SR 30.03.13:
    contains the table of unnormalized conformal block derivative vectors
    """
    def __init__(self, eps=0, dlist=[], nmax=0, mmax=0, prec = 212, FILE = ""):
        if FILE=="":
            self.eps = float(eps)
            self.nmax = nmax
            self.mmax = mmax   
            self.prec = prec
            self.dlist = dlist
            print("computing table...")
            self.compute()
        else: 
            print("loading table from file..."+FILE)
            self.read(FILE)    
        
    def compute(self):
        self.table = [ []  for l in range(len(self.dlist))];
        
        for l in range(0,len(self.dlist),2):
            print ("l=",l)
            self.table[l] = [CBEval.cb_ders(self.eps, l, opdim, 
                        self.nmax, self.mmax, prec=self.prec) 
                                for opdim in self.dlist[l]]

    def save(self, filename = ""):
        if filename == "":
            filename = "CBeps"+str(self.eps)+"n"+str(self.nmax)+"m"+str(self.mmax)+".txt"
        with open(filename,"w") as FILE:
            save_float(self.eps,FILE)
            save_int(self.nmax, FILE)
            save_int(self.mmax, FILE)
            save_int(self.prec, FILE)
            save_int(len(self.dlist), FILE)
            for l in range(len(self.dlist)):
                save_int( len(self.dlist[l]),  FILE)
                for opdim in self.dlist[l]:
                    save_float( opdim, FILE)
                    
            for l in range(len(self.dlist)):
                for h in self.table[l]:
                    for n in range(self.nmax+1):
                        for m in range(2*(self.nmax-n)+self.mmax+1):
                            save_pf(h[n][m], FILE)
        return filename
                    
    def read(self, filename):
        with open(filename,"r") as FILE:
            self.eps= read_float(FILE)
            self.nmax = read_int(FILE)
            self.mmax = read_int(FILE)
            self.prec = read_int(FILE)
            len_dlist = read_int(FILE)
            self.dlist = [[] for l in range(len_dlist)]
            for l in range(len(self.dlist)):
                len_l = read_int( FILE)
                self.dlist[l] = [read_float( FILE) for i in range(len_l)]
            
            self.table = [ []  for l in range(len(self.dlist))];
            for l in range(len(self.dlist)):
                #print ("l=",l)       
                self.table[l] = [
                    [[read_pf(FILE) for m in range(2*(self.nmax-n)+self.mmax+1) ]
                     for n in range(self.nmax+1)]
                       for opdim in self.dlist[l]]
                       

############################################################################
#  
# SIGMA_CBTABLE
# -------------
#
# The Sigma_Table class represents a NORMALIZED CB table convolved with a particular
# value of sigma.  So a new instance is needed for each value of sigma.
#
#
#############################################################################

class Sigma_Table:
   
    def __init__(self,ds=0.125, cbtab = [], FILE=""):
        # delta_simga (dimension of external scalar)
        if FILE == "":
            self.ds = ds
            # the cbtable read from disk
            self.cbtab = cbtab
            self.eps = cbtab.eps
            self.nmax = cbtab.nmax
            self.mmax = cbtab.mmax   
            self.prec = cbtab.prec
            self.dlist = cbtab.dlist
            
            coeff_tensor = CBEval.coeff_tensor(ds, self.nmax, self.mmax, self.prec)
            self.table = [[ CBEval.contract(coeff_tensor, h, self.nmax, self.mmax) 
                for h in self.cbtab.table[l]] for l in range(len(self.dlist))]
            self.CBlen = len(self.table[0][0])
            
            h0 = [[1 if m==0 and n==0 else 0 for m in range(2*(self.nmax-n)+self.mmax+1) ] for n in range(self.nmax+1)]
            self.unitCB = CBEval.contract(coeff_tensor, h0, self.nmax, self.mmax) 
        else: 
            self.read(FILE)
       
    def save(self, filename = ""):
        if filename == "":
            filename = "SCBeps"+str(self.eps)+"ds"+str(self.ds)+"n"+str(self.nmax)+"m"+str(self.mmax)+".txt"
        with open(filename,"w") as FILE:
            save_float(self.eps,FILE)
            save_float(self.ds,FILE)
            save_int(self.nmax, FILE)
            save_int(self.mmax, FILE)
            save_int(self.prec, FILE)
            save_int(self.CBlen, FILE)
            save_int(len(self.dlist), FILE)
            
            for l in range(len(self.dlist)):
                save_int( len(self.dlist[l]),  FILE)
                for opdim in self.dlist[l]:
                    save_float( opdim, FILE)
                    
            for l in range(len(self.dlist)):
                for h in self.table[l]:
                    for i in range(self.CBlen):
                        save_float(h[i], FILE)
            
            for i in range(self.CBlen):
                        save_float(self.unitCB[i], FILE)
                                
        return filename
                    
    def read(self, filename):
        with open(filename,"r") as FILE:
            self.eps= read_float(FILE)
            self.ds= read_float(FILE)
            self.nmax = read_int(FILE)
            self.mmax = read_int(FILE)
            self.prec = read_int(FILE)
            self.CBlen = read_int(FILE)
            len_dlist = read_int(FILE)
            self.dlist = [[] for l in range(len_dlist)]
            for l in range(len(self.dlist)):
                len_l = read_int( FILE)
                self.dlist[l] = [read_float( FILE) for i in range(len_l)]
            
            self.table = [ []  for l in range(len(self.dlist))];
            for l in range(len(self.dlist)):
                self.table[l] = [[read_float(FILE) for i in range(self.CBlen)]
                       for opdim in self.dlist[l]]
           
            self.unitCB = [read_float(FILE) for i in range(self.CBlen)]

# ------------- END SIGMA_CBTABLE
            
