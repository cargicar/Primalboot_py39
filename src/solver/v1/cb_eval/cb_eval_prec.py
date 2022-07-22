import numpy as np

from math import factorial as fac

import solver.v1.prec_float.prec_float as PF
from solver.v1.prec_float.prec_float import prec_float as pf

# TAG i refers to Mathematica notebook

def scalar_seed (eps0, opdim0, a = 0, b = 0, prec = 212):
    """ returns z-derivatives 0,1,2 along the diagonal at z=1/2 up to by summing 
    power series;
    it is implicit that eps """
    eps = pf(eps0,prec)
    opdim = pf(opdim0,prec)    
    root2 = PF.sqrt(pf(2,prec));
    r0 = 3-2*root2
    mmax= 2
    der_r = PF.zeros(mmax+1,prec) #for accumulating the first derivatives in the r coordinate
    deltader = PF.zeros(mmax+1,prec) #increments for der_r
    if (a == 0 and b == 0):
        # work with expansion coefficients in the variable r^2
        k = 0
        r02k = pf(1,prec) # r0^2k
        while k<1000:
            # c0, c1, c2 are correspodningly c[m], c[m-1] and c[m-2]
            if (k == 0):
                c0 = (r0)**opdim # reduced from (4*r0)**opdim
                c1 = 0
            else:
                # nb TAG 1  
                denom = k*(-1 - eps + k + opdim)*(-1 + 2*k + opdim)
                g1 = pf((-2 + 2*k + opdim)*(2 - eps - opdim + eps*opdim + 2*k*(-2 + k + opdim)),prec)/denom
                g2 = -pf((-1 + eps + k)*(-2 + k + opdim)*(-3 + 2*k + opdim),prec)/denom 
                c2 = c1
                c1 = c0
                c0 = g1 * c1 + g2 * c2
            
            #compute the increment
            deltader[0] = c0*r02k;
            for i in range(1,mmax+1):
                deltader[i] = deltader[i-1]*(opdim + 2*k - i+1)/r0
            
            #check if increment is still significant
            if k > 10 and deltader[mmax]/der_r[mmax] < pf(2**(-prec-3),prec):
                break # exit while cycle
            
            # accumulate    
            for i in range(mmax+1):
                der_r[i] = der_r[i] + deltader[i]
            r02k = r02k*r0*r0
            k+=1
        else:
            raise ValueError("k=1000 is reached")
            
        der_z = PF.zeros(mmax+1,prec)
        der_z[0] = der_r[0]
        # see nb TAG 2
        der_z[1] = (-8 + 6*root2)*der_r[1]
        der_z[2] = (32 - 22*root2)*der_r[1] + (136 - 96*root2)*der_r[2] 
  
    else:
        raise ValueError("Nonzero a,b not yet supported")
    return  der_z        
    
    
def spin_seed (eps0, l, opdim0, a = 0, b = 0, prec = 212 ):
    """ returns z-derivatives 0,1,2,3 along the diagonal at z=1/2 up to by summing 
    power series """
    eps = pf(eps0,prec)
    opdim = pf(opdim0,prec)    
    root2 = PF.sqrt(pf(2,prec));
    r0 = 3-2*root2
    mmax= 3
    der_r = PF.zeros(mmax+1,prec) #for accumulating the first derivatives in the r coordinate
    deltader = PF.zeros(mmax+1,prec) #increments for der_r
    if (a == 0 and b == 0):
        # work with expansion coefficients in the variable r^2
        k = 0
        r02k = pf(1,prec) # r0^2k
        while k<1000:
            # c0, c1, c2 are correspodningly c[m], c[m-1] and c[m-2]
            if k == 0:
                c0 = (r0)**opdim # reduced from (4*r0)**opdim
            elif k == 1:
                # nb TAG 4: need to compute k=1 separately to avoid problems near unitarity bound
                c1 = c0
                c0 = c0 * (2*eps*l*(eps*(-2 + opdim) + opdim) + l**2*(eps*(-2 + opdim) + opdim) - opdim*(1 - 2*eps + opdim)*(eps*(-1 + opdim) + opdim))/((-1 + 2*eps + l - opdim)*(-eps + opdim)*(1 + l + opdim))
                c2 = 0
            else:
                # nb TAG 3  
                denom = pf(k*(-1 - eps + k + opdim)*(-1 - 2*eps + 2*k - l + opdim)*(-1 + 2*k + l + opdim),prec)
                
                g1 = 12*k**4 + 8*k**3*(-7 - eps + 3*opdim) + (-2 + opdim)*(-17 + 2*l**2 + (11 - 2*opdim)*opdim) - 2*eps**2*(5 + l*(-1 + opdim) + (-4 + opdim)*opdim) - k**2*(-107 - 22*eps + 4*eps**2 + 3*l**2 + 3*(26 - 5*opdim)*opdim + 6*eps*(l + opdim)) + eps*(9 + 4*l*(-2 + opdim) - l**2*(-1 + opdim) + opdim*(-1 + (-3 + opdim)*opdim)) + k*(-97 - 23*eps + 2*eps**2*(7 - 3*opdim) + 2*eps*l*(7 + eps - 3*opdim) + l**2*(7 + eps - 3*opdim) + eps*opdim*(6 + opdim) + opdim*(93 + opdim*(-29 + 3*opdim)))
                
                g2 = -166 + 78*eps + 4*eps**2 - 12*k**4 + 4*eps*l*(5 + eps*(-2 + opdim) - 2*opdim) + 2*l**2*(5 + eps*(-2 + opdim) - 2*opdim) - 8*k**3*(-11 + eps + 3*opdim) + k**2*(-251 + 50*eps + 4*eps**2 + 3*l**2 + 6*eps*(l - 3*opdim) + 3*(42 - 5*opdim)*opdim) - 2*eps*opdim*(35 + (-10 + opdim)*opdim) + opdim*(143 + opdim*(-41 + 4*opdim)) + k*(329 - 107*eps + 2*eps**2*(-5 + opdim) + eps*(70 - 11*opdim)*opdim + 2*eps*l*(-11 + eps + 3*opdim) + l**2*(-11 + eps + 3*opdim) + opdim*(-229 + (49 - 3*opdim)*opdim))
                
                g3 = (-2 + eps + k)*(-3 + k + opdim)*(-5 + 2*k - l + opdim)*(-5 + 2*eps + 2*k + l + opdim)
               
                g1=g1/denom
                g2=g2/denom
                g3=g3/denom
                c3 = c2
                c2 = c1
                c1 = c0
                c0 = g1 * c1 + g2 * c2 + g3 * c3
                
            #compute the increment
            deltader[0] = c0*r02k;
            for i in range(1,mmax+1):
                deltader[i] = deltader[i-1]*(opdim + 2*k - i+1)/r0
            
            #check if increment is still significant
            if k > 10 and deltader[mmax]/der_r[mmax] < pf(2**(-prec-3),prec):
                break # exit while cycle
            
            # accumulate    
            for i in range(mmax+1):
                der_r[i] = der_r[i] + deltader[i]
            
            r02k = r02k*r0*r0
            k+=1
        else:
            raise ValueError("k=1000 is reached")
            
        der_z = PF.zeros(mmax+1,prec)
        der_z[0] = der_r[0]
        # see nb TAG 2
        der_z[1] = (-8 + 6*root2)*der_r[1]
        der_z[2] = (32 - 22*root2)*der_r[1] + (136 - 96*root2)*der_r[2] 
        der_z[3] = 6*(-32 + 23*root2)*der_r[1] + 24*(-65 + 46*root2)*der_r[2] \
                    + 16*(-140 + 99*root2)*der_r[3]
    else:
        raise ValueError("Nonzero a,b not yet supported")
    return der_z
                
def scalar_diag (eps0, opdim0, mmax, a = 0, b = 0, prec = 212 ):
    """ returns z-derivatives up to mmax along the diagonal at z=1/2 """
    eps = pf(eps0,prec)
    opdim = pf(opdim0,prec)    
    if opdim <= eps:
        raise ValueError("Scalar_diag call at or below unitarity bound")
    seed = scalar_seed (eps, opdim, a, b , prec)
    if (a != 0 or b != 0):
        raise ValueError("Nonzero a,b not yet supported")
        
    if mmax<= 2:
        return seed[:mmax+1]
    # mmax >2
    h = PF.zeros(mmax+1,prec)
    h[:3] = seed
    
    for m in range(3,mmax+1):
        # TAG 5
        c1 = 6*(2 + eps) - 2*m
        c2 = 4*(5 - 9*eps + 2*m*(-4 + eps + m) + 2*opdim*(-2*(1 + eps) + opdim))
        c3 = -8*(-((-3 + m)*(37 + 11*eps + 2*m*(-9 - 2*eps + m))) + 3*opdim*(-2*(1 + eps) + opdim))
        h[m] = c1*h[m-1]+c2*h[m-2]+c3*h[m-3]
        if m>= 4:
            c4 = -16*(-3 + m)*((-4 + m)*(3 - 9*eps + m*(-5 + 2*eps + m)) + (-9 + 2*m)*opdim*(-2*(1 + eps) + opdim))
            h[m]+=c4*h[m-4]   
        if m>=5:
            c5 = -32*(-5 + m)**2*(-4 + m)*(-3 + m)*(-5 - eps + m)
            h[m] += c5*h[m-5]
            
    return h           
    
def spin_diag (eps0, l, opdim0, mmax, a = 0, b = 0, prec = 53 ):
    """ returns z-derivatives up to mmax along the diagonal at z=1/2 """
    #print eps, l, opdim, mmax
    eps = pf(eps0,prec)
    opdim = pf(opdim0,prec)  
    if opdim < l+ 2*eps:
        raise ValueError("spin_diag call below unitarity bound: "+str(opdim))
    if (a != 0 or b != 0):
        raise ValueError("Nonzero a,b not yet supported")
    seed = spin_seed (eps, l, opdim, a, b , prec)        
    if mmax<= 3:
        return seed[:mmax+1]
    # mmax >2
    h = PF.zeros(mmax+1,prec)
    h[:4] = seed
    
    for m in range(4,mmax+1):
        # TAG 5-L
        c1 = -2*(-7 - 6*eps + m)
        c2 = 4*(22 - 29*eps - 9*eps**2 + 2*l**2 + m*(-17 + 3*m) + 4*eps*(l + m - opdim) + 2*(-2 + opdim)*opdim)
        c3 = -8*(318 + 3*(57 - 7*eps)*eps + 3*(-69 + (-33 + eps)*eps)*m + 2*(22 + 7*eps)*m**2 - 3*m**3 + (5 + 6*eps)*(2*eps*l + l**2 + opdim*(-2*(1 + eps) + opdim)))
        c4 = -16*((-4 + m)*(-126 + 262*eps + 58*eps**2 + (121 - 2*eps*(47 + 7*eps))*m + (-34 + 8*eps)*m**2 + 3*m**3) + 4*l*(2*eps + l)*(-1 + opdim)*(-1 - 2*eps + opdim) + (72 - eps + m*(-35 - 2*eps + 4*m))*(2*eps*l + l**2 + opdim*(-2*(1 + eps) + opdim)))
       
        h[m] = c1*h[m-1]+c2*h[m-2]+c3*h[m-3]+c4* h[m-4]

        if m>= 5:
            c5 = 32*(-4 + m)*(-((-5 + m)*(-570 - 226*eps - 2*eps**2*(-7 + m) + 2*eps*(48 - 5*m)*m + m*(299 + m*(-52 + 3*m)))) + 4*l*(2*eps + l)*(-1 + opdim)*(-1 - 2*eps + opdim) + (-29 + 6*eps*(-6 + m) + 5*m)*(2*eps*l + l**2 + opdim*(-2*(1 + eps) + opdim)))
            h[m]+=c5*h[m-5]   

        if m>=6:
            c6 = 64*(-5 + m)*(-4 + m)*(-6 - eps + m)*((-6 + m)*(11 - 31*eps + m*(-8 + 5*eps + m)) + (-13 + 2*m)*(2*eps*l + l**2 + opdim*(-2*(1 + eps) + opdim)))
            h[m] += c6*h[m-6]

        if m>=7:
            c7 = 128*(-7 + m)**2*(-6 + m)*(-5 + m)*(-4 + m)*(-7 - eps + m)*(-6 - eps + m)
            h[m] += c7*h[m-7] 
            
    return h           
        
def cb_ders (eps0, l, opdim0, nmax, mmax, a = 0, b = 0, prec = 212 ):  
    """ returns the derivatives h[n][m] with n<=nmax, m<= mmax+2*(nmax-n) """
    eps = pf(eps0,prec)
    opdim = pf(opdim0,prec)  

    h = [PF.zeros(mmax +2*(nmax-n)+1,prec) for n in range(nmax+1)]
    
    if l==0:
        h[0] = scalar_diag(eps, opdim, mmax+2*nmax, a, b, prec)
    else:
        h[0] = spin_diag(eps, l, opdim, mmax+2*nmax, a, b, prec)

    
    for n in range(1,nmax+1):
        for m in range(mmax+2*(nmax-n)+1):
            #h[n][m] = 0 already set to zero above
            # TAG 6
            if n>=2:
                h[n][m] = h[n][m] +  (4*(-1 + n)*(-2*(3 + eps) + 3*m + 4*n)*h[-2 + n][1 + m])/(-1 + 2*eps + 2*n) 
           
            if n>=2:
                h[n][m] = h[n][m] +  (2*(-1 + n)*h[-2 + n][2 + m])/(-1 + 2*eps + 2*n)

            if m>=1:
                h[n][m] = h[n][m] +  (4*m*(m*(-13 + 2*eps + 12*n) + 2*(11 + eps + n*(-17 - 2*eps + 6*n)) + m**2)*h[-1 + n][-1 + m])/(-1 + 2*eps + 2*n)

            h[n][m] = h[n][m] +  (2*((-5 + m)*m + 8*m*n + 4*eps*(l + m + n - opdim) - 2*(1 + 2*eps + n - (-2 + opdim)*opdim) + 2*l**2 + 4*n**2)*h[-1 + n][m])/(-1 + 2*eps + 2*n)

            h[n][m] = h[n][m] + (1 - (-7 + m + 6*n)/(-1 + 2*eps + 2*n))*h[-1 + n][1 + m]

            h[n][m] = h[n][m] + h[-1 + n][2 + m]/(2 - 4*eps - 4*n)

            if m>=3:
                h[n][m] = h[n][m] + 8*(-2 + m)*(-1 + m)*m*h[n][-3 + m]

            if m>=2:
                h[n][m] = h[n][m] + 4*(-1 + m)*m*h[n][-2 + m]

            if m>=1:
                h[n][m] = h[n][m] + -2*m*h[n][-1 + m]
    
    return h
    
def pochhammer (delta , n):
    if n == 0:
        return pf(1,delta.prec)
    else:
        ret = pf(1,delta.prec)
        for i in range(n):
            ret = ret * (delta+i)
        return ret  
        

def coeff_tensor (ds0, nmax, mmax, prec = 212 ):  
    """ returns the tensor t which converts derivatives h[n][m] with n<=nmax, m<= mmax+2*(nmax-n) 
    into coefficients to be fed into linear programming. Format of the tensor
    linprog[n][m] = t[n][m][n1][m1] h[n1][m1]
    Coefficients with even m are returned zero which is their theoretical value which is not computed, just set"""
    ds = pf(ds0,prec)
  
    
    a = [[ (-1)**(m+n)* pochhammer(ds-n+1,n)* pochhammer(2*ds-2*n-m+1,m)/pf(2,prec)**(2*ds-2*n-m) 
            for m in range(mmax +2*(nmax-n)+1) ] for n in range(nmax+1)] # these are derivatives of ((1/2-t)^2-s)^ds
            
        
    t = [[
        [[ 0 if m%2==0 or m1>m or n1>n else a[n-n1][m-m1]\
        /pf(str(fac(n-n1)*fac(m-m1)*fac(n1)*fac(m1)),prec) #*fac(n)*fac(m) if want to compute derivatives
        # current version computes Taylor coefficients
            for m1 in range(mmax +2*(nmax-n1)+1) ] for n1 in range(nmax+1)]
        for m in range(mmax +2*(nmax-n)+1) ] for n in range(nmax+1)]   
  
    return t

def contract (t, h, nmax, mmax):
    """ returns conformal block contracted with a coeff_tensor, flattened, 
    and rounded to double """
    return [ PF.to_double(sum([t[n][m][n1][m1]*h[n1][m1] 
    for n1 in range(n+1) for m1 in range(m+1) ])) 
            for n in range(nmax+1) 
                        for m in range(1,2*(nmax-n)+mmax+1,2)]
                     
                        
                        
                        
                        
                        
