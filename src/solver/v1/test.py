import tables.tables as cb
reload(cb)
import prec_float.prec_float as PF

dl = cb.make_dlist (eps = 1., lmax=10,dim_all_max=10, step=0.1, dim_scalar_min =2)

tab = cb.CB_Table(eps = 1., dlist = dl, nmax = 2, mmax = 1)

tab.save("a.txt")

tab1 = cb.CB_Table(FILE="a.txt")

#print tab1.eps, tab1.nmax, tab1.mmax, tab1.prec

#print max([abs(tab1.dlist[l][i] -tab.dlist[l][i]) for l in range(len(tab1.dlist))
#        for i in range(len(tab1.dlist[l])) ])

#print max([abs(tab1.dlist[l][i] -tab.dlist[l][i]) for l in range(len(tab.dlist))
#        for i in range(len(tab.dlist[l])) ])

#print max([abs(PF.to_double(tab1.table[l][i][n][m] -tab.table[l][i][n][m])) for l in range(len(tab1.table))
#        for i in range(len(tab1.table[l])) for n in range(tab.nmax+1) for m in range(2*(tab.nmax-n)+tab.mmax+1)])
#        
#print max([abs(PF.to_double(tab1.table[l][i][n][m] -tab.table[l][i][n][m])) for l in range(len(tab.table))
#        for i in range(len(tab.table[l])) for n in range(tab.nmax+1) for m in range(2*(tab.nmax-n)+tab.mmax+1)])

s = cb.Sigma_Table(ds=0.125, cbtab = tab)
s.save("b.txt")

s1= cb.Sigma_Table(FILE = "b.txt")
