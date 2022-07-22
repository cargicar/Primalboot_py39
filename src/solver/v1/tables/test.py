import cb_table as cb

dl = cb.make_dlist (eps = 1., lmax=10,dim_all_max=10, step=0.1, dim_scalar_min =2)

tab = cb.CB_Table(eps = 1., dlist = dlist, nmax = 2, mmax = 1)

tab.compute()

tab.save("a.txt")

tab1 = cb.CB_Table()

tab1.read("a.txt")






