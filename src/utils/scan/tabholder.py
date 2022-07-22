import utils.logging as log

import solver.v1.tables.tables as cb
import solver.v1.tables.lp_table as lp_table

# this class is just used to hold the table file reference to be used by the
# parallel python processes

tab={}
lptab={}

def get_tabfile(filename):
    if filename in tab.keys():
        return tab[filename]
    else:
        log.info("Loading table file: "+ filename)
        tab[filename] = cb.CB_Table(FILE = filename)
    return tab[filename]

def get_lptab(sigma, tabfilename):
    if (sigma,tabfilename) in lptab.keys():
        return lptab[(sigma, tabfilename)]
    else:
        cbtab = get_tabfile(tabfilename)
        sigmatab = cb.Sigma_Table(ds=sigma, cbtab = cbtab)
        lpt = lp_table.LP_Table(sigmatab)
        lptab[(sigma, tabfilename)]=lpt
    return lptab[(sigma, tabfilename)]

