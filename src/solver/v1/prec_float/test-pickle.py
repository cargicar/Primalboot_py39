import prec_float as PF
from prec_float import prec_float as pf
import pickle

a=pf("1.2345",212)

file = open ("data.txt","w");

a.to_file(file)

file.close()

