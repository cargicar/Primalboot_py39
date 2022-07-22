import utils.logging as logging

# true for MPFR, false for numpy
precision = False

# the branch&bound min as a multiplier of the current cost
threshold_multiplier=1.0e-4

# logging
logfile='../logs/log'
debug_level=logging.DEBUG

# set to None to not use a db
#dbfile=None
dbfile='../db/runs.db'

# parallelization of a single point 
point_parallel=False
point_poolsize=8

