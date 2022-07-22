import utils.logging as logging
import os

################# Carlos edit  #################
cores =8
env_cores=os.cpu_count()
if env_cores != '':
    cores = int(env_cores)
    
################# Carlos edit  #################

maxruntime=240000
# true for MPFR, false for numpy
precision = True

# Number of iterations to save data from iter
iters=2000



# Tables choice
nmax =6 # Carlos edit
mmax=1

# threshold to stop minimization of cost function
threshold=1e-12
#threshold=1e-3

# the branch&bound min as a multiplier of the current cost
threshold_multiplier=1.0e-6

# set to None to not use a db
#dbfile=None
dbfile='../db/runs_test.db'

# max number of steps before cutting off an LP
max_steps=10000

# parallelization of a single point 
point_parallel=True

################# Carlos edit  #################
if point_parallel==True:
    print(f"number of cores {cores}")
    
################# Carlos edit  #################
point_poolsize=cores

# logging/debugging
#log_path=
# logfile='../logs/log_ope_lxplus' # Carlos edit
#debug_level=logging.DEBUG2  # Carlos edit

# show status after 500 iterations
show_status_after=True


tabledir='../tables/'
fudge=1e-14
