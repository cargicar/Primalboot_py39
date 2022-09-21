# Primalboot_py39
Primal simplex implementation for the conformal bootstrap. This is an updated version from the original one created by S.Rychkov, S. El-Showk. 

In this version we have portated the whole code to Python 3.9 and activated multiprocessing features (paralellization). 

In order to understand this code, is convenient to look first at the original paper [Solving the 3d Ising Model with the Conformal Bootstrap II. c-Minimization and Precise Critical Exponents](https://arxiv.org/pdf/1403.4545.pdf)

## Simplest use

0. It is recommended to create a python==3.9 environment to install and isolate the code from the global environment.  The way of doing it dependent on the environment manager each user uses. My prefered environment and package manager is [conda](https://docs.conda.io/projects/conda/en/latest/index.html) . Follow the instructions there to install conda and to create local environments on your local machine. 

1. In a local terminal, activate your local environment and run *pip install -r requirements.txt* (Need to update requirements!)

2. Go to the directory */src* and run the command *./compile.sh*. This will create the neccesary Cython modules. (Need to include mpfr and mpfi into the requirements file!)

3. In a terminal go to  *py* directory and open *config* file. There is a number of parameters that can be set in config file, perhaps the most importants are: 
Number of *cores*  available in your local enviroment (or the amount you want to use) \
The highest order for the Taylor expansion of the conformal blocks *N*=(*n* + 1) (*n* + 2)/2 .Due to storage constraints on GitHub, this repository only contains the tables corresponding to *n*=6, 8, 10, 12.\
Boolean variable to enable or disable parallelization *point_parallel* \
Boolean variable to enable or disable printing of status during excecution *show_status_after*

4. The file *probspec.py* controls theory parameters, such as space time dimension *dim* and spin of the operator OPE coefficient to be optimized *opespin*, but modifying this file will requiere building up new tables. If you want to run with the preloaded tables here, don't modify  *probspec.py*.

5. In the file *run.py* it can be choose if the simplex algorithm wants to be run with hotstart or without hotstart (see paper for details), by setting the variable *hot* to *True* or *False* respectively.

6. Run *python run.py > 'Name_output.txt'* changing 'Name_output.txt' for whatever name you want to give to the output file.


