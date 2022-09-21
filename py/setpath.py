import sys
#from version import version

#if version==1:
#    path='../bootstrap/v1'
#elif version==2:
#    path='../bootstrap/v2'
#path='../bootstrap/merge'
path='../src'

if not(path in sys.path):
    sys.path = [path]+sys.path

