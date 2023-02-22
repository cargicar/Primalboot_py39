import numpy as np
s0=0.5179
s1=0.5187
sdelta=0.0001
siglist=np.arange(s0, s1, sdelta)
print(siglist)
for i in range(len(siglist)):
   with open('run' +str(i)+ '.py','w') as f:
       f.write(
          ('import time\n'
                'import sys\n'
                'import os\n'
                'import psutil\n'
                'hot=True\n'
                'start=time.perf_counter()\n'
                '# for f in os.listdir("jobdir/"):\n'
                '# os.remove(os.path.join("jobdir/", f)) \n'

                'jobdir="jobdir/"\n'
                'print("starting")\n'
                'if hot:\n'
                '    import reaper_funcs_hotstart as RFhs\n'
                '    RFhs.findbound({sig},jobdir)\n'
                'else:\n'
                '    import reaper_funcs as RF\n'
                '    RF.findbound({sig})\n'
                'process = psutil.Process(os.getpid())\n'
                'print(" Process memory %s " % process.memory_info().rss)\n'
                'end=time.perf_counter()\n'
                'seconds=(end-start)\n'
                'minutes=seconds/60\n'
                'print("Code took %s minutes" % minutes)\n'
                'sys.exit()\n'
                ).format(sig=siglist[i]))
