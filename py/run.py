import time
import sys
import os
import psutil


# Choose if hotstart or not-hotstart
hot=True

# Start timing
start=time.perf_counter()
 
# Create dir to save initial feasible solution
for f in os.listdir("jobdir/"):
    os.remove(os.path.join("jobdir/", f))                                                                                                   
jobdir="jobdir/"
print("starting")
if hot:
    import reaper_funcs_hotstart as RFhs
    RFhs.findbound(0.51815,jobdir)
else:
    import reaper_funcs as RF
    RF.findbound(0.51815)
    
#Record process to read memory usage
process = psutil.Process(os.getpid())
print(f" Process memory {process.memory_info().rss}")  # in bytes 

# end timing
end=time.perf_counter()
#milliseconds= 1000*(end-start)                                                                                      
seconds=(end-start)
minutes=seconds/60
print(f"Code took {minutes:.3f}m")

# Force exit in case the code did not exit
sys.exit()
