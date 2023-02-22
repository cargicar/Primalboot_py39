import matplotlib.pyplot as plt
import numpy as np

with open("agave_multiplerun_out_n20.txt") as file:
    lines20=file.readlines()
with open("agave_multiplerun_out_n16.txt") as file:
    lines16=file.readlines()
sigmas16=[]
eps16=[]
sigmas20=[]
eps20=[]
for line in lines20:
    for i in range(len(line)):
        if line[i]=="=":
            n=0
            while line[i+n]!=":":
                n+=1
            sigmas20.append(float(line[i+1:i+n]))
            eps20.append(float(line[i+n+4:len(line)-1]))

for line in lines16:
    for i in range(len(line)):
        if line[i]=="=":
            n=0
            while line[i+n]!=":":
                n+=1
            sigmas16.append(float(line[i+1:i+n]))
            eps16.append(float(line[i+n+4:len(line)-1]))
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(sigmas16, eps16, s=5, c='b', marker="s", label='n=16')
ax1.scatter(sigmas20, eps20, s=5, c='r', marker="o", label='n=20')
plt.legend(loc='upper left')

plt.title('$\Delta_{\epsilon}$ vs $\Delta{\sigma}$')
plt.xlabel('$\Delta_{\sigma}$')
plt.ylabel('$\Delta_{\epsilon}$')
plt.show()
#plt.savefig("sigma_vs_epsilon_n20.png")


