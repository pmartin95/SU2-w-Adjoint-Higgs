import numpy as np
import matplotlib.pyplot as plt
# plt.rcParams['text.usetex'] = True
f = open("plaqdata.txt", "r") 
contents = f.readlines()
f.close()
beta = []
Y = []
E = []
for x in contents:
    temp = x.split()
    beta.append(float(temp[0]))
    Y.append(float(temp[1]))
    E.append(float(temp[2]))

# fig, ax = plt.subplots()
plt.yscale("log")
# ax.set_yscale('log')
plt.errorbar(beta,Y,yerr=E,fmt='.',markersize=1,ecolor='red')
# ax.set_xlabel(r'$\Beta$')
plt.xlabel("{\Beta}")
plt.savefig("stringtension.png")