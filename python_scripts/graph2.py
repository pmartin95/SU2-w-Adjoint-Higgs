import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

def input_data(filename, num_col):
    data = []
    for i in range(num_col):
        data.append([])
    f = open(filename, "r")
    contents = f.readlines()
    f.close()
    for chunk in contents:
        line = chunk.split()
        for i in range(len(data)):
            data[i].append(float(line[i]))
    for x in data:
        yield x

beta,W22,W22error = input_data("rect2x2.txt",3)
_,W11,W11error = input_data("rect1x1.txt",3)
_,W21,W21error = input_data("rect2x1.txt",3)
_,W33,W33error = input_data("rect3x3.txt",3)
_,W32,W32error = input_data("rect3x2.txt",3)
beta = np.array(beta)/4
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
ax.errorbar(beta,W11,W11error,fmt='.',ecolor='red',label=r"$W(1,1)$")
ax.errorbar(beta,W21,W21error,fmt='.',ecolor='red',label=r"$W(1,2)$")
ax.errorbar(beta,W22,W22error,fmt='.',ecolor='red',label=r"$W(2,2)$")
ax.errorbar(beta,W32,W32error,fmt='.',ecolor='red',label=r"$W(3,2)$")
ax.errorbar(beta,W33,W33error,fmt='.',ecolor='red',label=r"$W(3,3)$")
ax.legend()
ax.set_yscale('log')
ax.set_xlabel(r'{$\beta$}')
ax.set_ylabel(r'{$W(M,N)$}')
plt.savefig("rectplot.png")
plt.clf()
##########################


fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
#X(1,1)
ax.errorbar(beta,-np.log(W11),np.divide(W11error,W11),fmt='.',ecolor='red',label=r"$\chi(1,1)$")
#X(2,2)
y = np.divide(np.array(W22)*np.array(W11), np.array(W21)**2)
error = np.sqrt( (np.divide(W22error,W22))**2 + (np.divide(W11error,W11))**2 + 2*(np.divide(W21error,W21))**2   )
mask1 = (y > 0.0) & (y <= 1.0)
mask2 = np.divide(error,y) < 0.1
mask = mask1 & mask2
ax.errorbar(np.array(beta)[mask],-np.log(y[mask]),error[mask],fmt='.',ecolor='red',label=r"$\chi(2,2)$")
#X(3,3)
y = np.divide(np.array(W33)*np.array(W22), np.array(W32)**2)
error = np.sqrt( (np.divide(W33error,W33))**2 + (np.divide(W22error,W22))**2 + 2*(np.divide(W32error,W32))**2   )
mask1 = (y > 0.0) & (y <= 1.0)
mask2 = np.divide(error,y) < 0.1
mask = mask1 & mask2
print(y)
print(error)
ax.errorbar(np.array(beta)[mask],-np.log(y[mask]),error[mask],fmt='.',ecolor='red',label=r"$\chi(3,3)$")

ax.set_yscale('log')
ax.legend()
ax.set_xlabel(r'{$\beta$}')
ax.set_ylabel(r'{$\chi(I,J)$}')
ax.set_yticks([0.01,0.1,1.0,10.0])
ax.set_yticklabels(['0.01','0.1','1.0','10'])
ax.set_xticks([0,0.25,0.5,0.75])
ax.set_xlim([0,0.75])
plt.savefig("creutz_ratio.png")


plt.clf()
mask = (y > 0.0) & (y < 1.0)
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
ax.errorbar(np.array(beta)[mask],y[mask],error[mask],fmt='.')
plt.savefig("test.png")