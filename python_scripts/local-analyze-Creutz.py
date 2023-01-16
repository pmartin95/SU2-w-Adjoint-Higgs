import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path
from os.path import exists
p = Path('.')
plt.rcParams['text.usetex'] = True

def input_data(filename):
    data = []
    with open(filename, "r") as f:
        contents = f.readlines()
        data = [float(x) for x in contents[0].split()]
    return data



rect_pattern = r"rect(\d)x(\d)beta(\d+\.?\d*).txt"
files = list(p.glob('**/*.txt'))

x = []
y = []
err = []

# 1x1 Creutz ratio
for file_name in files:
    m = re.search(rect_pattern,str(file_name))
    if m:
        # print(m.group(0))
        if(m.group(1) == '1' and m.group(2) == '1' and input_data(file_name)[1] > 0.0):

            x.append(input_data(file_name)[0]/4.0)
            y.append(-np.log(input_data(file_name)[1]))
            err.append(input_data(file_name)[2]/input_data(file_name)[1])

# 2x2 Creutz ratio
x2 = []
y2 = []
err2 = []
for file_name in files:
    m = re.search(rect_pattern,str(file_name))
    if m:
        if(m.group(1) == '2' and m.group(2) == '2' ):
            filename_template = "rect{}x{}beta{}.txt"
            W11 = input_data(filename_template.format(1,1,m.group(3)))[1]
            W22 = input_data(filename_template.format(2,2,m.group(3)))[1]
            W21 = input_data(filename_template.format(2,1,m.group(3)))[1]
            W11err = input_data(filename_template.format(1,1,m.group(3)))[2]
            W22err = input_data(filename_template.format(2,2,m.group(3)))[2]
            W21err = input_data(filename_template.format(2,1,m.group(3)))[2]
            temp = W11 * W22 / W21 **2
            
            if temp > 0.0:
                X22err = np.sqrt((W11err/W11)**2 + (W22err/W22)**2 + 2*(W21err/W21)**2)
                x2.append(input_data(file_name)[0]/4.0)
                y2.append(-np.log(temp))
                err2.append(X22err)

# 3x3 Creutz ratio
x3 = []
y3 = []
err3 = []
for file_name in files:
    m = re.search(rect_pattern,str(file_name))
    if m:
        if(m.group(1) == '3' and m.group(2) == '3' ):
            filename_template = "rect{}x{}beta{}.txt"
            W22 = input_data(filename_template.format(2,2,m.group(3)))[1]
            W33 = input_data(filename_template.format(3,3,m.group(3)))[1]
            W32 = input_data(filename_template.format(3,2,m.group(3)))[1]
            W22err = input_data(filename_template.format(2,2,m.group(3)))[2]
            W33err = input_data(filename_template.format(3,3,m.group(3)))[2]
            W32err = input_data(filename_template.format(3,2,m.group(3)))[2]
            temp = W33 * W22 / W32 **2
            
            if temp > 0.0:
                X33err = np.sqrt((W33err/W33)**2 + (W22err/W22)**2 + 2*(W32err/W32)**2)
                x3.append(input_data(file_name)[0]/4.0)
                y3.append(-np.log(temp))
                err3.append(X33err)


fig, ax = plt.subplots(figsize=(4, 4), tight_layout=True)
plt.yscale('log')
x1 = [ 2.3/4,2.3/4]
y1 = [0.0,10.0]
ax.errorbar(x,y,err,fmt='.')
ax.errorbar(x2,y2,err2,fmt='.')
ax.errorbar(x3,y3,err3,fmt='.')
ax.plot(x1,y1)
ax.set_xlabel(r'{$1/g^2$}')
ax.set_ylabel(r'{$\chi(I,J)$}')
ax.set_yticks([0.01,0.1,1.0,10.0])
ax.set_yticklabels(['0.01','0.1','1.0','10'])
ax.set_xticks([0,0.25,0.5,0.75])
# ax.set_xlim([0,0.75])
# p = Path('./../testdata5/')
# x = []
# y = []
# err = []
# files = list(p.glob('**/*.txt'))
# for file_name in files:
#     m = re.search(rect_pattern,str(file_name))
#     if m:
#         if(m.group(1) == '1' and m.group(2) == '1'):
#             x.append(input_data(file_name)[0]/4.0)
#             y.append(-np.log(input_data(file_name)[1]))
#             err.append(input_data(file_name)[2]/input_data(file_name)[1])
# print(x,y,err)
# ax.errorbar(x,y,err,fmt='.')
plt.savefig("testCreutzRatio.png")