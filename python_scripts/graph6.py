import os, glob
import numpy as np
import matplotlib.pyplot as plt
import re

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

polyakovLoops = {}

for filename in glob.glob('polyakov*.txt'):
    polyakovLoops[filename] =  list(input_data(os.path.join(os.getcwd(),filename),1))[0]

regex = r"polyakov(........).txt"

for x in polyakovLoops:
    p =re.compile(regex)
    m = p.match(x)
    beta = m.groups(0)[0]
    plt.hist(polyakovLoops[x],bins=40)
    print(len(polyakovLoops[x]))
    plt.title(r"$\beta = $"+ beta)
    plt.savefig("distr"+beta+".png")
    plt.clf()

# data = list(input_data("randsu2.dat",1))[0]
# plt.hist(data)
# # plt.title(r"$\beta = ")
# plt.savefig("distr_rand2.png")