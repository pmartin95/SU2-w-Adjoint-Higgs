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

hotAbsPoly = {}
for filename in glob.glob('hotAbsPoly*.txt'):
    hotAbsPoly[filename] =  list(input_data(os.path.join(os.getcwd(),filename),1))[0]

coldAbsPoly = {}
for filename in glob.glob('coldAbsPoly*.txt'):
    coldAbsPoly[filename] =  list(input_data(os.path.join(os.getcwd(),filename),1))[0]



regex1 = r"hotAbsPoly(........).txt"
regex2 = r"coldAbsPoly(........).txt"

for x in hotAbsPoly:
    p =re.compile(regex1)
    m = p.match(x)
    beta = m.groups(0)[0]
    coldName = "coldAbsPoly" + beta +".txt"
    plt.plot(hotAbsPoly[x])
    plt.plot(coldAbsPoly[coldName])
    plt.title(r"$\beta = $"+ beta)
    plt.ylabel(r"$\langle P \rangle$")
    plt.savefig("thermalize_poly"+beta+".png")
    plt.clf()

# data = list(input_data("randsu2.dat",1))[0]
# plt.hist(data)
# # plt.title(r"$\beta = ")
# plt.savefig("distr_rand2.png")