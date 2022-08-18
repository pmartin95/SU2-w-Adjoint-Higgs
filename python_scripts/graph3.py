import numpy as np
import matplotlib.pyplot as plt
# plt.rcParams['text.usetex'] = True

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

f = open("count.txt","r")
contents = f.readlines()
words = []
for val in contents:
    words.append(float(val)) 

plt.hist(words)
plt.savefig('test.png')
