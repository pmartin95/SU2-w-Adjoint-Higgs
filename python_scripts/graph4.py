from cProfile import label
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

cold = list(input_data('coldW33.txt',1))[0]
print(cold)
hot = list(input_data('hotW33.txt',1))[0]
# print(type(cold))
plt.plot(range(len(cold)),cold,label="cold")
plt.plot(range(len(hot)),hot,label='hot')
plt.legend()
plt.xlim([0,1000])
plt.savefig('coldhot.png')