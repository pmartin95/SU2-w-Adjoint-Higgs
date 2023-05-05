import numpy as np

def f(x):
    val = np.arange(x-4.3-.2,x-4.3+.2,0.1)
    val = val[val < 0.0]
    return val

betas = np.arange(2.0, 3.1, 0.1)


with open("beta_values.txt", "w")as beta_file,open("m2_values.txt","w") as m2_file:
    for beta in betas:
        m2s = f(beta)
        for m2 in m2s:
            beta_file.write(f"{round(beta,3)}\n")
            print(f"{round(beta,3)} {round(m2,3)}")
            m2_file.write(f"{round(m2,3)}\n")