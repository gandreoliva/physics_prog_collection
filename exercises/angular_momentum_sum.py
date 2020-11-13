import numpy as np

j1 = 1
j2 = 1/2
m1s = np.arange(-j1,j1+0.1,1)
m2s = np.arange(-j2,j2+0.1,1)


print("j1 =",j1, " j2 =",j2)
print("----------------------------")

j = j1 + j2
i = 0
while j >= abs(j1-j2):
    i = i + 1
    print("j =", j)
    ms = np.arange(-j,j+0.1,1)
    for m in ms:
        print("   m =", m)
        for m1 in m1s:
            for m2 in m2s:
                if m == m1 + m2:
                    print("      m1 =",m1,"  m2 =",m2)
    j = j1 + j2 - i
    print("----------------------------")
