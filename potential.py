import numpy as np 

import matplotlib.pyplot as plt 

N = 1000

x = np.linspace(-20,25,N)




def potential(x,n):
    pos = 0.0
    V = np.zeros(N)


    for i in range(0,n): 
        pos = pos + (1/n)*x[-1]

        V = V + 1.0*np.exp(-10.0*(x-pos)**2) 
    return V

plt.plot(x,potential(x,n=20))
plt.show()