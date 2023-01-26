import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm
import imageio
from matplotlib.pyplot import figure 

Nxpoints = 512
Nypoints = 512
x = np.linspace(-5,5,Nxpoints)
y = np.linspace(-5,5,Nypoints)

def matrix_of_pos(spacing,xfirst,xlast,ylength):

    ybot = -ylength / 2.0
    #yspacing = ylength / (Ny-1) 
    #xspacing = (xlast - xfirst) / (Nx-1)

    Nx = int((xlast-xfirst)/spacing) + 1
    Ny = int((ylength)/spacing) + 1

    posmatrixx = np.zeros([Nx])
    posmatrixy = np.zeros([Ny])

    for i in range(0,Ny):
        posmatrixy[i]= ybot + i*spacing
    for i in range(0,Nx):
        posmatrixx[i] = xfirst + i *spacing
    print(posmatrixy)
    print(posmatrixx)

    return posmatrixx,posmatrixy




def potential_term(x,y,xcentre,ycentre,width,strength):
    V = np.zeros([len(y),len(x)])

    Vxvector = np.exp(- (x-xcentre)**2 / width)
    Vyvector = np.array([np.exp(- (y-ycentre)**2 / width)])
    
    V = strength * np.matmul(np.matrix(Vyvector).T,np.matrix(Vxvector))
    return V

def create_potential_sea(x,y,width,strength,spacing,xfirst):
    V = np.zeros([len(y),len(x)])
    xpositions,ypositions = matrix_of_pos(spacing,xfirst,x[-1],y[-1]-y[0])

    for i in range(0,len(xpositions)):
        for j in range(0,len(ypositions)):
            V = V + potential_term(x,y,xpositions[i],ypositions[j],width,strength)

    return V


V = create_potential_sea(x,y,0.5,1.0,2,0.0)
#print(V)
fig, ax = plt.subplots()

x, y = np.meshgrid(x, y)
c = ax.pcolormesh(x, y, V, cmap=cm.coolwarm,vmin= 0.0,vmax = 1.0)
fig.colorbar(c, ax=ax) 
plt.show()

