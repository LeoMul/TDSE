from matplotlib import pyplot as plt
import imageio
import numpy as np 

data = np.loadtxt("fourierdatagif_abs_barrier.dat")
real_part = np.loadtxt("fourierdatagif_real_barrier.dat")
im_part = np.loadtxt("fourierdatagif_imag_barrier.dat")
plt.rcParams['text.usetex'] = True
x_array = data[0,:]
size = data.shape

time = range(1,size[0]-1,1)

square = [0.0,0.9]
square_base_1 = [4.0,4.0]
square_base_2 = [4.5,4.5]
square_top_x = [4.0,4.5]
square_top_y = [0.9,0.9]
def create_frame(time_index):
    fig = plt.figure()
    plt.plot(x_array,(data[time_index,:]),label = "$|\psi|$")
    plt.plot(square_base_1,square,'r')
    plt.plot(square_base_2,square,'r')
    plt.plot(square_top_x,square_top_y,'r')
    plt.plot(x_array,real_part[time_index,:],label = "Re $\psi$")
    plt.plot(x_array,im_part[time_index,:],label = "Im $\psi$")
    plt.title(str(time_index))
    plt.ylim([-1.1,1.1])
    plt.xlim([-15,20])
    #x = np.linspace(-50,50,1000)
    #v = 0.5*x*x
    #plt.plot(x,v,'k',label = "Harmonic Potential")
    plt.xlabel("x")
    plt.legend()
    #plt.ylabel("$|\psi|^2$", fontsize=16)
    plt.savefig(f'./img/img_{t}.png', 
                transparent = False,  
                facecolor = 'white'
               )
    plt.close()

for t in time:
    #print(t)
    create_frame(t)

frames = []
for t in time:
    image = imageio.v2.imread(f'./img/img_{t}.png')
    frames.append(image)
#print(size)

imageio.mimsave('./example.gif', # output gif
                frames,          # array of input frames
                fps = 40,loop = 0)         # optional: frames per second
