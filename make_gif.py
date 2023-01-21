from matplotlib import pyplot as plt
import imageio
import numpy as np 
#makes a gif.
data = np.loadtxt("fourierdatagif_abs_barrier.dat")
plt.rcParams['text.usetex'] = True
x_array = data[0,:]
size = data.shape
print(size)
delta_t = 0.025#data has only every 10th printed.
time = range(1,size[0],1)

square = [0.0,0.9]
square_base_1 = [5.0,5.0]
square_base_2 = [10.0,10.0]

def create_frame(time_index):
    fig = plt.figure()
    plt.plot(x_array,data[time_index,:],label = "explicit")

    plt.title("t = " + str(time_index*delta_t))
    plt.ylim([-1,1])
    plt.xlim([-20,20])
    #x = np.linspace(-50,50,1000)
    #v = 0.5*x*x
    #plt.plot(x,v,'k',label = "Harmonic Potential")
    plt.xlabel("$x$")
    plt.ylabel("$u(x,t)$")
    plt.legend()
    #plt.ylabel("$|\psi|^2$", fontsize=16)
    plt.savefig(f'./img/img_{t}.png', 
                transparent = False,  
                facecolor = 'white'
               )
    plt.close()

for t in time:
    print(t)
    create_frame(t)

frames = []
for t in time:
    image = imageio.v2.imread(f'./img/img_{t}.png')
    frames.append(image)
#print(size)

imageio.mimsave('./example.gif', # output gif
                frames,          # array of input frames
                fps = 30,loop = 0)         # optional: frames per second
