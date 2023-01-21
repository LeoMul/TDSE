from matplotlib import pyplot as plt
import imageio
import numpy as np 

#reading in basic data structure.
#first line is x_data
#each subsequent line is the wave function at that point in time. 
#maybe think about including delta_t somehow to get absolute time.
#read in as many files as you like, i plot abs, Re and Im.
data = np.loadtxt("psi_matrix_square_mod.dat")
real_part = np.loadtxt("psi_matrix_real.dat")
im_part = np.loadtxt("psi_matrix_imag.dat")
plt.rcParams['text.usetex'] = True
x_array = data[0,:]
size = data.shape

#the (discrete) time steps to be plotted: first, final, step
time = range(1,size[0]-1,5)

#for plotting a square well if you're so inclined.
square = [0.0,0.9]
square_base_1 = [5.0,5.0]
square_base_2 = [10.0,10.0]

def create_frame(time_index):
    fig = plt.figure()
    max_psi = np.max(data[1,:])
    plt.plot(x_array,np.sqrt(data[time_index,:]),label = "$|\psi|$")
    plt.plot(x_array,real_part[time_index,:],label = "Re $\psi$")
    plt.plot(x_array,im_part[time_index,:],label = "Im $\psi$")
    plt.title(str(time_index))
    plt.ylim([-max_psi,max_psi])
    plt.xlim([x_array[0],x_array[len(x_array)-1]])
    #x = np.linspace(-50,50,1000)
    #v = 0.5*x*x
    #plt.plot(x,v,'k',label = "Harmonic Potential")
    plt.xlabel("$x$")
    plt.legend()
    #plt.ylabel("$|\psi|^2$", fontsize=16)
    
    #maybe you could also use a pdf here, i haven't tested
    plt.savefig(f'./img/img_{t}.png', 
                transparent = False,  
                facecolor = 'white'
               )
    plt.close()

#makes the frames
for t in time:
    create_frame(t)

#makes the pngs into a gif. 
frames = []
for t in time:
    image = imageio.v2.imread(f'./img/img_{t}.png')
    frames.append(image)
#print(size)

imageio.mimsave('./test.gif', # output gif
                frames,          # array of input frames
                fps = 15,loop = 0)         # optional: frames per second
