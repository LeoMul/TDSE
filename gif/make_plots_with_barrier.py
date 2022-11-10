from matplotlib import pyplot as plt
import imageio
import numpy as np 

data = np.loadtxt("fourierdatagif.dat")
plt.rcParams['text.usetex'] = True
x_array = data[0,:]
size = data.shape

time = range(1,size[0]-1,1)

square = [0.0,0.9]
square_base_1 = [5.0,5.0]
square_base_2 = [10.0,10.0]

def create_frame(time_index):
    fig = plt.figure()
    plt.plot(x_array,data[time_index,:])
    plt.plot(square_base_1,square,'r')
    plt.plot(square_base_2,square,'r')
    plt.plot([5,10],[0.9,0.9],'r')

    plt.title(str(time_index))
    plt.ylim([0,1])
    plt.xlim([-3,25])

    plt.xlabel("x")
    plt.ylabel("$|\psi|^2$", fontsize=16)
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
                fps = 15,loop = 0)         # optional: frames per second
