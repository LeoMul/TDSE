import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm
import imageio




zeros = np.zeros([4,6])
print(zeros)
zeros[0:2,0:4] = 1
print(zeros)


def read_in_2d_schrodinger_data(path):

    misc_data = np.loadtxt(path,max_rows=3)
    raw_data = np.loadtxt(path,skiprows=7)
    information = misc_data[0,:]
  
    size_of_x = int(information[0])
    size_of_y = int(information[1])
    every = int(information[2])
    delta_t = information[3]
    time_steps = int(information[4])

    x_array = np.zeros(size_of_x)
    y_array = np.zeros(size_of_y)

    for j in range(0,size_of_x):
        x_array[j] = misc_data[1,j]
    for j in range(0,size_of_y):
        y_array[j] = misc_data[2,j]
    
    matrices = np.empty([int(time_steps),int(size_of_y),int(size_of_x)])
    print(matrices.shape)
    for j in range(0,int(time_steps)):
        matrices[j,:,:] = raw_data[(j*int(size_of_y)):((j+1)*int(size_of_y)),0:int(size_of_x)]
    return x_array,y_array,int(every),delta_t,matrices,int(time_steps)


def create_heatmap_frame(x_array,y_array,matrix,index,delta_t,every):
    fig, ax = plt.subplots()
    #print(matrix)
    x, y = np.meshgrid(x_array, y_array)
    c = ax.pcolormesh(x, y, matrix, cmap=cm.coolwarm,vmin= 0.0,vmax = 1.0)
    fig.colorbar(c, ax=ax)
    time = round(index * every * delta_t,10)
    plt.xlabel("x")
    plt.ylabel("y")
    c.set_label("|psi|")
    plt.title(str(time))
    plt.savefig(f'./img/img_{index}.png', 
                transparent = False,  
                facecolor = 'white'
               )
    plt.close()

raw_data = "fourierdata_2d_gif_abs.dat"
x_array,y_array,every,delta_t,matrices,time_steps = read_in_2d_schrodinger_data(raw_data)

time = range(0,time_steps,1)
for t in time:
    print(t)
    create_heatmap_frame(x_array,y_array,matrices[t,:,:],t,delta_t,every)

frames = []
for t in time:
    image = imageio.v2.imread(f'./img/img_{t}.png')
    frames.append(image)


imageio.mimsave('./2d.gif', # output gif
                frames,          # array of input frames
                fps = 10,loop = 0)         # optional: frames per second

