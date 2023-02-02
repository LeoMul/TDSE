import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm
import imageio
from matplotlib.pyplot import figure 

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 32})


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
    #fig, ax = plt.subplots()
    #print(matrix)
    fig = plt.figure()
    ax = plt.axes(projection = '3d')
    x, y = np.meshgrid(x_array, y_array)

    #c = ax.pcolormesh(x, y, matrix+0.2*potential, cmap=cm.coolwarm,vmin= 0.0,vmax = 1.0)
    #fig.colorbar(c, ax=ax)
    #ax.scatter3D(x,y,matrix)
    #ax.scatter3D(x,y,potential)
    surf2 = ax.plot_surface(x, y, potential,alpha = 0.1, cmap='Blues',linewidth=0, antialiased=False)
    surf1 = ax.plot_surface(x, y, 50.0*matrix + 0.1, cmap='Blues',linewidth=0, antialiased=False)
    ax.set_zlim(0, 50.0)

    

    fig.set_size_inches(16,12 )
    time = round(index * every * delta_t,10)
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    stepsize = 5.0
    #ticks_array = np.arange(-x[], 10.0+stepsize, stepsize)
    #ax.xaxis.set_ticks(ticks_array)
    #ax.yaxis.set_ticks(ticks_array)
    #plt.xlim([-10,10])
    #c.set_label("|psi|")
    plt.title(str(time))
    plt.savefig(f'./img/img_{index}.png', 
                transparent = False,  
                facecolor = 'white', dpi=200
               )
    plt.close()

def create_another_using_mayami(x_array,y_array,matrix,index,delta_t,every):
    from mayavi import mlab
    x, y = np.meshgrid(x_array, y_array)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf3 = mlab.surf(x, y, matrix, colormap='Blues')
    surf4 = mlab.surf(x, y, potential, colormap='Oranges')
    surf3.actor.property.opacity = 0.5
    surf4.actor.property.opacity = 0.5
    fig.scene.renderer.use_depth_peeling = 1
    mlab.savefig(f'./img/img_{index}.png')


potential = np.loadtxt("potential_gaussian_sea_strength10.000spacing2.000width0.250.dat")

raw_data = "abswavefunction_momentum10.000packetwidth2.000gaussian_sea_strength10.000spacing2.000width0.250.dat"

x_array,y_array,every,delta_t,matrices,time_steps = read_in_2d_schrodinger_data(raw_data)
#max_scale = max(matrices)

time = range(0,int(time_steps/2),1)
for t in time:
    print(t)
    create_heatmap_frame(x_array,y_array,matrices[t,:,:],t,delta_t,every)
    #create_another_using_mayami(x_array,y_array,matrices[t,:,:],t,delta_t,every)


frames = []
for t in time:
    image = imageio.v2.imread(f'./img/img_{t}.png')
    frames.append(image)


imageio.mimsave('./2d_test_pit.gif', # output gif
                frames,          # array of input frames
                fps = 10,loop = 0)         # optional: frames per second

