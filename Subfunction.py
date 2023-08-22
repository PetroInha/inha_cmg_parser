import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def read_out_grid(Ptime, filename = 'Forward_simulation.out', nx=64, ny=64, nz=20):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    ################### OSAT ##########################   

    grid_property = - np.ones((len(Ptime)+1,nz,ny,nx))
    t_flag = 0

    for i in range(len(lines)):
        # Read Line:
        line = lines[i].split()
        ## Read grid property:
        if [''.join(line)][0] == 'OilSaturation(fraction)':
            start = i + 2
            for j in range(start,len(lines)):
                if '*' * 70 in lines[j].split(' '):
                    end = j
                    break
                elif '*' * 131 + '\n' in lines[j].split(' '):
                    end = j
                    break
            val_lines = lines[start:end]

            if not any(['Plane K = ' in var for var in val_lines]):
                value = [''.join(val_lines[np.where(['All values are' in val for val in val_lines])[0][0]].split(' '))][0]
                value = value.split('Allvaluesare')[1].split('\n')[0]
                grid_property[t_flag] = value 
                t_flag += 1
                
            else:
                value_line_lst = list(np.where(['Plane K =' in  val for val in val_lines])[0])
                value_line_lst.append(len(val_lines))
                for v in range(len(value_line_lst)-1):
                    s, e = value_line_lst[v], value_line_lst[v+1]
                    values =  val_lines[s:e]
                    k_flag = int(val_lines[s:e][0].split(' Plane K = ')[-1].split(' ')[0]) - 1 
                    if any(['All values are' in val for val in values]):
                        value = [''.join(values[np.where(['All values are' in val for val in values])[0][0]].split(' '))][0]
                        value = value.split('Allvaluesare')[1].split('\n')[0]
                        grid_property[t_flag, k_flag] = float(value)

                    else:
                        i_index = np.where(['    I = ' in val for val in values])[0]
                        for tick in i_index: 
                            i_flags = np.array(values[tick].split('  I = ')[-1].split('     '),dtype = int) -1
                            f = lambda x: x if not x in ['', '\n'] else None
                            j_flags = [list(filter(f,val.split(' J= ')[-1].split(' ')))[0] for val in values[tick+1:tick+ny+1]]
                            j_flags = np.array(j_flags, dtype = int) - 1
                            value = [list(filter(f,val.split(' J= ')[-1].split(' ')))[1:] for val in values[tick+1:tick+ny+1]]
                            value = [ [val[:-1] if val[-1] in ['p','i'] else val for val in vals] for vals in value]
                            value = np.array(value, dtype = float)

                            for j_, val_ in zip(j_flags,value):
                                grid_property[t_flag, k_flag, j_, list(i_flags)[0]:list(i_flags)[-1]+1] = val_
                        
                t_flag += 1

    Oil_sat = grid_property  
    ################### WSAT ##########################   
    grid_property = - np.ones((len(Ptime)+1,nz,ny,nx))
    t_flag = 0

    for i in range(len(lines)):
        # Read Line:
        line = lines[i].split()
        ## Read grid property:
        if [''.join(line)][0] == 'WaterSaturation(fraction)':
            start = i + 2
            for j in range(start,len(lines)):
                if '*' * 70 in lines[j].split(' '):
                    end = j
                    break
                elif '*' * 131 + '\n' in lines[j].split(' '):
                    end = j
                    break
            val_lines = lines[start:end]

            if not any(['Plane K = ' in var for var in val_lines]):
                value = [''.join(val_lines[np.where(['All values are' in val for val in val_lines])[0][0]].split(' '))][0]
                value = value.split('Allvaluesare')[1].split('\n')[0]
                grid_property[t_flag] = value 
                t_flag += 1
                
            else:
                value_line_lst = list(np.where(['Plane K =' in  val for val in val_lines])[0])
                value_line_lst.append(len(val_lines))
                for v in range(len(value_line_lst)-1):
                    s, e = value_line_lst[v], value_line_lst[v+1]
                    values =  val_lines[s:e]
                    k_flag = int(val_lines[s:e][0].split(' Plane K = ')[-1].split(' ')[0]) - 1 
                    if any(['All values are' in val for val in values]):
                        value = [''.join(values[np.where(['All values are' in val for val in values])[0][0]].split(' '))][0]
                        value = value.split('Allvaluesare')[1].split('\n')[0]
                        grid_property[t_flag, k_flag] = float(value)

                    else:
                        i_index = np.where(['    I = ' in val for val in values])[0]
                        for tick in i_index: 
                            i_flags = np.array(values[tick].split('  I = ')[-1].split('     '),dtype = int) -1
                            f = lambda x: x if not x in ['', '\n'] else None
                            j_flags = [list(filter(f,val.split(' J= ')[-1].split(' ')))[0] for val in values[tick+1:tick+ny+1]]
                            j_flags = np.array(j_flags, dtype = int) - 1
                            value = [list(filter(f,val.split(' J= ')[-1].split(' ')))[1:] for val in values[tick+1:tick+ny+1]]
                            value = [ [val[:-1] if val[-1] in ['p','i'] else val for val in vals] for vals in value]
                            value = np.array(value, dtype = float)

                            for j_, val_ in zip(j_flags,value):
                                grid_property[t_flag, k_flag, j_, list(i_flags)[0]:list(i_flags)[-1]+1] = val_
                        
                t_flag += 1

    Water_sat = grid_property  

    ################### PRESSURE ##########################   
    grid_property = - np.ones((len(Ptime)+1,nz,ny,nx))
    t_flag = 0

    for i in range(len(lines)):
        # Read Line:
        line = lines[i].split()
        ## Read grid property:
        if [''.join(line)][0] == 'Pressure(psi)':
            start = i + 2
            for j in range(start,len(lines)):
                if '*' * 70 in lines[j].split(' '):
                    end = j
                    break
                elif '*' * 131 + '\n' in lines[j].split(' '):
                    end = j
                    break
            val_lines = lines[start:end]

            if not any(['Plane K = ' in var for var in val_lines]):
                value = [''.join(val_lines[np.where(['All values are' in val for val in val_lines])[0][0]].split(' '))][0]
                value = value.split('Allvaluesare')[1].split('\n')[0]
                grid_property[t_flag] = value 
                t_flag += 1
                
            else:
                value_line_lst = list(np.where(['Plane K =' in  val for val in val_lines])[0])
                value_line_lst.append(len(val_lines))
                for v in range(len(value_line_lst)-1):
                    s, e = value_line_lst[v], value_line_lst[v+1]
                    values =  val_lines[s:e]
                    k_flag = int(val_lines[s:e][0].split(' Plane K = ')[-1].split(' ')[0]) - 1 
                    if any(['All values are' in val for val in values]):
                        value = [''.join(values[np.where(['All values are' in val for val in values])[0][0]].split(' '))][0]
                        value = value.split('Allvaluesare')[1].split('\n')[0]
                        grid_property[t_flag, k_flag] = float(value)

                    else:
                        i_index = np.where(['    I = ' in val for val in values])[0]
                        for tick in i_index: 
                            i_flags = np.array(values[tick].split('  I = ')[-1].split('     '),dtype = int) -1
                            f = lambda x: x if not x in ['', '\n'] else None
                            j_flags = [list(filter(f,val.split(' J= ')[-1].split(' ')))[0] for val in values[tick+1:tick+ny+1]]
                            j_flags = np.array(j_flags, dtype = int) - 1
                            value = [list(filter(f,val.split(' J= ')[-1].split(' ')))[1:] for val in values[tick+1:tick+ny+1]]
                            value = [ [val[:-1] if val[-1] in ['p','i'] else val for val in vals] for vals in value]
                            value = np.array(value, dtype = float)

                            for j_, val_ in zip(j_flags,value):
                                grid_property[t_flag, k_flag, j_, list(i_flags)[0]:list(i_flags)[-1]+1] = val_
                        
                t_flag += 1

    Pressure = grid_property
        
    return Oil_sat, Water_sat, Pressure

def read_out(tstep,filename = 'Forward_simulation.out',num_prod=4, num_inj=1,) :
    f = open(filename, 'r')
    f_ = f.readlines()
    f.close()
    
    num_tstep = tstep.shape[0]
    
    WBHP = np.zeros((num_tstep, num_prod+num_inj))
    WGBP = np.zeros((num_tstep, num_prod+num_inj))
    
    WOPR = np.zeros((num_tstep, num_prod))
    WWPR = np.zeros((num_tstep, num_prod))
    WGPR = np.zeros((num_tstep, num_prod))
    
    WCOP = np.zeros((num_tstep, num_prod))
    WCWP = np.zeros((num_tstep, num_prod))
    WCGP = np.zeros((num_tstep, num_prod))
    
    WWIR = np.zeros((num_tstep, num_inj))
    WCWI = np.zeros((num_tstep, num_inj))
    
    flag = 0
    t_flag = 0
     
    for i in range(len(f_)):
    
        line = f_[i].split()
    
        if len(line)>20 and flag == 0:
            if [''.join(line[3 : -2])][0] == 'IMEXFIELDSUMMARY':
                flag = 1;
                
        if flag == 1:
            WBHP[t_flag, 0:4] = np.array(f_[i+17].split('+')[1:5],dtype = float)   #WBHP, psia
            WGBP[t_flag, 0:4] = np.array(f_[i+18].split('+')[1:5],dtype = float)   #WGBP, psia
            WOPR[t_flag, 0:4] = np.array(f_[i+27].split('+')[1:5],dtype = float)   #WOPR, STB/day
            WWPR[t_flag, 0:4] = np.array(f_[i+28].split('+')[1:5],dtype = float)   #WWPR
            WGPR[t_flag, 0:4] = np.array(f_[i+29].split('+')[1:5],dtype = float)   #WGPR, MSCF/day (Disolved Gas)
            WCOP[t_flag, 0:4] = np.array(f_[i+36].split('+')[1:5],dtype = float)   #WCOP, MSTB
            WCWP[t_flag, 0:4] = np.array(f_[i+37].split('+')[1:5],dtype = float)   #WCWP
            WCGP[t_flag, 0:4] = np.array(f_[i+38].split('+')[1:5],dtype = float)   #WCGP, MMSCF
            
           
            WBHP[t_flag, 4] = np.array(f_[i+82].split('+')[1],dtype = float)   #WBHP, psia
            WGBP[t_flag, 4] = np.array(f_[i+83].split('+')[1],dtype = float)            
            WWIR[t_flag,0] = np.array(f_[i+125].split('+')[1],dtype = float)   #WWIR, STB/day
            WCWI[t_flag,0] = np.array(f_[i+130].split('+')[1],dtype = float)   #WCWI, MSTB       
            
            t_flag += 1;
            flag = 0;
    # return tstep, WBHP, WGBP, WOPR, WWPR, WGPR, WCOP, WCWP, WCGP         

    return WBHP, WGBP, WOPR, WWPR, WGPR, WCOP, WCWP, WCGP, WWIR, WCWI       


def plot_Result(x, tstep = np.arange(50, 2001, 50), x_lab='NA', y_lab='NA',title='NA'):
    num = x.shape[-1]
    if num != 1:
        for i in range(num) :
            plt.plot(tstep,x[:,i], linewidth=3, label = 'Well #' + str(i+1))
        if x_lab != 'NA': plt.xlabel(x_lab)
        if y_lab != 'NA': plt.ylabel(y_lab)        
        if title != 'NA': plt.title(y_lab)        
        plt.legend()    
    else:
        plt.plot(tstep,x, linewidth=3)
        if x_lab != 'NA': plt.xlabel(x_lab)
        if y_lab != 'NA': plt.ylabel(y_lab)        
        if title != 'NA': plt.title(y_lab)        
        plt.legend()    
        
        
def load_real(filename = '_Real_Rulebased_Model.dat',v_min = 0, v_max = 1) :
    f=open(filename,'r')
    Real_res = f.read().split(",")
    f.close()
    Real_res = Real_res[:15680]
    Real_Res = np.asarray(Real_res, dtype=float).reshape(20,28,28)
    Real_Res = np.moveaxis(Real_Res, 2,0)
    Real_Res = np.moveaxis(Real_Res, 2,1)
    Real_Res = (Real_Res -Real_Res.min())/(Real_Res.max()-Real_Res.min())*2-1
    Real_Res_=Real_Res.reshape(1,28,28,20)
    Real_Res_[Real_Res_<-0.25] = -0.25
    Real_Res_var = (v_max - v_min) * ( Real_Res_ - Real_Res_.min() ) / (Real_Res_.max() - Real_Res_.min()) + v_min
    return Real_Res_var


def Tstep(time,filename = 'TSTEP.DAT'):
    f = open(filename, 'w')
    for i,j in enumerate(time):
        f.write('*TIME ' + str(j) +'\n')
    f.close()
    
def SetPerm(Real_Res_Perm, dir_ = 'null',nx = 28, ny = 28, nz = 20) :
    Real_Res_Perm = np.exp(Real_Res_Perm)
    if dir_ == 'null':
        f = open('PERM_X.DAT', 'w')
    else:
        f = open(dir_+'\\PERM_X.DAT', 'w')
    f.write('*PERMI *ALL \n')
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                f.write(str(Real_Res_Perm[i,j,k]) + '\n');    
    f.close()

    if dir_ == 'null':
        f = open('PERM_Y.DAT', 'w')
    else:
        f = open(dir_+'\\PERM_Y.DAT', 'w')
    f.write('*PERMJ *ALL \n')
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                f.write(str(Real_Res_Perm[i,j,k]) + '\n');
    f.close()

    if dir_ == 'null':
        f = open('PERM_Z.DAT', 'w')
    else:
        f = open(dir_+'\\PERM_Z.DAT', 'w')	
    f.write('*PERMK *ALL \n')
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                f.write(str(0.1*Real_Res_Perm[i,j,k]) + '\n');
    f.close()
    
def SetPoro(Real_Res_Poro, dir_ = 'null',nx = 28, ny = 28, nz = 20) :
    if dir_ == 'null':
        f = open('PORO.DAT', 'w')
    else:
        f = open(dir_+'\\PORO.DAT', 'w')
    f.write('*POR *ALL \n')
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                f.write(str(Real_Res_Poro[i,j,k]) + '\n');    
    f.close()

import matplotlib.cm as cm
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

def visualize_3D_full(cube, nx,ny,nz,angle=320, cmap_ = cm.afmhot):
    ### Note that cube should numpy array whose dimension is NX x NY x NZ
    cube = (cube-cube.min())/(cube.max()-cube.min())
    facecolors = cmap_(cube)
    facecolors[:,:,:,-1] = 1
    # this is to select where to colour
    fc = np.zeros((nx,ny,nz))
    fc[:,0,:] = 1
    fc[-1,:,:] = 1
    fc[:,:,-1] = 1
    facecolors[:,:,:,-1] = fc
    facecolors = explode(facecolors)

    filled = facecolors[:,:,:,-1] != 0
    x, y, z = expand_coordinates(np.indices(np.array(filled.shape) + 1))

    fig = plt.figure(figsize=(15,7))
    ax = plt.subplot(projection='3d')
    ax.view_init(30, angle)
    ax.set_xlim(right=nx*2)
    ax.set_ylim(top=ny*2)
    ax.set_zlim(top=nz*2)      
    im = ax.voxels(x, y, z, filled, facecolors=facecolors)

    ax.set_xlabel("X-axis, grid")
    ax.set_ylabel("Y-axis, grid")
    ax.set_zlabel("Z-axis, grid")
    ax.xaxis.set_pane_color([1.0,1.0,1.0,1.0])    # fig.colorbar(im, ax=ax)
    ax.yaxis.set_pane_color([1.0,1.0,1.0,1.0])    # fig.colorbar(im, ax=ax)
    ax.zaxis.set_pane_color([1.0,1.0,1.0,1.0])    # fig.colorbar(im, ax=ax)
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,1)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,1)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,1)
    
def visualize_3D_fence_diagram(cube, nx,ny,nz,angle=320, cmap_ = cm.afmhot):
    ### Note that cube should numpy array whose dimension is NX x NY x NZ
    cube = (cube-cube.min())/(cube.max()-cube.min())
    facecolors = cmap_(cube)
    
    # this is to select where to colour
    fc = np.zeros((nx,ny,nz))
    fc[:,-3,:] = 1
    fc[3,:,:] = 1
    fc[:,:,3] = 1
    facecolors[:,:,:,-1] = fc
    
    facecolors = explode(facecolors)
    filled = facecolors[:,:,:,-1] != 0
    x, y, z = expand_coordinates(np.indices(np.array(filled.shape) + 1))

    fig = plt.figure(figsize=(15,7))
    ax = plt.subplot(projection='3d')
    ax.view_init(30, angle)
    ax.set_xlim(right=nx*2)
    ax.set_ylim(top=ny*2)
    ax.set_zlim(top=nz*2)      
    im = ax.voxels(x, y, z, filled, facecolors=facecolors)

    ax.set_xlabel("X-axis, grid")
    ax.set_ylabel("Y-axis, grid")
    ax.set_zlabel("Z-axis, grid")
    ax.xaxis.set_pane_color([1.0,1.0,1.0,1.0])    # fig.colorbar(im, ax=ax)
    ax.yaxis.set_pane_color([1.0,1.0,1.0,1.0])    # fig.colorbar(im, ax=ax)
    ax.zaxis.set_pane_color([1.0,1.0,1.0,1.0])    # fig.colorbar(im, ax=ax)
    
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,1)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,1)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,1)
def explode(data):
    shape_arr = np.array(data.shape)
    size = shape_arr[:3]*2 - 1
    exploded = np.zeros(np.concatenate([size, shape_arr[3:]]), dtype=data.dtype)
    exploded[::2, ::2, ::2] = data
    return exploded

def expand_coordinates(indices):
    x, y, z = indices
    x[1::2, :, :] += 1
    y[:, 1::2, :] += 1
    z[:, :, 1::2] += 1
    return x, y, z         


## 2D Cross sectional area    
def plot_h_cs(Gen, n, save = False, colorbar = False, vmin = -1, vmax = 1):
    fig = plt.figure()
    ax = fig.gca()        
    plt.imshow(Gen[:,:,n].T,vmin=vmin,vmax=vmax,cmap='afmhot')

    # plt.imshow(Gen[:,:,n],vmin=vmin,vmax=vmax,cmap='afmhot',interpolation = 'kaiser')
    
	#fig.colorbar(im,ticks = np.array([0,0.25,0.5,0.75,1]))
    #plt.clim(0,1)
    ax.set_xlabel("X-axis, km")
    ax.set_ylabel("Y-axis, km")
    ax.set_xticks(np.array(range(0,8))*4-0.5)
    ax.set_yticks(np.array(range(0,8))*4-0.5)
    ax.set_xticklabels(['0','0.4','0.8','1.2','1.6','2.0','2.4','2.8'])
    ax.set_yticklabels(['0','0.4','0.8','1.2','1.6','2.0','2.4','2.8'])

    if save != False:
        fig.savefig("2D_plot_horizontal_" + str(n)+".png")     


def plot_v_cs(Gen, n, save = False, colorbar = False, vmin = -1, vmax = 1):
    fig = plt.figure()
    ax = fig.gca()    
    im= plt.imshow(np.flipud(Gen[:,n,:].T),vmin=vmin,vmax=vmax,cmap='afmhot')

    # im= plt.imshow(np.flipud(Gen[:,n,:].T),vmin=vmin,vmax=vmax,cmap='afmhot',interpolation = 'kaiser')

    if colorbar :
        fig.colorbar(im,ticks = np.array([0,0.25,0.5,0.75,1]))
        plt.clim(0,1)

    ax.set_xlabel("X-axis, km")
    ax.set_ylabel("Z-axis, m")
    ax.xaxis
    ax.set_xticks(np.array(range(0,8))*4-0.5)
    ax.set_yticks(np.array(range(0,6))*4-0.5)
    ax.set_xticklabels(['0','0.4','0.8','1.2','1.6','2.0','2.4','2.8'])
    ax.set_yticklabels(['10','8','6','4','2','0'])    
    if save != False:
        fig.savefig("2D_plot_vertical_" + str(n)+".png")     



## 2D Cross sectional area    
def plot_h_cs_(Gen, Mask, n, save = False, colorbar = False, vmin = -1, vmax = 1):
    tem = cm.afmhot(Gen)
    tem[:,:,:,-1] = Mask
    
    fig = plt.figure()
    ax = fig.gca()

    plt.imshow(tem[:,:,n,:], vmin=vmin,vmax=vmax)

    ax.set_xlabel("Y-axis, km")
    ax.set_ylabel("X-axis, km")
    ax.set_xticks(np.array(range(0,8))*4-0.5)
    ax.set_yticks(np.array(range(0,8))*4-0.5)
    ax.set_xticklabels(['0','0.4','0.8','1.2','1.6','2.0','2.4','2.8'])
    ax.set_yticklabels(['0','0.4','0.8','1.2','1.6','2.0','2.4','2.8'])


def plot_v_cs_(Gen, Mask,n, save = False, colorbar = False, vmin = -1, vmax = 1):
    tem = cm.afmhot(Gen)
    tem[:,:,:,-1] = Mask

    fig = plt.figure()
    ax = fig.gca()    
    tem = np.moveaxis(tem,2,0)
    tem = np.moveaxis(tem,1,2)
    plt.imshow(tem[::-1,n,:,:], vmin=vmin,vmax=vmax)

    # im= plt.imshow(np.flipud(Gen[:,n,:].T),vmin=vmin,vmax=vmax,cmap='afmhot',interpolation = 'kaiser')

    if colorbar :
        fig.colorbar(im,ticks = np.array([0,0.25,0.5,0.75,1]))
        plt.clim(0,1)

    ax.set_xlabel("X-axis, km")
    ax.set_ylabel("Z-axis, m")
    ax.xaxis
    ax.set_xticks(np.array(range(0,8))*4-0.5)
    ax.set_yticks(np.array(range(0,6))*4-0.5)
    ax.set_xticklabels(['0','0.4','0.8','1.2','1.6','2.0','2.4','2.8'])
    ax.set_yticklabels(['10','8','6','4','2','0'])    
    if save != False:
        fig.savefig("2D_plot_vertical_" + str(n)+".png")     

def null_grid(num_en = 100, nx = 100, ny = 100, nz = 10, tstep = np.array(range(100,3001,100))):
    return np.zeros((num_en, len(tstep), nz, ny, nx)), np.zeros((num_en, len(tstep), nz, ny, nx)), np.zeros((num_en, len(tstep), nz, ny, nx))

def null_prod(num_en = 100, num_prod=4, num_inj=1, tstep = np.array(range(100,3001,100))) :

    num_tstep = tstep.shape[0]
    WBHP = np.zeros((num_en,num_tstep, num_prod+num_inj))
    WGBP = np.zeros((num_en,num_tstep, num_prod+num_inj))
    
    WOPR = np.zeros((num_en,num_tstep, num_prod))
    WWPR = np.zeros((num_en,num_tstep, num_prod))
    WGPR = np.zeros((num_en,num_tstep, num_prod))
    
    WCOP = np.zeros((num_en,num_tstep, num_prod))
    WCWP = np.zeros((num_en,num_tstep, num_prod))
    WCGP = np.zeros((num_en,num_tstep, num_prod))
    
    WWIR = np.zeros((num_en,num_tstep, num_inj))
    WCWI = np.zeros((num_en,num_tstep, num_inj))
    return WBHP, WGBP, WOPR, WWPR, WGPR, WCOP, WCWP, WCGP, WWIR, WCWI

# def read_out(filename = 'HG_3D.out',num_prod=4, num_inj=1,tstep = np.array(range(100,3001,100))) :
#     f = open(filename, 'r')
#     f_ = f.readlines()
#     f.close()
    
#     num_tstep = tstep.shape[0]
    
#     WBHP = np.zeros((num_tstep, num_prod+num_inj))
#     WGBP = np.zeros((num_tstep, num_prod+num_inj))
    
#     WOPR = np.zeros((num_tstep, num_prod))
#     WWPR = np.zeros((num_tstep, num_prod))
#     WGPR = np.zeros((num_tstep, num_prod))
    
#     WCOP = np.zeros((num_tstep, num_prod))
#     WCWP = np.zeros((num_tstep, num_prod))
#     WCGP = np.zeros((num_tstep, num_prod))
    
#     WWIR = np.zeros((num_tstep, num_inj))
#     WCWI = np.zeros((num_tstep, num_inj))
    
#     flag = 0
#     t_flag = 0
     
#     for i in range(len(f_)):
    
#         line = f_[i].split()
    
#         if len(line)>20 and flag == 0:
#             if [''.join(line[3 : -2])][0] == 'IMEXFIELDSUMMARY':
#                 flag = 1;
                
#         if flag == 1:
#             WBHP[t_flag, 0:4] = np.array(f_[i+17].split('+')[1:5],dtype = float)   #WBHP, psia
#             WGBP[t_flag, 0:4] = np.array(f_[i+18].split('+')[1:5],dtype = float)     #WGBP, psia
#             WOPR[t_flag, 0:4] = np.array(f_[i+27].split('+')[1:5],dtype = float)   #WOPR, STB/day
#             WWPR[t_flag, 0:4] = np.array(f_[i+28].split('+')[1:5],dtype = float)   #WWPR
#             WGPR[t_flag, 0:4] = np.array(f_[i+29].split('+')[1:5],dtype = float)   #WGPR, MSCF/day (Disolved Gas)
#             WCOP[t_flag, 0:4] = np.array(f_[i+36].split('+')[1:5],dtype = float)   #WCOP, MSTB
#             WCWP[t_flag, 0:4] = np.array(f_[i+37].split('+')[1:5],dtype = float)   #WCWP
#             WCGP[t_flag, 0:4] = np.array(f_[i+38].split('+')[1:5],dtype = float)   #WCGP, MMSCF
            
           
#         #    WBHP[t_flag, 5] = np.array(f_[i+82].split('+')[1],dtype = float)   #WWIR, STB/day
#         #    WCWI[t_flag, 5] = np.array(f_[i+83].split('+')[1],dtype = float)    #WCWI, MSTB
            
#             WWIR[t_flag,0] = np.array(f_[i+125].split('+')[1],dtype = float)   #WWIR, STB/day
#             WCWI[t_flag,0] = np.array(f_[i+130].split('+')[1],dtype = float)   #WCWI, MSTB       
            
#             t_flag += 1;
#             flag = 0;
            
#     return WBHP, WGBP, WOPR, WWPR, WGPR, WCOP, WCWP, WCGP, WWIR, WCWI  


def plot_Result(x, tstep, x_lab='NA', y_lab='NA',title='NA'):
    num = x.shape[1]
    for i in range(num) :
        plt.plot(tstep,x[:,i], linewidth=3, label = 'Well #' + str(i+1))
    if x_lab != 'NA': plt.xlabel(x_lab)
    if y_lab != 'NA': plt.ylabel(y_lab)        
    if title != 'NA': plt.title(y_lab)        
    plt.legend()        