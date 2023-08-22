import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def read_out_grid(tstep, filename = 'Forward_simulation.out', nx=64, ny=64, ini_sw = 0.156, ini_press = 3100):
    f = open(filename, 'r')
    f_ = f.readlines()
    f.close()
    
    num_tstep = tstep.shape[0]
    
    
    Pressure = np.zeros((num_tstep+1,ny,nx)); Pressure[0] = ini_press
    Water_sat = np.zeros((num_tstep+1,ny,nx)); Water_sat[0] = ini_sw
    Oil_sat = np.zeros((num_tstep+1,ny,nx)); Oil_sat[0] = 1 - ini_sw
    
    flag = -1
    t_flag = 1
    
    for i in range(len(f_)):
        # Read Line:
        line = f_[i].split()
        ## Read Oil saturation:
        if flag == -1:
            if [''.join(line)][0] == 'OilSaturation(fraction)':
                flag += 1;
        elif flag == 0:
            if [''.join(line)][0] == 'OilSaturation(fraction)':
                flag += 1;
        
        if flag == 1:
            temp = f_[i+3:i+3+330]
            grid_list =[]
            for grid_line in range(len(temp)):
                if temp[grid_line].split() == []:
                    continue
                elif temp[grid_line].split()[0] == 'I': 
                    temp_ = temp[grid_line+1:grid_line+1+64]
                    temp_list = []
                    for grid_line_ in temp_:
                        grid_line_ = grid_line_.split()
                        grid_line_ = list(map(lambda x: x.replace("p", ""), grid_line_))
                        grid_line_ = list(map(lambda x: x.replace("i", ""), grid_line_))
                        grid_line_ = np.array(grid_line_)[2:].reshape(1,-1)
                        temp_list.append(grid_line_)
                    grid_list.append(np.array(temp_list,dtype = float).squeeze()[::-1,:])
            
            for flg in range(4):
                Oil_sat[t_flag,:,0+14*flg:14+14*flg] = grid_list[flg]
            Oil_sat[t_flag,:,-8:] = grid_list[-1]
            t_flag += 1;
            flag = 0;       

    flag = -1
    t_flag = 1
    for i in range(len(f_)):
        # Read Line:
        line = f_[i].split()
        ## Read Water saturation:
        if flag == -1:
            if [''.join(line)][0] == 'WaterSaturation(fraction)':
                flag += 1;
        elif flag == 0:
            if [''.join(line)][0] == 'WaterSaturation(fraction)':
                flag += 1;
        
        if flag == 1:
            temp = f_[i+3:i+3+330]
            grid_list =[]
            for grid_line in range(len(temp)):
                if temp[grid_line].split() == []:
                    continue
                elif temp[grid_line].split()[0] == 'I': 
                    temp_ = temp[grid_line+1:grid_line+1+64]
                    temp_list = []
                    for grid_line_ in temp_:
                        grid_line_ = grid_line_.split()
                        grid_line_ = list(map(lambda x: x.replace("p", ""), grid_line_))
                        grid_line_ = list(map(lambda x: x.replace("i", ""), grid_line_))
                        grid_line_ = np.array(grid_line_)[2:].reshape(1,-1)
                        temp_list.append(grid_line_)
                    grid_list.append(np.array(temp_list,dtype = float).squeeze()[::-1,:])
            
            for flg in range(4):
                Water_sat[t_flag,:,0+14*flg:14+14*flg] = grid_list[flg]
            Water_sat[t_flag,:,-8:] = grid_list[-1]
            t_flag += 1;
            flag = 0;
            
    flag = -1
    t_flag = 1
    for i in range(len(f_)):
        # Read Line:
        line = f_[i].split()
        ## Read Pressure:
        if flag == -1:
            if [''.join(line)][0] == 'Pressure(psia)':
                flag += 1;
        elif flag == 0:
            if [''.join(line)][0] == 'Pressure(psia)':
                flag += 1;
        
        if flag == 1:
            temp = f_[i+3:i+3+330]
            grid_list =[]
            for grid_line in range(len(temp)):
                if temp[grid_line].split() == []:
                    continue
                elif temp[grid_line].split()[0] == 'I': 
                    temp_ = temp[grid_line+1:grid_line+1+64]
                    temp_list = []
                    for grid_line_ in temp_:
                        grid_line_ = grid_line_.split()
                        grid_line_ = list(map(lambda x: x.replace("p", ""), grid_line_))
                        grid_line_ = list(map(lambda x: x.replace("i", ""), grid_line_))
                        grid_line_ = np.array(grid_line_)[2:].reshape(1,-1)
                        temp_list.append(grid_line_)
                    grid_list.append(np.array(temp_list,dtype = float).squeeze()[::-1,:])
            
            for flg in range(4):
                Pressure[t_flag,:,0+14*flg:14+14*flg] = grid_list[flg]
            Pressure[t_flag,:,-8:] = grid_list[-1]
            t_flag += 1;
            flag = 0;
            
    return Oil_sat, Water_sat, Pressure

def null_prod(tstep, num_en = 1, num_iter = 1, num_prod=4, num_inj=1) :
    num_tstep = tstep.shape[0]

    if num_en ==1:
        WBHP = np.zeros((num_en,num_tstep, num_prod+num_inj))
        WGBP = np.zeros((num_en,num_tstep, num_prod+num_inj))
        
        WOPR = np.zeros((num_tstep, num_prod))
        WWPR = np.zeros((num_tstep, num_prod))
        WGPR = np.zeros((num_tstep, num_prod))
        
        WCOP = np.zeros((num_tstep, num_prod))
        WCWP = np.zeros((num_tstep, num_prod))
        WCGP = np.zeros((num_tstep, num_prod))
        
        WWIR = np.zeros((num_tstep, num_inj))
        WCWI = np.zeros((num_tstep, num_inj))        
    elif num_iter==1:
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
    else:
        WBHP = np.zeros((num_iter, num_en, num_tstep, num_prod+num_inj))
        WGBP = np.zeros((num_iter, num_en, num_tstep, num_prod+num_inj))
        
        WOPR = np.zeros((num_iter, num_en, num_tstep, num_prod))
        WWPR = np.zeros((num_iter, num_en, num_tstep, num_prod))
        WGPR = np.zeros((num_iter, num_en, num_tstep, num_prod))
        
        WCOP = np.zeros((num_iter, num_en, num_tstep, num_prod))
        WCWP = np.zeros((num_iter, num_en, num_tstep, num_prod))
        WCGP = np.zeros((num_iter, num_en, num_tstep, num_prod))
        
        WWIR = np.zeros((num_iter, num_en, num_tstep, num_inj))
        WCWI = np.zeros((num_iter, num_en, num_tstep, num_inj))        
    # return WBHP, WGBP, WOPR, WWPR, WGPR, WCOP, WCWP, WCGP, WWIR, WCWI
    return WBHP, WOPR, WWPR, WGPR, WCOP, WCWP, WCGP, WWIR, WCWI

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