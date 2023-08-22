# =============================================================================
# Main Script for CMG Parser _ ver 1.0
# Author: Honggeun Jo
# Date: 08/21/2023 - 
# =============================================================================

# =============================================================================
# Import librarays & Initializing
# =============================================================================
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import Subfunction as sub

# Basic info for CMG
CMG =  "C:\\Program Files\\CMG\\IMEX\\2022.10\\Win_x64\\EXE\\mx202210.exe"
Current_Working_Dir = "C:/100 Research/cmg_parser"
os.chdir(Current_Working_Dir)
SIM =  "Forward_simulation.dat"
# Total_Number of realizations
No_realization = 100
# Reservoir extent
nx = 28; ny = 28; nz = 20
# number of wells 
No_well = 4
# time steps
Ptime = np.arange(50, 2501, 50)
## Run CMG : subprocess.call([CMG,"-f",SIM])

# =============================================================================
# Load ensemble
# =============================================================================
Realizations_LogPerm = np.load('ensemble_101.npz')['logPermeability'][1:]
Realizations_Porosity = np.load('ensemble_101.npz')['porosity'][1:]

# =============================================================================
# Make directories for ensemble
# =============================================================================
sub.Tstep(Ptime)
for i in range(No_realization) :
    if os.path.exists(f'realization_{str(i)}')==False: 
        os.mkdir(f'realization_{str(i)}')
    shutil.copy('Forward_simulation.dat', f'realization_{str(i)}/Forward_simulation.dat')
    shutil.copy('TSTEP.DAT', f'realization_{str(i)}/TSTEP.DAT')
    sub.SetPerm(Realizations_LogPerm[i],dir_ = f'realization_{str(i)}')
    sub.SetPoro(Realizations_Porosity[i],dir_ = f'realization_{str(i)}')
    
# =============================================================================
# Run simulation & save results:
# =============================================================================

WBHP_ini, WOPR_ini, WWPR_ini, WGPR_ini, WCOP_ini, WCWP_ini, WCGP_ini, WWIR_ini, WCWI_ini = sub.null_prod(Ptime,100,1,4,1)

for i in range(No_realization) :
    os.chdir(f'{Current_Working_Dir}/realization_'+str(i))
    subprocess.call([CMG,"-f",SIM])
    WBHP_ini[i], WOPR_ini[i], WWPR_ini[i], WGPR_ini[i], WCOP_ini[i], WCWP_ini[i], WCGP_ini[i], WWIR_ini[i], WCWI_ini[i] = sub.read_out(Ptime)
    print(f'Ensemble simulation for realization_{str(i)} is completed')

# =============================================================================
# Save and visualize results
# =============================================================================
# WOPR:
plt.figure(figsize =(10,10))
for i_ in range(1,5):
    well_num = i_ -1
    plt.subplot(2,2,i_)
    plt.plot(Ptime,WOPR_ini[i,:,well_num].T, '-', c = [0.8,0.8,0.8], label = 'Inital ensemble',linewidth=0.5);
    for i in range(1,100):    
        plt.plot(Ptime,WOPR_ini[i,:,well_num].T, '-', c = [0.8,0.8,0.8],linewidth=0.5);
    plt.legend()
    plt.xlabel('Time, days')
    plt.ylabel('Production rate, STB/day')
    
# WWPR:
plt.figure(figsize =(10,10))
for i_ in range(1,5):
    well_num = i_ -1
    plt.subplot(2,2,i_)
    plt.plot(Ptime,WWPR_ini[i,:,well_num].T, '-', c = [0.8,0.8,0.8], label = 'Inital ensemble',linewidth=0.5);
    for i in range(1,100):    
        plt.plot(Ptime,WWPR_ini[i,:,well_num].T, '-', c = [0.8,0.8,0.8],linewidth=0.5);
    plt.legend()
    plt.xlabel('Time, days')
    plt.ylabel('Production rate, STB/day')

WCOP_ini_ = WCOP_ini*2/3
# WCOP:
plt.figure(figsize =(10,10))
for i_ in range(1,5):
    well_num = i_ -1
    plt.subplot(2,2,i_)
    plt.plot(Ptime,WCOP_ini_[i,:,well_num].T, '-', c = [0.8,0.8,0.8], label = 'Inital ensemble',linewidth=0.5);
    for i in range(1,100):      
        plt.plot(Ptime,WCOP_ini_[i,:,well_num].T, '-', c = [0.8,0.8,0.8],linewidth=0.5);
    plt.legend()
    plt.xlabel('Time, days')
    plt.ylabel('Production amount, MSTB')


# WCWP:
plt.figure(figsize =(10,10))
for i_ in range(1,5):
    well_num = i_ -1
    plt.subplot(2,2,i_)
    plt.plot(Ptime,WCWP_ini[i,:,well_num].T, '-', c = [0.8,0.8,0.8], label = 'Inital ensemble',linewidth=0.5);
    for i in range(1,100):      
        plt.plot(Ptime,WCWP_ini[i,:,well_num].T, '-', c = [0.8,0.8,0.8],linewidth=0.5);
    plt.legend()
    plt.xlabel('Time, days')
    plt.ylabel('Production amount, MSTB')
        
    
# # Horizontal View 
# sub.plot_h_cs(Rule,1,vmin = -2.3,vmax=7)
# for i in range(1,5):
#     sub.plot_h_cs(en_Perm_[0,i-1],1,vmin = -2.3,vmax=7)

# sub.plot_h_cs(en_Perm_fin[0],1,vmin = -2.3,vmax=7)
# sub.plot_h_cs(en_Perm_fin[1],1,vmin = -2.3,vmax=7)
# sub.plot_h_cs(en_Perm_fin[2],1,vmin = -2.3,vmax=7)

# sub.plot_h_cs(en_Perm_[0,0],1,vmin = -2.3,vmax=7)
# sub.plot_h_cs(en_Perm_[0,1],1,vmin = -2.3,vmax=7)
# sub.plot_h_cs(en_Perm_[0,2],1,vmin = -2.3,vmax=7)

# sub.plot_cube_(Rule)
# sub.plot_cube_2(Rule)

# # sub.plot_cube_(en_Perm_[0,0])
# # sub.plot_cube_(en_Perm_[-1,0])

# sub.plot_cube_2(Rule)
# sub.plot_cube_2(en_Perm_[0,0])

# for temp in range(20): 
#     sub.plot_cube_(en_Perm_[0,temp])
#     plt.savefig('Plot_3D_'+str(temp+1))
#     sub.plot_cube_2(en_Perm_[-1,temp])
#     plt.savefig('Plot_CS_'+str(temp+1))

# for temp in range(20): 
#     sub.plot_cube_(en_Perm_[0,temp])
#     plt.savefig('Plot_3D_initial'+str(temp+1))
#     sub.plot_cube_(en_Perm_[-1,temp])
#     plt.savefig('Plot_3D_updated'+str(temp+1))
#     sub.plot_cube_2(en_Perm_[0,temp])
#     plt.savefig('Plot_CS_initial'+str(temp+1))
#     sub.plot_cube_2(en_Perm_[-1,temp])
#     plt.savefig('Plot_CS_updated'+str(temp+1))
    
    
    
# # Get P values:
# temp = WCOP_fin[:,-1,:]
# P_val = np.zeros((4,3))
# for i in range (4):
#     P_val[i,0] = np.percentile(temp[:,i],10)
#     P_val[i,1]= np.percentile(temp[:,i],50)
#     P_val[i,2]= np.percentile(temp[:,i],90)
    
    
# WCOP_r[-1,:]
    
    
    
    
    
    
    
    
    
    