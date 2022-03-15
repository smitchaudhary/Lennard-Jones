import numpy as np
from skeleton import *
from observables import *
import matplotlib.pyplot as plt
from time import time
import os

run_simulation = True
#-----------  INPUT  ----------
#------Thermodynamic input-----
density = 1
temp = np.linspace(0.5,2.5,21)

#------Simulation input -------
num_tsteps = 10000
timestep =  0.001
d_l_ratio = 2
#------------------------------
num_atoms = 4*d_l_ratio**3
box_dim = (num_atoms/density)**(1/3)
volume = box_dim**3

experiment_name = f'Phase_tr_d{density}_t_{temp[0]}-{temp[-1]}'
experiment_directory=os.path.join('Data', experiment_name)

pressures = np.zeros(temp.shape)
errors = np.zeros(temp.shape)

if run_simulation:
    for i,tem in enumerate(temp):
        # Run it and get the data
        simulation_setting = {'box_dim' : box_dim, 'num_tsteps': num_tsteps, 'density' : density, 'temp' : tem , 'd_l_ratio' : d_l_ratio, 'timestep' : timestep , 'num_atoms' : num_atoms, 'volume' : volume} 
        init_pos = np.array(fcc_lattice(d_l_ratio,box_dim))
        init_vel = init_velocity(num_atoms, tem)

        data_pos, data_vel, steps_before_equib, avg_temp = simulate(init_pos, init_vel, tem, num_tsteps, timestep, box_dim)

        simulation_setting['steps_before_equib'] = steps_before_equib
        #if the folder does not exist, create it
        post_equi_pos = data_pos[steps_before_equib:]
        post_equi_vel = data_vel[steps_before_equib:]

        p_value, p_array = pressure(post_equi_pos, simulation_setting)
        tau, err = auto_corr_func(p_array)

        temp[i] = avg_temp
        pressures[i] = p_value
        errors[i] = err

    if not os.path.exists(experiment_directory):
        os.makedirs(experiment_directory)
    
    np.save(f'{experiment_directory}/temp_v', temp)
    np.save(f'{experiment_directory}/pressures', pressures)
    np.save(f'{experiment_directory}/simulation_setting', simulation_setting)
    np.save(f'{experiment_directory}/errors', errors)
    
else:
    # Load the data from drive
    try:
        temp = np.load(f'{experiment_directory}/temp_v.npy')
        errors = np.load(f'{experiment_directory}/errors.npy')
        pressures = np.load(f'{experiment_directory}/pressures.npy')
        simulation_setting = np.load(f'{experiment_directory}/simulation_setting.npy', allow_pickle='TRUE').item()
    except:
        print(f'No directory found with the given system parameters.')

if(run_simulation):
    print(f'Max velocity is {np.amax(abs(np.array(data_vel)))}')
    print(f'Max position is {np.amax(abs(np.array(data_pos)))}')
print(simulation_setting['box_dim'])
arr1inds = temp.argsort()
temp = temp[arr1inds[::-1]]
pressures = pressures[arr1inds[::-1]]
errors = errors[arr1inds[::-1]]

#print(temp)
#print(pressures)
#plt.plot(temp, pressures)
#plt.show()
plt.errorbar(temp, pressures, yerr = errors, fmt='o', markersize = 3,  color='green', ecolor='lightgray')
plt.show()

