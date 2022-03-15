import numpy as np
from skeleton import *
from observables import *
import matplotlib.pyplot as plt
from time import time
import os

run_simulation = True
#-----------  INPUT  ----------
#------Thermodynamic input-----
density = np.linspace(0.25,0.4,50)
temp = 1.15

#------Simulation input -------
num_tsteps = 10000
timestep =  0.001
d_l_ratio = 3
#------------------------------
num_atoms = 4*d_l_ratio**3


experiment_name = f'Phase_tr_d{temp}_t_{density[0]}-{density[-1]}'
experiment_directory=os.path.join('Data', experiment_name)

tem =np.zeros(density.shape)
pressures = np.zeros(density.shape)
errors = np.zeros(density.shape)

if run_simulation:
    for i,den in enumerate(density):
        # Run it and get the data
        box_dim = (num_atoms/den)**(1/3)
        volume = box_dim**3
        timestep = 0.001/(den*temp)
        simulation_setting = {'box_dim' : box_dim, 'num_tsteps': num_tsteps, 'density' : den, 'temp' : temp , 'd_l_ratio' : d_l_ratio, 'timestep' : timestep , 'num_atoms' : num_atoms, 'volume' : volume} 
        init_pos = np.array(fcc_lattice(d_l_ratio,box_dim))
        init_vel = init_velocity(num_atoms, temp)

        data_pos, data_vel, steps_before_equib, avg_temp = simulate(init_pos, init_vel, temp, num_tsteps, timestep, box_dim)

        simulation_setting['steps_before_equib'] = steps_before_equib
        #if the folder does not exist, create it
        post_equi_pos = data_pos[steps_before_equib:]
        post_equi_vel = data_vel[steps_before_equib:]

        p_value, p_array = pressure(post_equi_pos, simulation_setting)
        tau, err = auto_corr_func(p_array)

        tem[i] = avg_temp
        pressures[i] = p_value
        errors[i] = err

    if not os.path.exists(experiment_directory):
        os.makedirs(experiment_directory)
    print(f'Saving the results...')
    np.save(f'{experiment_directory}/temp_v', temp)
    np.save(f'{experiment_directory}/den_v', density)
    np.save(f'{experiment_directory}/pressures', pressures)
    np.save(f'{experiment_directory}/simulation_setting', simulation_setting)
    np.save(f'{experiment_directory}/errors', errors)
    
else:
    # Load the data from drive
    try:
        tem = np.load(f'{experiment_directory}/temp_v.npy')
        density = np.load(f'{experiment_directory}/den_v.npy')
        errors = np.load(f'{experiment_directory}/errors.npy')
        pressures = np.load(f'{experiment_directory}/pressures.npy')
        simulation_setting = np.load(f'{experiment_directory}/simulation_setting.npy', allow_pickle='TRUE').item()
    except:
        print(f'No directory found with the given system parameters.')


#print(temp)
#print(pressures)
plt.plot(density, pressures/temp)
#plt.show()
#plt.errorbar(density, pressures/temp, yerr = errors/temp, fmt='o', markersize = 1,  color='green', ecolor='lightgray')
#plt.fill_between(density, y - err, y + error)
plt.show()

