import numpy as np
from skeleton import *
from observables import *
import matplotlib.pyplot as plt
from time import time
import os

run_simulation = False
#-----------  INPUT  ----------
#box_dim = 1
#------Thermodynamic input-----
density = 0.5
temp = 1.36
#------Simulation input -------
num_tsteps = 10000
timestep =  0.001/(density*temp)
d_l_ratio = 3
#------------------------------
num_atoms = 4*d_l_ratio**3
box_dim = (num_atoms/density)**(1/3)
volume = box_dim**3

simulation_setting = {'box_dim' : box_dim, 'num_tsteps': num_tsteps, 'density' : density, 'temp' : temp , 'd_l_ratio' : d_l_ratio, 'timestep' : timestep , 'num_atoms' : num_atoms, 'volume' : volume}
experiment_name = f'T{temp}_density{density}_numsteps{num_tsteps}_numatoms{num_atoms}'
experiment_directory=os.path.join('Data', experiment_name)

print('The parameters of the simulation are:')
print(simulation_setting)

if(run_simulation):
    # Run it and get the data
    init_pos = np.array(fcc_lattice(d_l_ratio,box_dim))
    init_vel = init_velocity(num_atoms, temp)
    tic = time()
    data_pos, data_vel, steps_before_equib, avg_temp = simulate(init_pos, init_vel, temp, num_tsteps, timestep, box_dim)
    toc = time()
    print(f'It took {toc - tic} seconds to simulate.')

    simulation_setting['steps_before_equib'] = steps_before_equib
    simulation_setting['avg_temp'] = avg_temp

    #if the folder does not exist, create it
    if not os.path.exists(experiment_directory):
        os.makedirs(experiment_directory)

    np.save(f'{experiment_directory}/pos', data_pos)
    np.save(f'{experiment_directory}/vel', data_vel)
    np.save(f'{experiment_directory}/simulation_setting', simulation_setting)

else:
    # Load the data from drive
    try:
        data_pos = np.load(f'{experiment_directory}/pos.npy')
        data_vel = np.load(f'{experiment_directory}/vel.npy')
        simulation_setting = np.load(f'{experiment_directory}/simulation_setting.npy', allow_pickle='TRUE').item()
        locals().update(simulation_setting)  #load data from the dictionary and update local variables
    except:
        print(f'No directory found with the given system parameters.')

post_equi_pos = data_pos[steps_before_equib:]
post_equi_vel = data_vel[steps_before_equib:]


plot_energy(data_pos,data_vel,simulation_setting)


tic = time()
p_value, p_array = pressure(post_equi_pos, simulation_setting)
print(f'Pressure is {p_value}')
print(f'Compressibility factor is {p_value/(temp*density)}')

toc = time()
print(f'It took {toc - tic} seconds to calculate pressure')

tic = time()
tau, err = auto_corr_func(p_array)
print(f'The correlation time is {tau*(timestep)}, the pressure error is {err}')
print(f'The compressibility factor error is {err/(temp*density)}')
toc  = time()
print(f'It took {toc - tic} seconds to calculate the auto correlation function')


tic = time()
disc_steps=100
g = pair_corr_func(post_equi_pos, simulation_setting, disc_steps=disc_steps)
toc = time()
print(f'It took {toc - tic} seconds to calculate the pair correlation function.')

plt.plot(np.linspace(0,box_dim,disc_steps), g, '-b')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.show()
