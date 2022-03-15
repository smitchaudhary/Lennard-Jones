from skeleton import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def function_to_fit(t, c, tau):
    return c*np.exp(-t/tau)

def auto_corr_func(A):
    """
    Calculates the auto-correlation function.

    parameters
    ----------
    A : np.array
        Observable value at each timestep.

    returns
    ------
    x : np.ndarray
        The value of the auto correlation function at each timestep.
    """
    num_timesteps = len(A)

    x = np.zeros(num_timesteps-1)

    t_array = np.arange(num_timesteps-1)

    for t in t_array:
        A_t = A[t:] #len=N_t

        N_t = num_timesteps - t

        sum_A  = np.sum(A[:N_t])
        sum_A2 = np.sum((A*A)[:N_t])

        sum_A_t = np.sum(A_t[:N_t])

        numerator = N_t * np.sum( A_t*( A[:N_t] ) ) - sum_A * sum_A_t

        denominator = np.sqrt( N_t * sum_A2 - sum_A**2 ) * np.sqrt( N_t * np.sum(A_t*A_t) - sum_A_t**2 )

        x[t] = numerator/denominator

    #Take the well behaved part. By our convention 3/5 of our data.
    num_wb = int(num_timesteps*3/5)
    x = x[:num_wb]
    t_array = t_array[:num_wb]

    opt_params, params_cov = curve_fit(function_to_fit, t_array, x)
    tau = opt_params[1]
    ideal = [function_to_fit(t, 1, tau) for t in t_array]

    error = np.sqrt( (2*tau/num_timesteps)*(np.average(A[:num_wb]*A[:num_wb]) - np.average(A[:num_wb])**2  ) )
    return tau, error

def pressure(pos, simulation_setting):
    """
    Calculates the Pressure

    parameters
    ----------
    pos : list
        The positions of the particles in cartesian space for every time step
    simulation_setting : dictionary
        Contains all the system parameters

    returns
    ------
    p_value : float
        The value of the pressure at each timestep.
    p_array : np.ndarray
        The value of the pressure at each timestep in all directions.
    """

    box_dim = simulation_setting['box_dim']
    temp = simulation_setting['temp']
    num_atoms = simulation_setting['num_atoms']
    density = simulation_setting['density']
    num_timesteps = len(pos)

    p_array = np.zeros(num_timesteps)
    pair_force = lambda r : 4*( 6*(1/(r))**7 - 12*(1/r)**13) #pair_force=du/dr

    for t in range(num_timesteps):

        _, rel_dist = atomic_distances(pos[t], box_dim)

        np.fill_diagonal(rel_dist, val=2**(1/6))
        pair_f = 4*( 6*(1/(rel_dist))**7 - 12*(1/rel_dist)**13)
        summation = 1/2*np.sum(np.multiply(rel_dist,pair_f))

        p_array[t] = density*temp - density/(3*num_atoms) * summation

    p_value = np.average(p_array,axis=0)

    return p_value, p_array


def pair_corr_func(pos, simulation_setting, disc_steps = 100):
    """
    Calculates the pair-correlation function.

    parameters
    ----------
    post_equi_pos : list
        List storing the positions of the particles in cartesian space for every time step after equilibrium is reached.
    disc_steps : int
        Number of discrete steps the box is divided in.
    returns
    ------
    g : np.ndarray
        The value of the pair-correlation function for every discrete step.
    """

    box_dim = simulation_setting['box_dim']
    volume = simulation_setting['volume']
    num_atoms = simulation_setting['num_atoms']
    num_timesteps = len(pos)
    n = np.zeros(disc_steps)
    g = np.zeros(disc_steps)
    dr = box_dim/disc_steps
    r_val = np.arange(dr, box_dim, dr)
    equilib_time_steps = len(pos)

    for t in range(num_timesteps):
        _, rel_dist = atomic_distances(pos[t], box_dim)

        for j in range(len(r_val)):
            r = r_val[j]
            a = r < rel_dist
            b = rel_dist < (r+dr)
            counts = np.count_nonzero(a*b)/2
            n[j] +=  counts

    n = n/equilib_time_steps

    for i in range(len(r_val)):
        g[i] = ( 2*volume/((num_atoms)*(num_atoms-1))*(n[i]/ (4 * np.pi * r_val[i]**2 * dr) ) )

    return g



def plot_energy(data_pos, data_vel, simulation_setting):
    '''
    Plots energies for some data.

    parameters
    ----------
    data_pos : list
        List storing the positions of the particles in cartesian space for every time step.
    data_vel : list
        List storing the velocities of the particles in cartesian space for every time step.
    simulation_setting : dictionary
        Contains all the system parameters
    '''

    box_dim = simulation_setting['box_dim']
    steps_before_equib = simulation_setting['steps_before_equib']
    timestep = simulation_setting['timestep']
    num_tsteps = len(data_pos)
    num_steps_toplot = steps_before_equib + 500

    kin_energy = np.zeros(num_steps_toplot)
    pot_energy = np.zeros(num_steps_toplot)
    tot_energy = np.zeros(num_steps_toplot)

    total_sim_time_plot = num_steps_toplot*timestep
    sim_time_equil = steps_before_equib*timestep

    for i in range(0, num_steps_toplot):

            vel = data_vel[i]
            rel_pos, rel_dist = atomic_distances(data_pos[i], box_dim)

            kin_energy[i] = ( kinetic_energy(vel) )
            pot_energy[i] = ( potential_energy(rel_dist) )
            tot_energy[i] = ( kin_energy[i]+ pot_energy[i] )

    tau, error_pot = auto_corr_func(pot_energy)
    tau, error_kin = auto_corr_func(kin_energy)
    tau, error_tot = auto_corr_func(tot_energy)

    x = np.linspace(0, total_sim_time_plot, num_steps_toplot)

    plt.plot(x, kin_energy, 'r-', label = 'Kinetic energy')
    plt.fill_between(x, kin_energy - error_kin, kin_energy + error_kin, color = 'grey')

    plt.plot(x, pot_energy, 'b-', label = 'Potential energy' )
    plt.fill_between(x, pot_energy - error_pot, pot_energy + error_pot, color = 'grey')

    plt.plot(x, tot_energy, 'g-', label = 'Total energy')
    plt.fill_between(x, tot_energy - error_tot, tot_energy + error_tot, color = 'grey')

    plt.axvline(x=sim_time_equil)
    plt.xlabel('Simulation time')
    plt.ylabel('Energy')
    plt.legend(loc = 'upper right')
    plt.show()
