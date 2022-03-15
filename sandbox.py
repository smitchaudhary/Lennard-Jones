import numpy as np
from skeleton import *
import matplotlib.pyplot as plt
from time import time
from observables import*


#POSITIONS
p0_pos = np.array([0.,0,0])
p1_pos = np.array([6.,0,0]) #2.12246204831
p2_pos = np.array([12.,0,0])
box_dim = 5
pos = np.array([p0_pos, p1_pos, p2_pos])


rel_pos,rel_dist = atomic_distances(pos,box_dim)

print(rel_dist)

init_pos = np.array(fcc_lattice(5,1))

#print(np.shape(init_velocity(3,1)))
#print(init_velocity(3,1))

#VELOCITIES
#p0_vel = np.array( [2,0,0])
#p1_vel = np.array( [0,0,0])
#p2_vel = np.array( [-2,0,0])

#vel = np.array([p0_vel, p1_vel, p2_vel])
#kinetic_energy(init_velocity(10,2))

rel_pos, rel_dist = atomic_distances(init_pos,1)
potential_energy(rel_dist)

lj_force(rel_pos, rel_dist)

def normal_autocorr(mu, sigma, tau, N):
    """Generates an autocorrelated sequence of Gaussian random numbers.
    
    Each of the random numbers in the sequence of length `N` is distributed
    according to a Gaussian with mean `mu` and standard deviation `sigma` (just
    as in `numpy.random.normal`, with `loc=mu` and `scale=sigma`). Subsequent
    random numbers are correlated such that the autocorrelation function
    is on average `exp(-n/tau)` where `n` is the distance between random
    numbers in the sequence.
    
    This function implements the algorithm described in
    https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
    
    Parameters
    ----------
    
    mu: float
        mean of each Gaussian random number
    sigma: float
        standard deviation of each Gaussian random number
    tau: float
        autocorrelation time
    N: int
        number of desired random numbers
    
    Returns:
    --------
    sequence: numpy array
        array of autocorrelated random numbers
    """
    f = np.exp(-1./tau)
    
    sequence = np.zeros(shape=(N,))
    
    sequence[0] = np.random.normal(0, 1)
    for i in range(1, N):
        sequence[i] = f * sequence[i-1] + np.sqrt(1 - f**2) * np.random.normal(0, 1)
    
    return mu + sigma * sequence


#test = normal_autocorr(0, 1, 50  , 20000) # normal_autocorr(mu, sigma, tau, N)

#print(auto_corr_func(test))


def week4():
    #-----------  INPUT  ----------

    #------Thermodynamic input-----
    density = 1.2
    temp = 0.5
    #------Simulation input -------
    num_tsteps = 5000
    timestep =  0.001
    d_l_ratio = 2
    #------------------------------

    num_atoms = 4*d_l_ratio**3
    box_dim = (num_atoms/density)**(1/3)
    volume = box_dim**3

    #rel_pos, rel_dist = atomic_distances(init_pos, box_dim)
    #print('Potential energy: ', potential_energy(rel_dist) )

    init_pos = np.array(fcc_lattice(d_l_ratio,box_dim))
    num_atoms = init_pos.shape[0]
    #num_atoms = 10000
    print('Number of atoms: ',num_atoms)
    init_vel = init_velocity(num_atoms, temp)

    # Checking Boltzmann distribution
    #norm_init_vel = np.linalg.norm(init_vel, axis = 1)
    #plt.hist(norm_init_vel, bins = 75)
    #plt.savefig('Plots/Boltzmann_distribution')

    # Simukeleation
    tic = time()
    data_pos, data_vel, g, r_val = simulate(init_pos, init_vel, temp, num_tsteps, timestep, box_dim)
    toc = time()
    print(f'Time taken {toc - tic} seconds')
    #print(g, r_val)
    print(np.shape(g), np.shape(r_val))
    plt.plot(r_val, g)
    plt.savefig('Pair-correlation_function')
    plt.show()

    kin_energy = []
    pot_energy = []
    tot_energy = []
    dist = []

    for i in range(num_tsteps):
        
        vel = data_vel[i]
        rel_pos, rel_dist = atomic_distances(data_pos[i], box_dim)
        dist.append(rel_dist[0][1])
        
        kin_energy.append( kinetic_energy(vel) )    
        pot_energy.append( potential_energy(rel_dist) )
        tot_energy.append( kin_energy[-1]+ pot_energy[-1] )
        
    plt.plot(kin_energy, 'r-', label = 'Kinetic energy')
    plt.plot(pot_energy, 'b-', label = 'Potential energy' )
    plt.plot(tot_energy, 'g-', label = 'Total energy')
    plt.legend()
    plt.show() 
    
    
    xpos0 = [data_pos[i][0][0] for i in range(0, num_tsteps)]
    ypos0 = [data_pos[i][0][1] for i in range(0, num_tsteps)]
    zpos0 = [data_pos[i][0][2] for i in range(0, num_tsteps)]

    #plt.figure(figsize=(15, 8))
    #plt.plot(xpos0, 'k*', label = 'X Position of p0')
    #plt.plot(ypos0, 'm*', label = 'Y Position of p0')
    #plt.plot(zpos0, 'y*', label = 'Z Position of p0')
    #plt.ylim([-argon_lattice_const, argon_lattice_const])


    #plt.show()
    #plt.savefig('Plots/3_particles_simulation')


    #plt.plot(dist, 'r-', label = 'Inter atomic distance')
    #plt.title("Inter atomic distance for colliding particles")
    #plt.xlabel("Iteration number")
    #plt.ylabel("Distance in natural units")
    #plt.legend()
    #plt.show()



def week1_test():
    #POSITIONS
    p0_pos = np.array([0.5,1])
    p1_pos = np.array([0.5,2.02246204831]) #2.12246204831

    pos = np.array([p0_pos, p1_pos])

    #VELOCITIES
    p0_vel = np.array( [0.6,0.8])
    p1_vel = np.array( [0.6,0.8])

    vel = np.array([p0_vel, p1_vel])

    #TEST OF ATOMIC DISTANCES FUNCTION
    rel_pos, rel_dist = atomic_distances(pos, 100)
    force = lj_force(rel_pos, rel_dist)
    print(rel_pos[1][0][:])
    print()
    print(rel_dist)

    #TEST OF LJ FORCE
    print()
    print(force)

    #TEST OF THE ENERGIES

    potential_energy = potential_energy(rel_dist)
    kinetic_energy = kinetic_energy(vel)
    total_energy = total_energy(rel_dist, vel)
    print()
    print(f'Potential energy is {potential_energy}')
    print()
    print(f'Kinetic energy is {kinetic_energy}')
    print()
    print(f'Total energy is {kinetic_energy}')


'''
def fcc_lattice(num_atoms, lat_const):
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system
    lattice_const : float
        The lattice constant for an fcc lattice

    Returns
    -------
    pos_vec : np.ndarray
        Array of particle coordinates
    """
    box_dimension = 1
    num_unit_cells = box_dimension/lat_const
    dim = 3
    pos_vec = np.array(num_atoms, dim)

    for l in range(num_unit_cells):
        for i in range(num_atoms):
            for k in range(dim):
                pos_vec[i,l] = (box_dimension/lat_const) * l

    return pos_vec

'''
