import numpy as np
from numpy.linalg import norm
from itertools import product

def simulate(init_pos, init_vel, temp, num_tsteps, timestep, box_dim, euler_method = False ):
    """
    Molecular dynamics simulation using the Euler or Verlet's algorithms
    to integrate the equations of motion.

    Parameters
    ----------
    init_pos : np.ndarray
        The initial positions of the atoms in Cartesian space in units of sigma
    init_vel : np.ndarray
        The initial velocities of the atoms in Cartesian space in units of sqrt(epsilon\mass)
    temp : float
        The (unitless) temperature of the system.<
    num_tsteps : int
        The total number of simulation steps
    timestep : float
        Duration of a single simulation step
    box_dim : float
        Dimensions of the simulation box in units of sigma
    euler_method : bool (by default False)
        The simulation uses by default the Verlet algorithm. To use Euler's algorithm set to True.

    Returns
    -------
    data_pos : list
        Position of all particles over time
    data_vel : list
        Velocity of all particles over time.
    steps_before_equilib : int
        Number of steps the simulation has run before reaching equilibrium. Will be -1 if equilibrium has not been reached.
    avg_temp : float
        Temperature of the system based on the kinetic energy.
    """

    #CONVENTION: pos is defined as a (number of particles)xdim numpy array where each row represents a particle
    #and the columns store the coordinates of the particle in each dimension

    #lists to store particles positions and velocities at each time step
    print('Starting simulation')
    data_pos = []
    data_vel = []

    new_pos = init_pos
    new_vel = init_vel

    data_pos.append(new_pos)
    data_vel.append(new_vel)

    #variables for rescaling
    volume = box_dim**3
    rescale_flag = True
    rescale_flag2 = True
    #rescale_step= np.ceil(0.4/timestep) # let it run 0.4 seconds at least to rescale
    rescale_step = 20
    kin_ener = 0

    steps_before_equilib = -1
    equilib_time_steps = 0
    num_particles = init_pos.shape[0]

    for i in range(num_tsteps):
        # Pos and vel in previous step
        pos = np.copy(new_pos)
        vel = np.copy(new_vel)
        kin_ener += (2*kinetic_energy(vel)) #Is actually not Kinetic energy but twice of kinetic energy.

        if ((i+1)%rescale_step == 0 and rescale_flag2 ):
            rescale_factor = np.sqrt( (num_particles-1)*3*temp / (kin_ener/rescale_step) )
            kin_ener = 0
            vel = vel*rescale_factor
            if (abs(rescale_factor - 1) < 0.01):
                if rescale_flag:
                    rescale_flag = False
                else:
                    rescale_flag2= False
                    steps_before_equilib = i+1
                    print(f'Equilibrium Reached in {i+1} steps!!')
            else:
                rescale_flag = True
        # Calculate force
        rel_pos, rel_dist = atomic_distances(pos, box_dim)
        force = lj_force(rel_pos, rel_dist)

        if euler_method: #Euler's algorithm
            new_pos = pos + vel*timestep
            new_vel = vel + force*timestep

        else:        # Verlet Algorithm
            new_pos = pos + timestep*vel + 0.5*(timestep**2)*force

            new_rel_pos, new_rel_dist = atomic_distances(new_pos, box_dim)
            new_force = lj_force(new_rel_pos, new_rel_dist)

            new_vel = vel + (timestep/2)* ( new_force + force)

        data_pos.append(new_pos)
        data_vel.append(new_vel)

    if (rescale_flag):
        print(f'Equilibrium NOT reached!')
    else:
        avg_kin_ener = kin_ener/(2*(num_tsteps-steps_before_equilib))
        avg_temp = avg_kin_ener*2/(3*(num_particles))
    return data_pos, data_vel, steps_before_equilib, avg_temp


def atomic_distances(pos, box_dim):
    """
    Calculates relative positions and distances between particles according to minimal image convention

    parameters
    ----------
    pos : np.ndarray
        The positions of the particles in cartesian space
    box_dim : float
        The dimension of the simulation box

    returns
    -------
    rel_pos : np.ndarray
        Relative positions of particles according to minimal image convention
    rel_dist : np.ndarray
        The distance between particles according to minimal image convention
    """

    part_num = pos.shape[0]
    dimensions = pos.shape[1]

    # array[i][j] is pos[i] - pos[j]
    rel_pos = np.zeros( ( part_num, part_num, dimensions)  )
    rel_dist = np.zeros( ( part_num, part_num)  )

    rel_pos = pos

    posp =   np.concatenate( [pos[np.newaxis,:,:] ]* part_num , axis=0 )
    posm =   np.concatenate( [pos[:,np.newaxis,:] ]* part_num , axis=1 )
    rel_pos = -posp + posm

    rel_pos = (rel_pos+box_dim/2) % box_dim - box_dim/2

    rel_dist = norm(rel_pos, axis=2)

    return rel_pos, rel_dist


def lj_force(rel_pos, rel_dist):
    """
    Calculates the net forces on each atom.

    Parameters
    ----------
    rel_pos : np.ndarray
        Relative particle positions as obtained from atomic_distances
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    np.ndarray
        The net force acting on particle i due to all other particles
    """

    dimensions = rel_pos.shape[2]

    rel_dist1 = np.copy(rel_dist)
    np.fill_diagonal(rel_dist1, val=2**(1/6)) #Fill the diagonal with 2^(1/6) since the diagonal is 0. Chosen to be 2^(1/6) since Force it zero at that distance.

    pair_force = 4*( 6*(1/rel_dist1)**8 - 12*(1/rel_dist1)**14)
    pair_for3d = pair_force[:,:,np.newaxis]
    pair_for3d = np.concatenate( [pair_for3d]*dimensions,axis=2 )
    total_force = pair_for3d*rel_pos
    force = np.sum(total_force,axis=0)

    return force


def fcc_lattice(d_l_ratio, box_dim):
    """
    Initializes a system of atoms on an fcc lattice.

    Parameters
    ----------
    d_l_ratio : int
        Ratio of box dimension to lattice constant
    box_dim : float
        The dimension of the simulation box

    Returns
    -------
    pos : np.ndarray
        Array of particle coordinates
    """
    lattice_const = box_dim/ d_l_ratio
    pos = []

    pos = [[[x, y, z], [x, y + 0.5, z + 0.5], [x + 0.5, y + 0.5, z], [x + 0.5, y, z + 0.5]] for x, y, z in product(range(d_l_ratio), range(d_l_ratio), range(d_l_ratio))]

    pos = np.array(pos).reshape(-1,3)
    pos = lattice_const*pos

    return pos


def kinetic_energy(vel):
    """
    Computes the kinetic energy of an atomic system.

    Parameters
    ----------
    vel: np.ndarray
        Velocity of particle

    Returns
    -------
    total_kinetic : float
        The total kinetic energy of the system.
    """

    total_kinetic = np.sum( 0.5*np.square(vel) )
    return total_kinetic


def potential_energy(rel_dist):
    """
    Computes the potential energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances

    Returns
    -------
    total_potential : float
        The total potential energy of the system.
    """

    rel_dist1 = np.copy(rel_dist)
    np.fill_diagonal(rel_dist1, val=1) #Fill the diagonal with 1 since the diagonal is 0. Chosen to be 1 since PE it zero at that distance.
    total_potential =   1/2 * np.sum(  4*( (1/rel_dist1)**12 - (1/rel_dist1)**6 ))
    return total_potential


def total_energy(rel_dist, vel):
    """
    Computes the total energy of an atomic system.

    Parameters
    ----------
    rel_dist : np.ndarray
        Relative particle distances as obtained from atomic_distances
    vel: np.ndarray
        Velocity of particle
    Returns
    -------
    total_energy = kinetic_energy + potnetial_energy : float
        The total energy of the system.
    """

    return kinetic_energy(vel) + potential_energy(rel_dist)


def init_velocity(num_atoms, temp):
    """
    Initializes the system with Gaussian distributed velocities.

    Parameters
    ----------
    num_atoms : int
        The number of particles in the system.
    temp : float
        The (unitless) temperature of the system.

    Returns
    -------
    vel_vec : np.ndarray
        Array of particle velocities
    """

    init_vel = np.random.normal(scale = np.sqrt(temp), size = 3*num_atoms).reshape(-1,3)
    cm_vel = np.average(init_vel, axis = 0) # velocity of the center of mass
    init_vel -= cm_vel
    return init_vel
