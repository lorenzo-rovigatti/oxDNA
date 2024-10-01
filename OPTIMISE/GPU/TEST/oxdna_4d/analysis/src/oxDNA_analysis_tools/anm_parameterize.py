import numpy as np
import os
import json
import argparse
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.pca import align_positions
from oxDNA_analysis_tools.UTILS.RyeReader import describe, linear_read, inbox, get_confs

# calculate the position of the super particle
def mean_pos(conf,idxs):
    return np.sum(conf.positions[idxs], axis=0) / idxs.shape[0]

# calculate the particle positions
def get_superparticle_positions(conf, particles_array):
    return np.array(
        list(map(lambda idx: mean_pos(conf,idx),particles_array))
    )

# calculate deviation from mean position
def calculate_deviations(positions, reference_configuration):
    d = np.subtract(positions, reference_configuration) # Positions subtracted (dx, dy, dz)
    devs = np.sqrt(np.sum(np.square(d), axis=1)) # sqrt(dx**2 + dy**2 + dz**2)
    return devs

def anm_parameterize(particles_array:np.ndarray, trajectory:str, ref_conf:Configuration) -> np.ndarray:
    """
        Computes the coarse-grained RMSF for a given trajectory.

        Parameters:
            particles_array (np.array): An array containing the indices of the particles in each super particle.
            trajectory (str): The path to the trajectory to evaluate.
            ref_conf (Configuration): The reference configuration.

        Returns:
            np.array : The RMSF for each super particle.
    """
    
    ref_conf = inbox(ref_conf, center=True)
    ref_particles = get_superparticle_positions(ref_conf, particles_array)
        # Get trajectory information
    top_info, traj_info = describe(None, trajectory)

    # to collect the distance data of the superparticles 
    trajectory_devs = []
    for chunk in linear_read(traj_info, top_info):
        for conf in chunk:
            conf = inbox(conf, center=True)
            cur_conf_particles = get_superparticle_positions(conf, particles_array)

            # align the superparticles
            cur_conf_particles = align_positions(ref_particles, cur_conf_particles)

            current_devs = calculate_deviations(cur_conf_particles, ref_particles)
            trajectory_devs.append(current_devs)

    # make it into numpy as that is easier to compute with
    trajectory_devs = np.array(trajectory_devs)

    devs_sqrd = np.square(trajectory_devs) # Square all devs

    mean_devs_sqrd = np.mean(devs_sqrd, axis=0) # Mean of Squared Devs

    #deviations are in sim units, convert to nm in next step

    # Convert from sim units to nm
    rmsf = np.multiply(np.sqrt(mean_devs_sqrd),0.8518)

    return rmsf

def cli_parser(prog="anm_parameterize.py"):
    parser = argparse.ArgumentParser(prog = prog, description="compute par file for DNA-ANM model")
    parser.add_argument('index_file', type=str, help="Index file describing what Bases belong to which Super particle.")
    parser.add_argument('mean_file', type=str,  help="Reference configuration")
    parser.add_argument('trajectory', type=str,  help="Trajectory to evaluate")
    parser.add_argument('out_file', type=str, help="Output par file name")
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    particles_array = []
    # import index file
    with open(args.index_file) as file:
        index_lines = file.readlines()
    # parse the file
    for line in index_lines:
        particles_array.append(
                np.array(line.split(" "),dtype=int)
        )
    # parsed index file
    particles_array = np.array(particles_array,dtype=object)

    # Get the reference configuration
    top_info, ref_info = describe(None, args.mean_file)
    ref_conf = get_confs(top_info, ref_info, 0, 1)[0]

    rmsf = anm_parameterize(particles_array, args.trajectory, ref_conf)

    #Format for json output
    dict = {"RMSF (nm)": rmsf.tolist()}

    with open(args.out_file, 'w') as f:
        json.dump(dict, f)

if __name__ == '__main__':
    main()