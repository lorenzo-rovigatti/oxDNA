# reintegrate trajectory
from json import loads
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, get_input_parameter, write_conf

import argparse

if __name__ == "__main__":
    #handle commandline arguments
    parser = argparse.ArgumentParser(description="Retrive temperature trajectories from an REMD run.")
    parser.add_argument('input_file', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('history_file', type=str, nargs=1, help="The json history file keeping the remd run information.")
    parser.add_argument('output_file', type=str, nargs=1, help="name of the file to save the output trajectories to")
    
    args = parser.parse_args()
    input_file = args.input_file[0]
    history_file = args.history_file[0]
    out_file = args.output_file[0]


    # parse the input file to retrieve the temperature settings
    pt_temp_list = None
    with open(input_file) as file:
        for line in file:
            if "pt_temp_list" in line and not line.startswith('#'):
                pt_temp_list = line.replace(" ","").strip().split("=")
                pt_temp_list = pt_temp_list[1].split(",")
    
    with open(history_file) as file:
        history = loads(file.read())

    top_name = get_input_parameter(input_file, 'topology')
    traj_name = "trajectory"#get_input_parameter(input_file, 'trajectory_file')
    #traj_base = traj_name.split('.')[:-1]
    traj_ext = "dat"#traj_name.split('.')[-1]
    tis = []
    dis = []
    for i in range(len(pt_temp_list)):
        ti, di = describe(top_name, f'{traj_name}_{i}.{traj_ext}')
        tis.append(ti)
        dis.append(di)
        
    out_files = [f"./{out_file}_{T}.dat" for T in pt_temp_list]

    read_poses = [0 for _ in range(len(pt_temp_list))]
    for i,locations in enumerate(history):
        for j,l in enumerate(locations):
            configuration = get_confs(tis[l], dis[l], read_poses[l], 1)[0]
            read_poses[l] += 1
            write_conf(out_files[l], configuration, append=True, include_vel=dis[0].incl_v)
