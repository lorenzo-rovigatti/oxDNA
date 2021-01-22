# reintegrate trajectory
import sys
from json import loads
try:
    sys.path.append("../")
    from oxdna_analysis_tools.UTILS.readers import ErikReader
except :
    print("Make sure oxdna_analysis_tools is on your python path.")
    sys.exit()

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
    
    readers= [ErikReader(f"./trajectory_{i}.dat") for i in range(len(pt_temp_list))]
    out_files = [f"./{out_file}_{T}.dat" for T in pt_temp_list]

    for i,locations in enumerate(history):
        print(i)
        for j,l in enumerate(locations):
            configuration = readers[j].read()
            configuration.write_append(out_files[l])