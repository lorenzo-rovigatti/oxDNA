from oxDNA_analysis_tools.external_force_utils.force_reader import read_force_file, write_force_file
from sys import argv

force_file = argv[1]
index_file = argv[2]

with open(index_file, 'r') as f:
            indexes = f.readline().split()
            try:
                indexes = [int(i) for i in indexes]
            except:
                raise RuntimeError("The index file must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button")

force_list = read_force_file(force_file)

for i, f in enumerate(force_list):
    if f["particle"] in indexes:
        del force_list[i]

write_force_file(force_list, "out.txt" )