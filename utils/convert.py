import os, sys

def prefix_path(prefix, path, input_dir):
    base, filename = os.path.split(path)

    if input_dir:
        return os.path.join(base, prefix + filename)
    else:
        return prefix + filename

def print_inverted_configuration(old_filename, strand_lengths, new_filename):
    with open(old_filename) as old_conf, open(new_filename, "w") as new_conf:
        lines = list(old_conf.readlines())
        
        # split the lines array into chunks, one for each strand
        strand_lines = []
        start = 3
        for length in strand_lengths:
            strand_lines.append(lines[start : start + length])
            start += length
        
        new_conf_content = "".join(lines[0:3]) # headers
        for strand in strand_lines:
            new_conf_content += "".join(list(reversed(strand)))
        
        new_conf.write(new_conf_content)

def old_to_new(topology, configuration, prefix, input_dir):
    print("INFO: converting from old to new", file=sys.stderr)
    
    # convert the topology
    with open(topology) as f:
        N, N_strands = [int(x) for x in f.readline().split()]
        strands = {str(i) : [] for i in range(1, N_strands + 1)}
        for line in f.readlines():
            spl = line.split()
            strands[spl[0]].append(spl[1:])
            
    with open(prefix_path(prefix, topology, input_dir), "w") as new_t:
        print(N, N_strands, "5->3", file=new_t)

        strand_lengths = []        
        for k, v in strands.items():
            circular = True
            sequence = []
            for nucleotide in v:
                if nucleotide[1] == "-1" or nucleotide[2] == "-1":
                    circular = False
    
                try:
                    base = str(int(nucleotide[0])) # yields an exception if nucleotide[0] is a letter                
                    sequence.append("(" + nucleotide[0] + ")")
                except ValueError:
                    sequence.append(nucleotide[0])
                
            line = "".join(sequence[::-1])
            strand_lengths.append(len(sequence))
            if circular:
                line += " circular=True"
            else:
                line += " circular=False"
            print(line, file=new_t)
            
    print_inverted_configuration(configuration, strand_lengths, prefix_path(prefix, configuration, input_dir))

def new_to_old(topology, configuration, prefix, input_dir):
    print("INFO: converting from new to old", file=sys.stderr)
    
    # convert the topology
    with open(topology) as f:
        N, N_strands = [int(x) for x in f.readline().split()[0:2]]
        sequences = []
        is_circular = []
        for line in f.readlines():
            spl = line.split()
            circular = False
            for option in spl[1:]:
                key, value = [str.lower() for str in option.split("=")]
                if key == "circular":
                    circular = value == "true"
            
            is_circular.append(circular)
            
            open_parenthesis = False
            parenthesis_token = ""
            sequence = []
        
            for c in spl[0]:
                if c == '(':
                    open_parenthesis = True;
                    parenthesis_token = ""
                elif c == ')':
                    if not open_parenthesis:
                        print(f"ERROR: unbalanced parenthesis in sequence '{spl[0]}'", file=sys.stderr)
                    open_parenthesis = False
                    sequence.append(parenthesis_token)
                else:
                    if open_parenthesis:
                        parenthesis_token += c
                    else:
                        sequence.append(c)
                        
            sequences.append(sequence[::-1])
                        
    with open(prefix_path(prefix, topology, input_dir), "w") as old_t:
        print(N, N_strands, file=old_t)
        current_idx = 0
        for i, sequence in enumerate(sequences):
            strand = i + 1
            first = current_idx
            last = current_idx + len(sequence) - 1
            for nucl in sequence:
                prev = current_idx - 1
                next = current_idx + 1
                if current_idx == first:
                    if is_circular[i]:
                        prev = last
                    else:
                        prev = -1
                if current_idx == last:
                    if is_circular[i]:
                        next = first
                    else:
                        next = -1
                current_idx += 1
                
                print(strand, nucl, prev, next, file=old_t)

    strand_lengths = list(map(len, sequences))
    print_inverted_configuration(configuration, strand_lengths, prefix_path(prefix, configuration, input_dir))

if __name__ == '__main__': 

    import argparse
    
    parser = argparse.ArgumentParser(description="A utility to convert oxDNA topology and configuration files from the old to the new format (and viceversa).")
    parser.add_argument("topology", help="The source topology file")
    parser.add_argument("configuration", help="The source configuration file")
    parser.add_argument("-p", "--prefix", default="", required=False, help="A string that will prepended to the output filenames")
    parser.add_argument("-i", "--input-dir", default=False, action="store_true", required=False, help="Output files will be written to the directory containing the input files instead of the current working directory")
    
    args = parser.parse_args()
    if args.prefix == "":
        args.prefix = "converted_"

    with open(args.topology) as top:
        spl = top.readline().split()
        if len(spl) == 2:
            old_to_new(args.topology, args.configuration, args.prefix, args.input_dir)
        elif len(spl) == 3 and spl[2] == "5->3":
            new_to_old(args.topology, args.configuration, args.prefix, args.input_dir)
        else:
            print("ERROR: Invalid topology file", file=sys.stderr)
            exit(1)
        