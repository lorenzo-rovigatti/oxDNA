from oxDNA_analysis_tools.external_force_utils import forces

def read_force_file(file):
    """
    Read a force file into a list of force dictionaries.

    Parameters:
        file (str): file to read from

    Returns:
        force_list (list): a list of force dictionaries
    """
    force_list = []
    with open(file, 'r') as f:
        l = f.readline()
        while l:
            if "{" in l:  #new force
                l = f.readline()
                args = {}
                while "}" not in l: #read until end of description
                    l = l.split("=")
                    if l[0].strip() == "type":
                        t = l[1].strip()
                    else:
                        value = l[1].strip()
                        if len(value.split(' ')) != 1:
                            value = [float(v) for v in value.split(' ')]
                        else:
                            value = float(value)
                        args[l[0].strip()] = value
                    l = f.readline()
                force_list.append(getattr(forces, t)(**args)) #calls the function "t" from module "forces"
            l = f.readline()
    print("read {} forces".format(len(force_list)))
    return(force_list)


def write_force_file(force_list, filename, mode='w'):
    """
    Write a list of forces out to a file.

    Parameters:
        force_list (list): a list of force dictionaries
        filename (str): file to write out to
        <optional> mode (str): the mode of the python open funciton.  Defaults to 'w'
            change mode to 'w+' if you want to append instead of overwrite.
    """
    with open(filename, mode=mode) as f:
        out = ""
        for force in force_list:
            out += "{\n"
            for k in force.keys():
                out += "{} = ".format(k)
                if isinstance(force[k], list):
                    out += ", ".join(force[k])
                else:
                    out += str(force[k])
                out += "\n"
            out += "}\n\n"
        f.write(out)