from oxDNA_analysis_tools.UTILS.RyeReader import Configuration, TopInfo
from IPython.display import display, HTML
from random import randint
from json import dumps
from typing import List, Union,Dict
    
def __compress_dat(conf):
    """ generate me a compressed double escaped dat file """
    dat_string = f"t = {conf.time}\\nb = {conf.box[0]} {conf.box[1]} {conf.box[2]}\\nE = {conf.energy[0]} {conf.energy[1]} {conf.energy[2]}\\n"
    dat_string += "\\n".join([f"{p[0]:.3f} {p[1]:.3f} {p[2]:.3f} {a1[0]:.8f} {a1[1]:.8f} {a1[2]:.8f} {a3[0]:.8f} {a3[1]:.8f} {a3[2]:.8f}"
                                                                                    for p,a1,a3 in zip(conf.positions, conf.a1s, conf.a3s)])  
    # we are also throwing out the velocity data 
    return dat_string

def __fetch_file_from_path(path):
    """ fetch the file being double escaped """
    with open(path) as file:
        top_string = file.read()
        return top_string.replace("\n", "\\n")  

def display_files(system_list,  inbox_settings: List[str] =  ["Monomer", "Origin"], oxview_src:str = "https://sulcgroup.github.io/oxdna-viewer/", height:int = 500):
    """
        Generate an iframe displaying the provided files in oxview

        Parameters:
            system_list (tuple) : a list of lists of tuples (file_string, extension) 
            inbox_settings (list[str]) : a list of strings, the inbox settings to use
            oxview_src (str) : the url of the oxview source
            height (int) : height of the iframe

        Returns: 
            The iframe-id for reuse
    """
    #generate a unique id for our iframe
    frame_id = str(randint(0,1000000))
    # buffer where our html will go to
    out_lines = []
    a = out_lines.append

    a("<script>")
    a("function handle(){")
    a("let t_files = [];")
    a("let t_ext = [];")    
    # get the reference to the iframe
    a(f"const frame = document.getElementById('oxview-frame-{frame_id}');")    
    for files_with_ext in system_list: 
        # push only 1 system per postMessage 
        a("t_files = [];")
        a("t_ext = [];")
        # now let's create all those files and extensions we want to pass to view
        for file_string, ext in files_with_ext:
            a(f"t_files.push(new Blob(['{file_string}'], {{type : 'text/plain'}}));")
            a(f't_ext.push("{ext}");')
        # and forward the files
        a(f"frame.contentWindow.postMessage({{message : 'iframe_drop',files: t_files, ext: t_ext, inbox_settings : {inbox_settings} }}, \"{oxview_src}\");")
        
    a("}")
    a("</script>")
    a(f'<iframe width="99%" height="{height}"  src="{oxview_src}" id="oxview-frame-{frame_id}" onload="handle()">')
    a("</iframe>")
    
    display(HTML( "".join(out_lines) ))
    
def from_path(*args:Union[List[str],List[List[str]]] , **kwargs):
    """ 
        Display oxview frame based on the string path provided 

        Parameters:
            args (Union[List[str],List[List[str]]]) : contains the paths to the files or a list of lists of paths
            kwargs (dict) : the properties to oxview defaults to = {"inbox_settings":["Monomer", "Origin"], "oxview_src" : "https://sulcgroup.github.io/oxdna-viewer/"}
        
        Usage:
            from_path("conf.top", "conf.dat",**{"inbox_settings":["Monomer", "Origin"]})
    """

    # we have two possible caseses 
    # args is a list of lists or a list of strings
    if all(isinstance(i, str) for i in args):
        file_list = [[(__fetch_file_from_path(path),path) for path in args]]
        # print("single")
        # print (file_list)
        # #file_list = [(__fetch_file_from_path(path),path) for path in args]
    elif all(isinstance(i, list) for i in args):
        # print("multi")
        file_list = file_list = [[(__fetch_file_from_path(path), path)  for path in lst] for lst in args]
        # print (file_list)
    else:
        raise ValueError("args must be a list of strings or a list of lists of strings")

    
    #make sure we have some default view
    if not "inbox_settings" in kwargs:
        inbox_settings = ["Monomer", "Origin"]
    else:
        inbox_settings = kwargs["inbox_settings"]
    if not "oxview_src" in kwargs:
        oxview_src = "https://sulcgroup.github.io/oxdna-viewer/"
    else:
        oxview_src = kwargs["oxview_src"]    
    display_files(file_list, inbox_settings, oxview_src)
    
    
    
    
    
def oxdna_conf(top: TopInfo, conf:Configuration, overlay:Dict[str,List] = None, forces_path:str = None, par_file_path:str = None , script_file_path:str = None, inbox_settings:List[str] = ["Monomer", "Origin"], oxview_src:str = "https://sulcgroup.github.io/oxdna-viewer/", height:int = 500):
    """
        Display an oxDNA configuration in oxview

        Parameters:
            top (TopInfo) : the top file data
            conf (Configuration) : the configuration data
            overlay Dict[str:List] : (optional) dictionary with the color overlay
            forces_path (str) : (optional) the path to the forces file
            par_file_path (str) : (optional) the path to the par file
            script_file_path (str) : (optional) the path to the script file (js)
            inbox_settings (List[str]) : (optional) a list of strings, the inbox settings to use
            oxview_src (str) : (optional) the url of the oxview source
            height (int) : (optional) the height of the view
    """
    # compress the dat
    dat_string = __compress_dat(conf)
    # load up the top file
    top_string = __fetch_file_from_path(top.path)

    # location to store the things we want to past to js
    file_list = [(top_string, "top"), (dat_string, "dat")]

    # handle script files
    if script_file_path:
        script_string = __fetch_file_from_path(script_file_path)
        file_list.append((script_string, "js"))

    # if we have an overlay it's supposedly a json-like object
    if overlay:
        overlay_string = dumps(overlay)
        file_list.append((overlay_string, "json"))
    # handle par files    
    if par_file_path:
        parfile_string = __fetch_file_from_path(par_file_path)
        file_list.append((parfile_string, "par"))
    
    #show forces
    if forces_path:
        forces_string = __fetch_file_from_path(forces_path)
        file_list.append((forces_string, "forces.txt"))      
        
    # hack cause we support multiple systems in from_path 
    display_files([file_list], inbox_settings, oxview_src, height=height)
    
def loro_patchy_conf(top_path:str, conf:Configuration,  matrix_path:str, inbox_settings:List[str] = ["Monomer", "Origin"], oxview_src:str = "https://sulcgroup.github.io/oxdna-viewer/"):
    """
        Display a loro patchy configuration in oxview
        
        Parameters:
            top (str) : the top file path
            conf (Configuration) : the configuration data
            matrix_path (str) : the path to the matrix file
            inbox_settings (list[str]) : (optional) a list of strings, the inbox settings to use
            oxview_src (str) : (optional) the url of the oxview source
    """
    # compress the dat
    dat_string = __compress_dat(conf)
    # load up the top file and the matrix
    top_string = __fetch_file_from_path(top_path)
    matrix_string = __fetch_file_from_path(matrix_path)

    # locationr to store the things we want to past to js
    file_list = [
        (top_string,   "LORO.top"), 
        (dat_string,   "LORO.dat"),
        (matrix_string,"matrix")
    ]

    display_files([file_list], inbox_settings, oxview_src)

def flro_patchy_conf(top_path:str,  conf:Configuration, particles_path:str, inbox_settings:List[str] = ["Monomer", "Origin"], oxview_src:str = "https://sulcgroup.github.io/oxdna-viewer/"):
    """
        Display a flro patchy configuration in oxview

        Parameters:
            top (str) : the top file path
            conf (Configuration) : the configuration data
            patricles (str) : the path to the particles file
            inbox_settings (list[str]) : (optional) a list of strings, the inbox settings to use
            oxview_src (str) : (optional) the url of the oxview source
    """
    # compress the dat
    dat_string = __compress_dat(conf)
    # load up the top file
    top_string = __fetch_file_from_path(top_path)
    particles_string = __fetch_file_from_path(particles_path)
    # locationr to store the things we want to past to js
    file_list = [
        (top_string,   "top"), 
        (dat_string,   "dat"),
        (particles_string,"particles")
    ]
    display_files([file_list], inbox_settings, oxview_src)
    
