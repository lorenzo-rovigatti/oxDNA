from oxDNA_analysis_tools.UTILS.RyeReader import Configuration, TopInfo
from IPython.display import display, HTML
from random import randint
from json import dumps
    
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

def display_files(files_with_ext,  inbox_settings =  ["Monomer", "Origin"], oxview_src = "https://sulcgroup.github.io/oxdna-viewer/"):
    """
        generate an iframe displaying the provided files in oxview
        returns the iframe-id for reuse
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
    # now let's create all those files and extensions we want to pass to view
    for file_string, ext in files_with_ext:
        a(f"t_files.push(new Blob(['{file_string}'], {{type : 'text/plain'}}));")
        a(f't_ext.push("{ext}");')
    # get the reference to the iframe
    a(f"const frame = document.getElementById('oxview-frame-{frame_id}');")
    # and forward the files
    a(f"frame.contentWindow.postMessage({{message : 'iframe_drop',files: t_files, ext: t_ext, inbox_settings : {inbox_settings} }}, \"{oxview_src}\");")
    a("}")
    a("</script>")
    a(f'<iframe width="99%" height="500"  src="{oxview_src}" id="oxview-frame-{frame_id}" onload="handle()">')
    a("</iframe>")
    
    display(HTML( "".join(out_lines) ))
    
def from_path(*args:[str], **kwargs):
    """ 
        display oxview frame based on the string path provided 
        args - contains the paths to the files
        kwargs the properties to oxview defaults to = {"inbox_settings":["Monomer", "Origin"], "oxview_src" : "https://sulcgroup.github.io/oxdna-viewer/"}
        usage:
        raw("conf.top", "conf.dat",**{"inbox_settings":["Monomer", "Origin"]})
    """
    file_list = [(__fetch_file_from_path(path),path) for path in args]
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
    
    
    
    
    
def oxdna_conf(top: TopInfo, conf:Configuration, overlay = None, forces_path = None, par_file_path = None , inbox_settings =  ["Monomer", "Origin"], oxview_src = "https://sulcgroup.github.io/oxdna-viewer/"):
    # compress the dat
    dat_string = __compress_dat(conf)
    # load up the top file
    top_string = __fetch_file_from_path(top.path)

    # location to store the things we want to past to js
    file_list = [(top_string, "top"), (dat_string, "dat")]
   
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
        
    display_files(file_list, inbox_settings, oxview_src)
    
def loro_patchy_conf(top_path:str, conf:Configuration,  matrix_path:str, inbox_settings =  ["Monomer", "Origin"], oxview_src = "https://sulcgroup.github.io/oxdna-viewer/"):
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

    display_files(file_list, inbox_settings, oxview_src)

def flro_patchy_conf(top_path:str,  conf:Configuration, particles_path:str, inbox_settings =  ["Monomer", "Origin"], oxview_src = "https://sulcgroup.github.io/oxdna-viewer/"):
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
    display_files(file_list, inbox_settings, oxview_src)
    
