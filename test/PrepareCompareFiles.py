#!/usr/bin/env python

from TestSuite import *

class CompareFileGenerator(object):
    def __init__(self, folder, executable, level):
        self.log_prefix = "CompareFileGenerator:"
        
        self.base_folder = os.getcwd()
        self.folder = folder
        self.executable = os.path.abspath(executable)
        self.executable_name = os.path.basename(self.executable)
        self.level = level
            
        input_name = level + SUFFIX_INPUT
        input_path = os.path.join(folder, input_name)
        if os.path.exists(input_path) and os.path.isfile(input_path):
            self.system = System(folder, level, self.executable_name) 
                
    def launch(self):
        os.chdir(self.folder)
        
        log_file = get_log_name(self.level)
        to_execute = "%s %s log_file=%s no_stdout_energy=0" % (self.executable, self.system.input_name, log_file)
        p = sp.Popen(to_execute, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        p.wait()
        
        self.system.simulation_done(p, False)
        
        os.chdir(self.base_folder)
        
    def print_compare_file(self):
        self.system.analyser.print_compare_file()
    

def main():
    def print_usage():
        print("USAGE:")
        print("\t%s path executable test_level [-d|--debug] [-h|--help] [-v|--version]" % sys.argv[0])
        exit(1)

    def print_version():
        print("oxDNA Test Suite v 0.1")
        print("Copyright (C) 2015 oxDNA")
        print("This is free software; see the source for copying conditions.  There is NO")
        print("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")
        exit(1)

    shortArgs = 'dhv'
    longArgs = ['debug', 'help', 'version']
    
    try:
        import getopt
        args, files = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-d' or k[0] == '--debug': Logger.debug_level = 0
            if k[0] == '-h' or k[0] == '--help': print_usage()
            if k[0] == '-v' or k[0] == '--version': print_version()
            
    except Exception as e:
        print(e)
        print_usage()
        
    if len(sys.argv) < 4:
        print_usage()
        exit(1)
        
    generator = CompareFileGenerator(sys.argv[1], sys.argv[2], sys.argv[3])
    generator.launch()
    generator.print_compare_file()

if __name__ == '__main__':
    main()