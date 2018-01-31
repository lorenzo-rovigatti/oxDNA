#!/usr/bin/env python

'''
This utility requires python-plex, which can be downloaded from http://pypi.python.org/pypi/plex/
''' 

import sys
import plex as px
from glob import iglob
import os
from textwrap import wrap

try:
    from collections import OrderedDict
    ORDERED_DICT = True
except Exception:
    ORDERED_DICT = False
    print >> sys.stderr, "WARNING: OrderedDict not available, the options will be output unordered."

'''
Here we assign source files to categories in order to tidy up the output
'''
if ORDERED_DICT: CATEGORIES = OrderedDict()
else: CATEGORIES = {}

CATEGORIES["Core"] = []
CATEGORIES["MD"] = []
CATEGORIES["MC"] = []
CATEGORIES["Interactions"] = []
CATEGORIES["CUDA"] = []
CATEGORIES["Analysis"] = []
CATEGORIES["Observables"] = []
CATEGORIES["External Forces"] = []
CATEGORIES["Forward Flux Sampling (FFS)"] = []

CATEGORIES["Core"] = [
                      "Backends/SimBackend.h",
                      "Managers/SimManager.h",
                      "Lists/*.h"
                      ]

CATEGORIES["MD"] = [
                    "Backends/MDBackend.h",
                    "Backends/MD_CPUBackend.h",
                    "Backends/FFS_MD_CPUBackend.h",
                    "Backends/Thermostats/*.h"
                    ]

CATEGORIES["MC"] = [
                    "Backends/MCBackend.h",
                    "Backends/MC_CPUBackend.h",
                    "Backends/MC_CPUBackend2.h",
                    "Backends/VMMC_CPUBackend.h",
                    "Backends/MCMoves/*.h"
                    ]

CATEGORIES["Interactions"] = [
                              "Interactions/*.h"
                              ]

CATEGORIES["CUDA"] = [
                      "CUDA/CUDAForces.h",
                      "CUDA/CUDAForces.h",
                      "CUDA/*/*.h"
                      ]

CATEGORIES["Analysis"] = [
                          "Managers/AnalysisManager.h",
                          "Backends/AnalysisBackend.h"
                          ]

CATEGORIES["Observables"] = [
                             "Observables/*.h",
                             "Observables/*/*.h"
                             ]
                             
CATEGORIES["External Forces"] = [
                             "Forces/*.h"
                             ]

CATEGORIES["Forward Flux Sampling (FFS)"] = [
                     "Backends/FFS_MD_CPUBackend.h",
                     "CUDA/Backends/FFS_MD_CUDAMixedBackend.h"
]

# Should the options belonging to a certain category be joined together or left separated?
JOIN = {
        "Core" : True,
        "MD" : True,
        "MC" : True,
        "Interactions" : False,
        "CUDA" : True,
        "Analysis" : True,
        "Observables" : False,
        "External Forces" : False,
        "Forward Flux Sampling (FFS)" : True
        }

LR_TAB = '    '

def warning(position, text):
    filename = position[0].rpartition("src")[2]
    print >> sys.stderr, "WARNING: '%s', line %d: %s" % (filename, position[1], text)
    

def indent(my_str):
    return '\n'.join((LR_TAB + x) for x in my_str.split('\n'))

class Option(object):
    def __init__(self, key, value, description, optional):
        object.__init__(self)
        self.key = key.strip()
        self.value = value.strip().replace('\\n', '\n')
        self.description = description.strip().replace('\\n', '\n')
        self.optional = optional
        
    def text(self, wiki=False):
        if not wiki:
            spl = wrap(self.description, 70)
            description = "\n".join(LR_TAB + s for s in spl)
            if self.optional: return "[%s = %s]\n%s" % (self.key, self.value, description)
            else: return "%s = %s\n%s" % (self.key, self.value, description)
        else:
            key_value = "%s = %s" % (self.key, self.value)
            key_value = key_value.replace("\n", "<br />")
            if self.optional: return ";[%s]\n: %s" % (key_value, self.description)
            else: return ";%s\n: %s" % (key_value, self.description)
        
        
class OptionScanner(px.Scanner):
    def set_as_optional(self, text):
        self.current_is_optional = True
        return None
    
    def option_found(self, text):
        spl = [x.strip() for x in text.strip().split("\n")]
        for line in spl:
            if len(line) == 0: continue
            position = list(self.position())
            position[1] += spl.index(line)
            
            optional = (line[0] == "[")
            if optional and line[-1] != "]":
                warning(position, "unmatched '['")
                line += "]"
            elif not optional and line[-1] == "]":
                warning(position, "unmatched ']', treating the option as 'optional'")
                optional = True
                line = "[" + line
                
            if optional:
                type = "optional"
                line = line[1:-1]
            else: type = "required"
            
            spl_line = line.split("=", 1)
            if len(spl_line) != 2:
                warning(position, "invalid format, missing '='")
                continue
            
            key = spl_line[0].strip()
            parts = spl_line[1].partition("(")
            if parts[1] == "":
                warning(position, "invalid format, the option description should be preceeded by a '('")
                continue
            
            value = parts[0].strip()
            if parts[2][-1] != ")":
                warning(position, "invalid format, the option description should be followed by a ')'")
                continue
            
            description = parts[2][:-1]
                
            self.produce(type, [key, value, description])
    
    space = px.Any(" \t")
    lineterm = px.Str("\n") | px.Eof
    option = px.Rep(px.AnyBut("\n")) + lineterm
    begin_option_section = px.Str("@verbatim") + px.Rep(px.AnyBut("\n")) + lineterm
    end_option_section = px.Str("@endverbatim") + px.Rep(px.AnyBut("\n")) + lineterm
    options = px.Rep1(px.AnyBut("@"))
    
    lexicon = px.Lexicon([
                         (begin_option_section, px.Begin("option_section")),
                         (px.AnyChar, px.IGNORE),
                          px.State("option_section", [
                                           (end_option_section, px.Begin('')),
                                           (options, option_found),
                                           ]),
                         ])
    
    def __init__(self, filename):
        self.file = open(filename)
        px.Scanner.__init__(self, self.lexicon, self.file, filename)
        self.current_is_optional = False
        
        
class Options(object):
    def __init__(self, id):
        object.__init__(self)
        self.id = id
        if ORDERED_DICT: self.options = OrderedDict()
        else: self.options = {}

    def add_option(self, new_option):
        key = new_option.key
        if key in self.options.keys() and self.options[key].value != new_option.value:
            self.options[key].value += "/" + new_option.value
            self.options[key].description += "/" + new_option.description
        else:
            self.options[key] = new_option

    def print_options(self, wiki):
        N = self.get_N_options()
        if N == 0: return
        
        print >> sys.stderr, "Number of '%s' options: %d" % (self.id, N)
        
        if not wiki:
            print self.id + " options:\n"
            for option in self.options.itervalues():
                print indent(option.text())
            print ""
        else:
            print "===%s options===\n" % self.id
            for option in self.options.itervalues():
                print option.text(True)
            print ""

    def get_N_options(self):
        return len(self.options)

        
cwd = os.path.dirname(os.path.realpath(__file__))
base_path = os.path.join(cwd, "..", "src")

wiki = len(sys.argv) > 1 and sys.argv[1] == "wiki"
    
for cat, cat_paths in CATEGORIES.iteritems():
    join_cat = JOIN[cat]
    cat_options = Options(cat)
    for cat_path in cat_paths:
        tot_path = os.path.join(base_path, cat_path)
        for cat_file in iglob(tot_path):
            base_file = cat_file.replace(base_path + os.path.sep, "")
            file_options = Options(base_file)
            
            scanner = OptionScanner(cat_file)

            token = scanner.read()
            while token[0] != None:
                new_option = Option(token[1][0], token[1][1], token[1][2], token[0] == "optional")
                file_options.add_option(new_option)
                cat_options.add_option(new_option)
                token = scanner.read()
                
            if not join_cat: file_options.print_options(wiki)
                
    if join_cat: cat_options.print_options(wiki)
    if not wiki: print "-------------------------------------------------------------------------------\n"
