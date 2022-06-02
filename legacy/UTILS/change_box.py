#!/usr/bin/env python

import sys
import math
import base
import shutil
import subprocess

MAX_SEED = 2**32 - 1
MAX_DELTA_E = 0.01
MAX_TRIES = 20

def change_box(new_box, inp, out):
    inp_file = open(inp)
    out_file = open(out, "w")
    inp_file.readline()
    old_box = float(inp_file.readline().split()[2])
    ratio = new_box / old_box
    print >> out_file, "t = 0"
    print >> out_file, "b = %f %f %f" % (new_box, new_box, new_box)
    print >> out_file, inp_file.readline()[:-1]
    i = 0
    for line in inp_file.readlines():
        coords = [ratio*float(x) for x in line.split()]
        out_file.write(" ".join(str(x) for x in coords) + "\n")

    out_file.close()
    inp_file.close()

def print_new_input(options, inp, out):
    inp_file = open(inp)
    out_file = open(out, "w")
    
    for line in inp_file.readlines():
        sp_line = [x.strip() for x in line.split("=")]
        if len(sp_line) == 1 or sp_line[0] not in options.keys(): out_file.write(line)

    out_file.write("# AUTOMATICALLY CHANGED/ADDED BY CHANGE_BOX.PY\n")
    for k, v in options.iteritems():
        out_file.write(str(k) + " = " + str(v) + "\n")

    out_file.close()
    inp_file.close()

if len(sys.argv) < 4:
    base.Logger.log("Usage is %s configuration input target_box [delta=0.01]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(sys.argv) > 4:
    delta = math.fabs(float(sys.argv[4]))
else: delta = 0.01

conf = sys.argv[1]
inp = sys.argv[2]
target_box = float(sys.argv[3])

conf_file = open(conf)
conf_file.readline()
old_box = float(conf_file.readline().split()[2])
conf_file.readline()
coord_lines = conf_file.readlines()
conf_file.close()

if old_box > target_box: delta = -delta
steps = int(math.ceil((target_box - old_box) / delta))
delta = (target_box - old_box) / steps
delta_E = MAX_DELTA_E / steps

options = {
    "lastconf_file" : "last_new_box.dat",
    "conf_file" : "init_new_box.dat",
    "sim_type" : "MC",
    "ensemble" : "nvt",
    "delta_translation" : "0.1",
    "delta_rotation" : "0.1",
    "print_energy_every" : "100",
    "backend" : "CPU",
    "backend_precision" : "double",
    "steps" : "1",
    "no_stdout_energy" : "0",
    "restart_step_counter" : "1",
    "verlet_skin" : "0.3"}

print_new_input(options, inp, "input_new_box")
shutil.copy(conf, "init_new_box.dat")
shutil.copy(conf, "last_new_box.dat")

args = ["/home/lorenzo/programmi/oxDNA_new/build/bin/oxDNA", "input_new_box"]
args = ["/users/randisif/tep/CURRENToxdna-code/oxDNA/build/bin/oxDNA", "input_new_box"]
#~/tep/CURRENToxdna-code/oxDNA/build/bin/oxDNA
p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
output = p.stdout.readline()
initial_E = float(output.split()[1])

new_box = old_box

options["steps"] = "10"
print_new_input(options, inp, "input_new_box")

for i in range(steps):
    new_box += delta
    print "Trying box of size %f" % new_box

    change_box(new_box, "last_new_box.dat", "init_new_box.dat")
    done = False
    tries = 0
    while not done:
        tries += 1
        options['seed'] = str(base.np.random.randint(MAX_SEED))
        print_new_input(options, inp, "input_new_box")
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = p.stdout.readlines()[-1]
        last_E = float(output.split()[1])
        print "New vs old: %f %f" % (last_E, initial_E)
        if last_E < initial_E or math.fabs((last_E - initial_E) / initial_E) < delta_E or tries > MAX_TRIES: done = True
        else: shutil.copy("last_new_box.dat", "init_new_box.dat")

    initial_E = last_E
