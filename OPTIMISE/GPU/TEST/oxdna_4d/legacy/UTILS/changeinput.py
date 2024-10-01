#!/usr/bin/env python

import getopt
import sys
import shutil
import os

input_options = ['T','sym_type','dt', 'steps', 'conf_file', 'topology', 'delta_translation', 'seed', 'lambda_or', 'lambda_tr', 'maxclust', 'restart_step_counter', 'init_hist_file', 'verlet_skin', 'use_average_seq']

input_options.append ('F0')
input_options.append ('stiff')
input_options.append ('diff_coeff')

# remove duplicates http://docs.python.org/faq/programming.html#how-do-you-remove-duplicates-from-a-list
input_options = list (set (input_options))

longArgs=[e + "=" for e in input_options]
shortArgs=[]

try:
    args, files = getopt.getopt(sys.argv[1:], shortArgs, longArgs)
except:
    print >> sys.stderr, "Wrong usage. Aborting"
    sys.exit (-2)

if len(args) < 1:
    print >> sys.stderr, "No point in changing nothing..."
    sys.exit (-1)

'''
print "ARGUMENTS"
for k in args:
    print k
    
print "FILES:"
for f in files:
    print f
'''

tochange = [[e[2:], v] for e, v in args]
print "## (changeinput) tochange: ", tochange
for file in files:
    inp = open (file, 'r')
    #tmpfile = './.' + file + '_new'
    tmpfile = 'input_new'
    out = open (tmpfile, 'w')
    changed = []

    print >> sys.stderr, "## (changeinput) modifying `%s\'" % file
    inlines = inp.readlines ()
    for inline in inlines:
        written = False
        for e, v in tochange:
            if inline.strip().startswith (e):
                out.write("%s = %s\n" % (e, v))
                #print >> sys.stderr, inline, "changed to", "%s = %s\n" % (e, v)
                written = True
                changed.append ([e,v])
        
        if not written:
            out.write(inline)
    
    for e, v in tochange:
        if [e, v] not in changed:
            #print >> sys.stderr, "key %s not found in original file, appending %s = %s" % (e, e, v)
            out.write("%s = %s\n" % (e, v))
    
    inp.close ()
    out.close ()
    shutil.copyfile(tmpfile, file)
    os.remove (tmpfile)

