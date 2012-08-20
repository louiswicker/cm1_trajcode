#!/usr/bin/env python

import sys
import os
import glob
import string
from optparse import OptionParser
subfolders = ['./']
# append path to the folder
sys.path.extend(subfolders)

preprocess   = "-DF2PY_REPORT_ON_ARRAY_COPY"
preprocess   = ""
fopts        = "--f90flags='-openmp'" 
fopts        = "--f90flags='-O3'"
libs         = "-liomp5"
libs         = ""
fortran_only = []
f2py_only    = ["trajsubs.f90"]


parser = OptionParser()
parser.add_option("-f","--fc",dest="compiler",type="string", help = "fortran compiler to be used valid compiler: [gfortran, ifort_x86_64]")

(options, args) = parser.parse_args()

objects = ""

# First compiler the fortran only files
 
for item in fortran_only:
    print "\n=====================================================\n"
    print "  Compiling file: %s " % item
    print "\n====================================================="
    if options.compiler == 'gfortran':
        ret = os.system('gfortran -O2 -c %s ' % (item))
    else:
        ret = os.system('ifort -O3 -c %s ' % (item))
    if ret == 0:
        objects = objects + " %s.o" % item.split(".")[0]

print "\nFortran object files compiled:  %s \n" % objects    
print "\n=====================================================\n"

for item in f2py_only:
    if options.compiler == 'gfortran':
        ret = os.system('f2py --fcompiler="gfortran" %s -c -m %s %s %s ' % (preprocess, item.split(".")[0], item, objects))
    else:
        ret = os.system('f2py --fcompiler="intelem" %s %s -c -m %s %s %s %s' % (fopts, preprocess, item.split(".")[0], item, 
                        objects,libs))
    if ret == 0:
        print "\n=====================================================\n"
        print "   Successfully compiled file: %s " % item
        print "   Object file is: %s " % (item.split(".")[0] + ".so")
        print "\n======================================================"
