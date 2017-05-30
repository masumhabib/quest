#!/usr/bin/env python

import os, sys, subprocess

current_dir = os.path.dirname(os.path.realpath(__file__))
work_dir = current_dir + "/work"
test_command = "../test_dirac_transport.py"
log_file = "run.log"


print ("-D-: " + work_dir)
if not os.path.isdir (work_dir):
    os.makedirs (work_dir)

os.chdir (work_dir)
subprocess.call(["python", test_command])
#subprocess.call(["python", test_command], stdout=log_file, stderr=log_file)










