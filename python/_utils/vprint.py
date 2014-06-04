"""
    Print in MPI environment with given verbosity level.
    Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
"""

import sys

""" 
    MPI print function: only prints if this is master.
"""
def mpiprint(msg):
    if (IAmMaster):
        sys.stdout.write(str(msg))
"""
    Vetobosity print: prints according to verbosity level.
"""
def vprint(msg, verb):
    if (verbosity >= verb):
        mpiprint(msg)

""" 
    Print with normal messages.
"""
def nprint(msg):
    vprint(msg, MSG_NORMAL)
    
""" 
    Print with debug messages.
"""
def dprint(msg):
    vprint(msg, MSG_DEBUG)

""" 
    Print with everything.
"""
def eprint(msg):
    vprint(msg, MSG_DUMP)


# Verbosity levels
MSG_QUIET  = 0             # Do not print any message
MSG_NORMAL = 1             # Print some messages
MSG_DEBUG  = 10            # Print debug message
MSG_DUMP   = 100           # Print everyting
# Set default
verbosity    = MSG_NORMAL # Default is normal.        

# MPI master
IAmMaster = False