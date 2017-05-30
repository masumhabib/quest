#!/usr/bin/python

#------------------------------------------------------------------------------
# Include local lib
#------------------------------------------------------------------------------
import sys, os
from os.path import dirname, realpath
REGRESSION_ROOT = dirname(dirname(realpath(__file__)))
sys.path.append (REGRESSION_ROOT+"/lib")


#------------------------------------------------------------------------------
# Imports
#------------------------------------------------------------------------------
import argparse, textwrap
import pod

# Globals
#opts = pod.POD()


#------------------------------------------------------------------------------
"""
 The main() function.
"""
#------------------------------------------------------------------------------
def main(argv = None):
    if argv is None:
        argv = sys.argv

    opts = set_default_opts()
    parse_args(argv, opts)
    perform_actions(opts)

    return 0

def set_default_opts():
    opts = pod.POD()
    return opts 

#------------------------------------------------------------------------------
# Parse commandline arguments
#------------------------------------------------------------------------------
def parse_args(argv, opts):
    parser = argparse.ArgumentParser(
            usage="%(prog)s [options] <control_file.cntl> command",
            description="This script runs regression tests.",
            epilog="",
            formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument("-d", "--debug", action="store_true", 
            help="Run the script in debug mode.")
    #parser.add_argument("-r", "--release", action="store_true", default=True,
    #        help="build QUEST with release mode")
    #parser.add_argument("--doc", action="store_true", 
    #        help="build documentation")
    #parser.add_argument("--boost", action="store_true", 
    #        help="build Boost library from source")
    #parser.add_argument("--lapack", action="store_true", 
    #        help="build LAPACK library from source")
    #parser.add_argument("--arpack", action="store_true", 
    #        help="build ARPACK library from source")
    #parser.add_argument("--superlu", action="store_true", 
    #        help="build SuperLU library from source")
    #parser.add_argument("--armadillo", action="store_true", 
    #        help="build Armadillo library from source")
    #parser.add_argument("--mpich", action="store_true", 
    #        help="build MPICH2 library from source")
    #parser.add_argument("--dest", type=str, metavar="/path/to/destination",
    #        help="destination path where the action will be performed.")
    #parser.add_argument("-n", "--num_cpus", type=int, metavar="num_cpus", default=8,
    #        help="number of CPUs that will be used for compilation")

    commands = parser.add_argument_group("commands") 
    commands.add_argument("commands", type=str, nargs='+', default="all",
            choices=["clean", "prepare", "run", "analyze", "report", 
                     "update", "all", "help"],
            
            help=textwrap.dedent(
'''
clean       Cleans the run area.
prepare     Prepares the run area.      
run         Runs tests.
analyze     Analyzes results.
report      Generates reports.
update      Updates the golden.
all         Perform all the steps.

help:       Show this help message and exit.''')
            )

    args = parser.parse_args()

    if "help" in args.commands:
        parser.print_help ()
        exit (0)

    if "all" in args.commands:
        args.commands = ["clean", "prepare", "run", "analyze", "report", "update"]

    opts.commands = args.commands
    opts.debug = args.debug

    return opts

def perform_actions(opts):
    for command in opts.commands:
        if command == "clean":
            pass
        elif command == "prepare":
            pass
        elif command == "run":
            pass
        elif command == "analyze":
            pass
        elif command == "report":
            pass
        elif command == "update":
            pass
        else:
            print ("Uknown command " + command + ", run '-h' or 'help' for usage details.")
            exit (1)


#------------------------------------------------------------------------------
"""
 Entry point.
"""
#------------------------------------------------------------------------------
if __name__ == "__main__":
    sys.exit(main())



