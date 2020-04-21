#!/opt/apps/intel18/python3/3.7.0/bin/python3.7
"""

        Author: Munier Salem - early-med 2010's
        Modified: Andrew Emerick

	FRAGGLE

	Designed to repeatedly run enzo from within a
	single batch submission. Originally created to
        avoid memory fragmentation issues, this restarts
        Enzo up to a user-defined number of times to attempt
        to maximize use of a single job submission without
        re-queing.

        Run using:

        python3 ./fraggle.py

        Currently assumes the executable is 'enzo.exe' and
        running on Stampede2.

	HOW IT WORKS:

          - The parameter file (if it exists) is found

          - The most recent data dump (determined by sim. time)
            is found, for both regular (e.g. DD) and redshift
            data dumps. By defualt, interpolated outputs are
            ignored.

          - Enzo is called from the most recent of the above,
            or the initial parameter file if no dumps exist.

          - If Enzo exits (for any reason) this process repeats
            for up to the number of desired times OR until
            the same output is restarted from after a specified
            number of times. This can be controlled with the
            N_TOTAL_RESTART and N_SAME_RESTART parameters

          - Code exists (but is not yet used) to modify the
            output cadence and Courant timestep to try and increase
            dump frequency or slow down the Hydro solver in
            case there are issues here. Will be added in soon (Apr 2020)
"""
# troubleshooting
#import sys
#print sys.path

import os, sys, glob
import numpy as np

from subprocess import call,check_output,CalledProcessError

MIN_TIME_STEP = 0.1 # does nothing now
MAX_TIME_STEP = 5.0 # does nothing now

#
# Dictionary of system-dependent run commands
#

RUN        = {'Stampede2' : 'ibrun mem_affinity '}

# default exectable filename
EXEC_FNAME = "./enzo.exe"

SYSTEM     = 'Stampede2' # what are we on?

VERBOSE    = True   # use -d command


def get_parameter_file():
    """
    Find the .enzo parameter file
    """

    par_name = glob.glob("*.enzo")

    if np.size(par_name) > 1:
        print("WARNING: More than one parameter file found")
        print(par_name)
        par_name = par_name[0]
        print("Using: ", par_name)

    elif np.size(par_name) == 0:
        print("WARNING: No parameter file found")
        par_name = None

    else:
        par_name = par_name[0]
        print("Found parameter file " + par_name)

    return par_name
    

def resetdtDataDump(fName,stretch):

    # open and read parameter file
    file = open(fName,'rb')
    params = file.read().split('\n')
    file.close()
    # alter run time
    key = 'dtDataDump'
    idx = None
    for i,p in enumerate(params):
        line = p.split('=')
        if line[0].strip() == key:
            val = min(float(line[1])*stretch,MAX_TIME_STEP)
            if val < MIN_TIME_STEP:
                raise RuntimeError('Run has stalled!!!')
                params[i] = " = ".join([key,str(val)])
                break
    # Print updated params ...
    file = open('tmp.par','wb')
    file.write("\n".join(params))
    file.close()

    return


def get_all_restart_dumps(DataDumpName = "DD", DataDumpDir = "DD", RedshiftDumpName = "RD",
                          RedshiftDumpDir = "RD", parameter_file = None,
                          ignore_interpolated = True):

    """
    Gathers a time-ordered (most recent sim-time first) list of available restart
    data dumps.
    """
    if not (parameter_file is None):
        if (os.path.isfile(parameter_file)):
            # 
            with open(parameter_file) as f:
                names_list = ['DataDumpDir','DataDumpName','RedshiftDumpDir',
                              'RedshiftDumpName']
                for line in f:
                    if "=" in line:
                        split_line = line.split("=")
                        if split_line[0] == 'DataDumpDir':
                            DataDumpDir = split_line[-1]
                        elif split_line[0] == 'DataDumpName':
                            DataDumpName = split_line[-1]
                        elif split_line[0] == 'RedshiftDumpName':
                            RedshiftDumpName = split_line[-1]
                        elif split_line[0] == 'RedshiftDumpDir':
                            RedshiftDumpDir  = split_line[-1]

        else:
            print(parameter_file, " path cannot be found")
            raise RuntimeError

    # look for all files

    DD_pattern = DataDumpDir + "????/" + DataDumpName + "????"
    RD_pattern = RedshiftDumpDir + "????/" + RedshiftDumpName + "????"

    DD_files = np.sort(glob.glob(DD_pattern))
    RD_files = np.sort(glob.glob(RD_pattern))

    all_files = np.array(list(DD_files) + list(RD_files))

    if np.size(all_files) == 0:
        print("No files found matching the pattern " + DD_pattern + " or " + RD_pattern)
        return None

    all_times = np.zeros(np.size(all_files))
    interpolated = np.zeros(np.size(all_files)) * False 

    # check sim time for all files and whether or not they
    # are interpolated dumps
    for fnum, fname in enumerate(all_files):
        with open(fname) as f:
            for line in f:
                if "WARNING! Interpolated output" in line: # do not restart from interpolated dumps
                    interpolated[fnum] = True
                if "InitialTime" in line:
                    all_times[fnum] = float(line.split("=")[-1])

    # reverse sort by time (most recent first)
    sort      = np.argsort(all_times)
    all_files = all_files[ sort ][::-1]
    all_times = all_times[ sort ][::-1]
    interpolated = interpolated[ sort ][::-1]

    # ignore interpolated
    if ignore_interpolated:
        return_files = all_files[ np.logical_not(interpolated) ]
    else:
        return_files = all_files


    return return_files


if __name__ == "__main__":

    if VERBOSE:
        verbose_str = ' -d'
    else:
        verbose_str = ' -'


    if np.size(sys.argv) > 1:

        if len(sys.argv) >= 2:
            SYSTEM = str(sys.arv[1])
        if len(sys.argv) >= 3:
            EXEC_FNAME = str(sys.argv[2])

    if not (os.path.isfile(EXEC_FNAME)):
        print(EXEC_FNAME + " not found. Are we in the right place?")
        raise RuntimeError

    if not SYSTEM in RUN.keys():
        print(SYSTEM + " not specified in the RUN dictionary."
              " fraggle does not know how to execute the job. Please include the command")
        raise RuntimeError

    print("Attempting to run " + EXEC_FNAME + " on system " + SYSTEM)

    #
    # find parameter file name
    #
    parfile = get_parameter_file()

    #
    # figure out the most recent data dump (DD????/DD????)
    #
    all_dumps = get_all_restart_dumps(parameter_file = parfile)

    if (all_dumps is None) and (parfile is None):
        print("Cannot find anything to start from")
        raise RuntimeError

    #
    # If there are no data dumps, start from the beginning
    #
    if all_dumps is None:
        start_command = RUN[SYSTEM] + EXEC_FNAME + verbose_str + ' ' + parfile
        call(start_command, shell=True)

    #
    # Otherwise, keep running
    # change count_max to try and do multiple restarts in a single run
    #    this can be helped by forcing Enzo to exit on output
    count_max = 1
    # 
    nlast = None; count = 0; n = 0
    prev_restart_name = None
    restart_count = 1
    while (count < count_max):

        if count > 0:
            all_dumps = get_all_restart_dumps(parameter_file = parfile)

        restart_name = all_dumps[n]

        if restart_name == prev_restart_name:
            restart_count = restart_count + 1
            print("FRAGGLE-WARN: Restarting from the same data dump again (try number %i): "%(restart_count) + restart_name)
        else:
            print("Restarting simulation from: " + restart_name)
            restart_count = 1

        restart_command = RUN[SYSTEM] + EXEC_FNAME + verbose_str + 'r ' + restart_name
        print("Using command: " + restart_command)
        call(restart_command, shell=True)

        prev_restart_name = restart_name
        count = count + 1

