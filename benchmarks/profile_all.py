#!/usr/local/python3
from subprocess import check_output, call
from numpy import median, nan, zeros
from pandas import read_csv
from sys import platform
from timeit import repeat


# Input variables (review before running)

#CPU_string = 'Intel(R) Core(TM) i5-4258U CPU @ 2.40GHz'    # Macbook
CPU_string = 'Intel(R) Core(TM) i7-3770K CPU @ 3.50GHz'     # Proteome
#CPU_string = 'Intel(R) Xeon(R) E5620 @ 2.40GHz'         # Rosenberg Mac
#CPU_string = 'Intel(R) Core(TM) i5-3470 CPU @ 3.20GHz'         # quill
#CPU_string = 'Intel(R) Core(TM) i5-M540 @ 2.53GHz'      # Jettabook


#CPU_string = 'Intel(R) Core(TM) i7-3770K CPU @ 3.50GHz'     # From /proc/cpuinfo or similar
#compiler = "gcc-4.9"                                        # or 'clang', etc
#compiler = "clang"                                        # or 'clang', etc
compiler = "gcc-4.8"
#compiler = "gcc"
optimization_level = "-O2"                                  
N_trials = 5

#Default files to use (do not change)
csv_file = 'timings.csv'
profiler = 'profile.c'
compiled_program = 'test.out'

#Routine 
column_name = '{CPU_string} ({compiler} on {platform})'.format(**locals())
tests = read_csv(csv_file, index_col=0)
#tests[column_name] = nan                                    # Delete any previous profiling

for name, program in tests['Program'].items():
    if 'MATLAB' in name or "MARSAGLIA" not in program:                                    # MATLAB runtimes are profiled by hand
        continue
    if ' ' in program:
        program, extra_def = program.split()                # Extract program (for profile.c) and extra definitions
    else:
        extra_def = '-DNA'                                  # No extra compiler definitions
    if 'Fast_prng' not in name and 'Numpy' not in name:     # C/C++ profiling
        compilation_command = [compiler, optimization_level, "-D"+program, extra_def, '-o', compiled_program, profiler, '-lm']
        compilation_str = '"'+" ".join(compilation_command) 
        assert call(compilation_command, shell=platform=='win32') != 1, compilation_str + '" failed to compile.'
        startup_time   = zeros(N_trials)
        execution_time = zeros(N_trials)
        for trial in range(N_trials):
            out = check_output("./"+compiled_program).decode(encoding='UTF-8')
            out_lines = out.split('\n')
            assert len(out_lines) == 4, 'Funky output:\n' + out
            startup_time[trial]   = float(out_lines[1].split()[2])
            execution_time[trial] = float(out_lines[2].split()[5])
    else:                                                   # Python profiling
        continue
        startup_time = nan
        execution_time = repeat(stmt=program+'(size={:})'.format(1000000000), setup='import numpy' if 'numpy' in program else 'import fast_prng', number=1, repeat=N_trials)

    print('{name}: {:g} (median us, startup); {:g} (median ns, per PRN)'.format(median(startup_time), median(execution_time), **locals()))
    tests.loc[name, column_name] = median(execution_time)

tests.to_csv(csv_file)

