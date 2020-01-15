Calculate dynamical correlation functions from LAMMPS text dump files.

Usage: dcorr 1000 1 0.002 T1000K.dump dcorr.dat [options].

```
usage: dcorr [-h] [--itype ITYPE] [--qmax QMAX] [--rtol RTOL]
             [--maxframes MAXFRAMES]
             ncorr nshift dt dumpfile outputfile

Calculate time-dependent sisf sisfX4 msd alpha2 Qt and QtX4 values.

positional arguments:
  ncorr                 window width to calculate correlations
  nshift                shift time origin for this many frames
  dt                    timestep size of MD run
  dumpfile              name of LAMMPS dump files
  outputfile            name of output file

optional arguments:
  -h, --help            show this help message and exit
  --itype ITYPE         type of atoms to include in calculation, 0 for all
                        types
  --qmax QMAX           |q| value
  --rtol RTOL           distance cutoff to check atom overlap
  --maxframes MAXFRAMES
                        number of frames to read from dump, 0 for all frames
```

