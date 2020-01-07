Calculate dynamical correlation functions from LAMMPS text dump files.

Usage: dcorr 5000 1000 0.002 T1000K.dump dcorr.dat [options].

```
usage: dcorr [-h] [--itype ITYPE] [--qmax QMAX] [--nq NQ] [--rtol RTOL]
             [--maxframes MAXFRAMES]
             nt ndt dt dumpfile outputfile

Calculate time-dependent sisf sisfX4 msd alpha2 Qt and QtX4 values.

positional arguments:
  nt                    window width to calculate correlations
  ndt                   shift time origin for this many frames
  dt                    timestep size of MD run
  dumpfile              name of LAMMPS dump files
  outputfile            name of output file

optional arguments:
  -h, --help            show this help message and exit
  --itype ITYPE         type of atoms to include in calculation, 0 for all
                        types
  --qmax QMAX           |q| value
  --nq NQ               LL grid size: [6, 14, 26, 38, 50, 74]
  --rtol RTOL           distance cutoff to check atom overlap
  --maxframes MAXFRAMES
                        number of frames to read from dump, 0 for all frames
```

