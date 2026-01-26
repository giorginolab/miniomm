# MiniOMM

```
            _         _                              
 _ __ ___  (_) _ __  (_)  ___   _ __ ___   _ __ ___  
| '_ ` _ \ | || '_ \ | | / _ \ | '_ ` _ \ | '_ ` _ \ 
| | | | | || || | | || || (_) || | | | | || | | | | |
|_| |_| |_||_||_| |_||_| \___/ |_| |_| |_||_| |_| |_|
                                                     
```


A simple, HPC-friendly OpenMM wrapper for common-case MD runs, with minimal dependencies.

It works if it works for me.

[![PyPI version](https://badge.fury.io/py/miniomm.svg)](https://badge.fury.io/py/miniomm)

## Rationale

Developed to launch OpenMM runs on recent GPU-endowed machines,
including those with ppc64le architecture. This is because OpenMM can
currently be installed (also via Conda), but several related packages
such as `mdtraj` cannot. MiniOMM aims to provide a "minimal working"
environment for MD runs without requiring C++ or Python coding.


## Installation

In principle, the usual `pip install -e .` should be sufficient
after checking out this repository. The version on PyPI currently
lags behind. 

However, you will need additional software:

 * OpenMM and its Python interface
 * Optionally, Plumed, openmm-plumed and its Python interface

See the [Wiki](https://github.com/giorginolab/miniomm/wiki)
for (somewhat site-specific) installation instructions.


## Features

Supports

 * NVT (constant volume) production simulations with PME electrostatics and explicit solvent
 * NPT (constant pressure) equilibration
 * Runs pre-built systems in **AMBER** (prmtop) and **CHARMM** (psf) formats
 * Checkpoints and restarts are enabled out of the box
 * Plumed, if openmm-plumed is installed
 * Tested on Linux x64 and ppc64le, with and without NVIDIA GPUs

Does not support

 * Restraints of any kind (tentative implementation is in place)

May support in the future, given sufficient interest

 * Any feature provided by OpenMM, e.g. custom potentials.


The script is designed to be idempotent, that is, you may stop and
restart it repeatedly, and it will progress until the end of the
simulation. This may be convenient for time-limited batch systems.


# Usage 

Go to a directory containing structure files and input, then just

    miniomm

Useful options are

* `--help`
* `--input` (reads an input file named anything different than `input`)
* `--platform` (to run without a GPU)
* `--device` (to select a different GPU).


## Input format

MiniOMM requires the input file to:
- be called just "input"
- key and values must be on the same line
- key and values must be space/tab separated
- trajectory is saved as .dcd con nom fisso "output"
- restart (chk) files are automatically searched for and generated
- if bolean on = true, off = false

# Syntax


For now, the syntax is only documented in the examples files. It
is a simple list of "keyword value" pairs, with sensible
defaults, and largely compatible with ACEMD.



```
#--- System Structure and Parameters, Initial Conditions ---
structure                 # Path to structure file (.pdb, .psf, .prmtop)
system                    # rPath to pre-built systems from .prmtop (AMBER) or .psf (CHARMM)
parameters                # Path to CHARMM parameters to resume simulation (.prm file)
parmfile                  # Path to AMBER  parameters to resume simulation ()
openmmsystem              # Path to openMM parameters to start/resume simulation (.xml)
openmmstate               # Path to resumes simulation with latest velocities, position and time saved at (.xml)
         # if no chk or no openmmstate in the same folder, simulation starts from the beginning 
coordinates               # Path to .pdb to read atomic coordinates  (mandatory if no bincoordinates)
bincoordinates            # Resumes simulation from NAMBD file .coor (mandatory if no coordinates)
                          # priority: .chk > .coor > .pdb
		 # COORDINATES FROM SYSTEM (PDB) ARE STILL CONSIDERED  
boxsize                   # Box area in Angstrom. xyz in single line, space-separated (mandatory if no extendedsystem)
extendedsystem            # PAth to boxsize from NAMBD file .xsc (mandatory if no boxsize)
                          # priority: .chk > .xsc > boxsize > .pdb + coordinates

minimize                  # Number of steps of the minimization 

binvelocities             # Path to velocity NAMD files (.vel)
temperature               # Temperature for initial velocity assignment (K) 
          # priority: .chk > .vel > randomized

trajectoryperiod          # Number of steps to save in dcd/xtc
timestep                  # Dt of each step in femptoseconds (10-6 ns) 
run                       # Number of steps of the simulation (if a restart file is saved after equilibration, its length will be counted as part of the production step)
                          # total runtime = run * timestep

#--- Non-bonded Interactions ---
PME                       # Use Particle Mesh Ewald for electrostatics (bool)
cutoff                    # Non-bonded cutoff distance (Angstroms)
switchdistance            # Distance at which switching function starts (Angstroms)

#--- Temperature Control ---
thermostat                # Enable Langevin thermostat (bool)
thermostattemperature     # Target temperature (K)
thermostatdamping         # Thermostat damping (1/ps)

#--- Pressure Control ---
barostat                  # Enable barostat (bool)
barostatpressure          # Target pressure (bar)


#--- PLUMED ---
plumed                    # Path to PLUMED input file               

atomrestraint             # String or restraints specifications. Requires: atomic selection string within "", axes and width (ignored), setpoints force@time (specify time units)
                          # it is possible to add multiple atomrestraints (ref to ACEMD extforces)



# priority minim.coor > equilib.coor (NPT) > NVT ####se specificati?
```


