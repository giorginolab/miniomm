# MiniOMM

```
            _         _                              
 _ __ ___  (_) _ __  (_)  ___   _ __ ___   _ __ ___  
| '_ ` _ \ | || '_ \ | | / _ \ | '_ ` _ \ | '_ ` _ \ 
| | | | | || || | | || || (_) || | | | | || | | | | |
|_| |_| |_||_||_| |_||_| \___/ |_| |_| |_||_| |_| |_|
                                                     
```


A simple, supercomputer-friendly OpenMM wrapper for common-case MD runs, with minimal dependencies.

It works if it works for me.

[![PyPI version](https://badge.fury.io/py/miniomm.svg)](https://badge.fury.io/py/miniomm)

## Rationale

Developed to launch OpenMM runs on recent GPU-endowed machines,
including those with ppc64le architecture. This is because OpenMM can
currently be installed (also via Conda), but several related packages
such as `mdtraj` cannot. MiniOMM aims to provide a "minimal working"
environment for MD runs without requiring C++ or Python coding.


## Installation

In principle, the usual `pip install miniomm` should be sufficient.
However, you will need additional software:

 * OpenMM and its Python interface
 * Optionally, Plumed, openmm-plumed and its Python interface

See the [Wiki](https://github.com/giorginolab/miniomm/wiki)
for (somewhat site-specific) installation instructions.


## Syntax

For now, the syntax is only documented in the examples files. It
is a simple list of "keyword value" pairs, with sensible
defaults, and largely compatible with ACEMD.


## Features

Supports

 * NVT (constant volume) production simulations with PME electrostatics and explicit solvent
 * NPT (constant pressure) equilibration
 * Runs pre-built systems in **AMBER** (prmtop) and **CHARMM** (psf) formats
 * Checkpoints and restarts are enabled out of the box
 * Plumed, if openmm-plumed is installed
 * Tested on Linux x64 and ppc64le, with and without NVIDIA GPUs

Does not support

 * Restraints of any kind

May support in the future, given sufficient interest

 * Any feature provided by OpenMM, e.g. custom potentials.


The script is designed to be idempotent, that is, you may stop and
restart it repeatedly, and it will progress until the end of the
simulation. This may be convenient for time-limited batch systems.



