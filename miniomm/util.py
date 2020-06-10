import sys

from simtk.openmm import app
from simtk.openmm.openmm import Platform
import simtk.openmm as mm
import simtk.unit as u

from miniomm.reporters import *
from miniomm.namdbin import NAMDBin
from miniomm.namdxsc import *


_cachedPdb = {}
def get_pdb(n):
    if n not in _cachedPdb:
        _cachedPdb[n] = app.PDBFile(n)
    return _cachedPdb[n]



def getBanner():
    return """
            _         _                              
 _ __ ___  (_) _ __  (_)  ___   _ __ ___   _ __ ___  
| '_ ` _ \ | || '_ \ | | / _ \ | '_ ` _ \ | '_ ` _ \ 
| | | | | || || | | || || (_) || | | | | || | | | | |
|_| |_| |_||_||_| |_||_| \___/ |_| |_| |_||_| |_| |_|
                                                     
A minimalistic OpenMM MD frontend.   
https://github.com/giorginolab/miniomm
"""





def parse_boxsize_units(txt):
    box = [float(x) for x in txt.split(" ")]
    boxa = mm.Vec3(box[0], 0., 0.) * u.angstrom
    boxb = mm.Vec3(0., box[1],  0.) * u.angstrom
    boxc = mm.Vec3(0., 0., box[2]) * u.angstrom
    return (boxa, boxb, boxc)


def get_box_size(inp):
    if 'extendedsystem' in inp:
        print("Reading box size from "+inp.extendedsystem)
        (boxa, boxb, boxc) = parse_xsc_units(inp.extendedsystem)
    elif 'boxsize' in inp:
        print("Using boxsize from input string "+inp.boxsize)
        (boxa, boxb, boxc) = parse_boxsize_units(inp.boxsize)
    else:
        print("Last resort: PDB CRYST1...")
        try:
            import warnings
            print(f"Reading positions from PDB: ")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
            pdb = get_pdb(inp.coordinates)
            (boxa, boxb, boxc) = pdb.topology.getPeriodicBoxVectors()
        except:
            raise ValueError("Failed to load CRYST1 information")
    print("Using this cell:\n   " + str(boxa) +
          "\n   " + str(boxb) + "\n   " + str(boxc))
    return (boxa, boxb, boxc)


def get_coords(inp):
    if 'bincoordinates' in inp:
        print(f"Reading positions from NAMDBin: "+inp.bincoordinates)
        coords = NAMDBin(inp.bincoordinates).getPositions()
    else:
        import warnings
        print(f"Reading positions from PDB: ")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pdb = get_pdb(inp.coordinates)
        coords = pdb.positions
    return coords


def plumed_parser(fn):
    # Hack to workaround https://github.com/openmm/openmm-plumed/pull/27
    # Joins continuation lines and strips comments
    out = []
    continuing = False
    with open(fn) as f:
        for l in f:
            l = l.strip()
            if "#" in l:
                l = l[:l.find("#")]  # Strip comments
            dots = l.find("...")
            if not continuing:
                if dots == -1:
                    out.append(l)
                else:
                    out.append(l[:dots])
                    continuing = True
            else:
                if dots == -1:
                    out[-1] = out[-1] + " " + l
                else:
                    out[-1] = out[-1] + " " + l[:dots]
                    continuing = False
    return "\n".join(out)


def every(T, t):
    f = T/t
    assert f.is_integer() is True
    return int(f)


def get_best_platform():
    num = Platform.getNumPlatforms()
    pp_s = {}
    for i in range(num):
        pp = Platform.getPlatform(i)
        pn = pp.getName()
        ps = float(pp.getSpeed())
        pp_s[pn] = ps

    so = sorted(pp_s.items(), key=lambda x: x[1], reverse=True)
    sel = "->"
    for i in so:
        print(sel, i[0], i[1])
        sel = "  "
    return so[0][0]


def check_openmm():
    version = Platform.getOpenMMVersion()
    plugindir = Platform.getDefaultPluginsDirectory()
    print(" OpenMM details:")
    print("  Version     : OpenMM " + str(version))
    print("  Plugin dir  : " + str(plugindir))

    # Try loading the plugins and checking for errors
    Platform.loadPluginsFromDirectory(plugindir)
    errs = Platform.getPluginLoadFailures()
    if len(errs):
        print("Some errors were found loading plugins. Some platforms may not be available: \n")
    for x in errs:
        print(x)


def add_reporters(simulation, basename, log_every, save_every,
                  total_steps, continuing, checkpoint_file):
    print(
        f"Reporting every {log_every} steps and checkpointing on {checkpoint_file} every {save_every} steps.")

    fp=open(f"{basename}.log", "a" if continuing else "w")
    simulation.reporters.append(app.DCDReporter(f"{basename}.dcd", save_every,
                                                append=continuing, enforcePeriodicBox=False))
    simulation.reporters.append(app.CheckpointReporter(checkpoint_file,
                                                       save_every))
    simulation.reporters.append(StdoutLogReporter(log_every, total_steps))
    simulation.reporters.append(app.StateDataReporter(fp,
                                                      log_every,
                                                      step=True,
                                                      time=True,
                                                      potentialEnergy=True,
                                                      kineticEnergy=True,
                                                      totalEnergy=True,
                                                      temperature=True,
                                                      volume=True,
                                                      progress=True,
                                                      remainingTime=True,
                                                      speed=True,
                                                      totalSteps=total_steps,
                                                      separator='\t'))
