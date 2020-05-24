from simtk.openmm import app
import sys
from simtk.openmm.openmm import Platform
from reporters import *


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


def parse_xsc(xsc):
    with open(xsc) as f:
        for l in f:
            if '#' in l:
                continue
            ls = l.split(" ")
            if len(ls) >= 10:
                boxx = float(ls[1])
                boxy = float(ls[5])
                boxz = float(ls[9])
    return (boxx, boxy, boxz)


def remove_barostat(system):
    """Remove MonteCarloBarostat if present"""
    fs = system.getForces()
    for i, f in enumerate(fs):
        if type(f) == simtk.openmm.openmm.MonteCarloBarostat or \
           type(f) == simtk.openmm.openmm.MonteCarloAnisotropicBarostat:
            system.removeForce(i)
            return


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
        Platform.loadPluginsFromDirectory( plugindir )
        errs = Platform.getPluginLoadFailures()
        if len(errs):
            print("Some errors were found loading plugins. Some platforms may not be available: \n")
        for x in errs:
            print(x)


def add_reporters(simulation, basename, log_every, save_every,
                  total_steps, continuing, checkpoint_file):
    print(
        f"Reporting every {log_every} steps and checkpointing on {checkpoint_file} every {save_every} steps.")
    simulation.reporters.append(app.DCDReporter(f"{basename}.dcd", save_every,
                                                append=continuing))
    simulation.reporters.append(app.CheckpointReporter(checkpoint_file,
                                                       save_every))
    """                                                       
    simulation.reporters.append(app.StateDataReporter(sys.stdout,
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
    """
    simulation.reporters.append(StdoutLogReporter(log_every, total_steps))
    simulation.reporters.append(app.StateDataReporter(f"{basename}.log",
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
