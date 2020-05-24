from simtk.openmm import app
import sys

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
    for i,f in enumerate(fs):
        if type(f) == simtk.openmm.openmm.MonteCarloBarostat or \
           type(f) == simtk.openmm.openmm.MonteCarloAnisotropicBarostat:
            system.removeForce(i)
            return

def every(T,t):
    f = T/t
    assert f.is_integer() is True
    return int(f)


def add_reporters(simulation, basename, log_every, save_every,
                  total_steps, continuing, checkpoint_file):
    print(f"Reporting every {log_every} steps and checkpointing on {checkpoint_file} every {save_every} steps.")
    simulation.reporters.append(app.DCDReporter(f"{basename}.dcd", save_every,
                                                append=continuing))
    simulation.reporters.append(app.CheckpointReporter(checkpoint_file,
                                                       save_every))
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

