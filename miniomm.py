

import sys
from sys import stdout

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit


def main():

    prmtop_file="system.prmtop"
    coord_file="input.coor"

    out = "output"
    output_file=f"{out}.dcd"
    log_file=f"{out}.log"
    checkpoint_file=f"{out}.chk"
    
    dt = 5.0 * unit.femtosecond
    T = 300 * unit.kelvin
    energyFreq = 1 * unit.picosecond
    trajectoryFreq = 100 * unit.picosecond
    restartFreq = trajectoryFreq
    equilibrationTime = 100 * unit.picosecond
    runTime = 1 * unit.nanosecond

    nonbondedCutoff = 1.0*unit.nanometers
    frictionCoefficient = 1 / unit.picosecond

    pdb = app.PDBFile('input.pdb')
    #forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    """
    system = forcefield.createSystem(pdb.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0*unit.nanometers,
                                     constraints=app.AllBonds,
                                     rigidWater=True,
                                     ewaldErrorTolerance=0.0005)
    """

    prmtop = AmberPrmtopFile(prmtop_file)
    system = prmtop.createSystem(nonbondedMethod=app.PME,
                                 nonbondedCutoff=nonbondedCutoff,
                                 constraints=app.AllBonds)
    
    integrator = mm.LangevinIntegrator(T,
                                       frictionCoefficient,
                                       dt)
    integrator.setConstraintTolerance(0.00001)

    """
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    """

    simulation = app.Simulation(pdb.topology, system, integrator)
    # ... , platform, properties

    try:
        with open(checkpoint_file, "rb") as cf:
            simulation.context.loadCheckpoint(cf.read())
        print(f"Successfully loaded {checkpoint_file}")
    except FileNotFoundError:
        print(f"Starting from initial positions")
        simulation.context.setPositions(pdb.positions)

    print('Minimizing...')
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(T)
    print('Equilibrating...')
    simulation.step(equilibrationTime/dt)

    simulation.reporters.append(app.DCDReporter(output_file,
                                                energyFreq/dt,
                                                append=True))
    simulation.reporters.append(app.CheckPointReporter(checkpoint_file,
                                                       trajectoryFreq/dt)
    simulation.reporters.append(app.StateDataReporter(stdout, trajectoryFreq/dt,
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
                                                      totalSteps=runTime/dt,
                                                      separator='\t'))
    simulation.reporters.append(app.StateDataReporter(log_file, trajectoryFreq/dt,
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
                                                      totalSteps=runTime/dt,
                                                      separator='\t'))

    print('Running Production...')
    simulation.step(runTime/dt)

    print('Done!')



 if __name__ == "__main__":
     main()


