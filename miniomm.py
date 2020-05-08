

import sys
import os.path

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit


def main(runtype="NVT"):

    system_basename = "structure"

    output_basename = "output"

    dt = 5.0 * unit.femtosecond
    T = 300 * unit.kelvin
    energyFreq = 1 * unit.picosecond
    trajectoryFreq = 100 * unit.picosecond
    # restartFreq = trajectoryFreq
    equilibrationTime = 100 * unit.picosecond
    runTime = 1 * unit.nanosecond

    checkpoint_file = "restart.chk"
    nonbondedCutoff = 1.0*unit.nanometers
    frictionCoefficient = 1 / unit.picosecond

    pdb = app.PDBFile(f'{system_basename}.pdb')
    # forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    """
    system = forcefield.createSystem(pdb.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0*unit.nanometers,
                                     constraints=app.AllBonds,
                                     rigidWater=True,
                                     ewaldErrorTolerance=0.0005)
    """

    if os.path.exists(f"{system_basename}.prmtop"):
        prmtop = app.AmberPrmtopFile(f"{system_basename}.prmtop")
        system = prmtop.createSystem(nonbondedMethod=app.PME,
                                     nonbondedCutoff=nonbondedCutoff,
                                     constraints=app.AllBonds)
        topology = prmtop.topology
    elif os.path.exists(f"{system_basename}.psf"):
        psf = app.CharmmPsfFile('{system_basename}.psf')
        params = app.CharmmParameterSet("parameters")
        system = psf.createSystem(params,
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=nonbondedCutoff,
                                  constraints=app.AllBonds)
    else:
        raise FileNotFoundError("No AMBER nor CHARMM run files found")

    integrator = mm.LangevinIntegrator(T,
                                       frictionCoefficient,
                                       dt)
    integrator.setConstraintTolerance(0.00001)

    """
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    """

    simulation = app.Simulation(topology, system, integrator)
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

    simulation.reporters.append(app.DCDReporter(f"{output_basename}.dcd",
                                                energyFreq/dt,
                                                append=True))
    simulation.reporters.append(app.CheckPointReporter(checkpoint_file,
                                                       trajectoryFreq/dt))
    simulation.reporters.append(app.StateDataReporter(sys.stdout,
                                                      trajectoryFreq/dt,
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
    simulation.reporters.append(app.StateDataReporter(f"{output_basename}.log",
                                                      trajectoryFreq/dt,
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
    return


if __name__ == "__main__":
    main()
