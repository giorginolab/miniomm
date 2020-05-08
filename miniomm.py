

import sys
import os.path

from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as u


# Order of priority for the operation to perform
#   Minimization if minim.coor absent
#   NPT if equilib.coor absent
#   NVT otherwise

# Order of priority for box size
#   Checkpoint file
#   Box size (input.xsc)

# Order of priority for starting coordinates
#   Checkpoint file
#   Input coordinates (input.coor)
#   PDB coordinates (structure.pdb)


# Order of priority for starting velocities
#   Checkpoint file
#   Input velocities (input.vel)
#   Randomized velocities

checkpoint_file = "restart.chk"
    

def add_reporters(simulation, basename, log_every, save_every, total_steps, continuing):

    simulation.reporters.append(app.DCDReporter(f"{basename}.dcd", save_every, append=continuing))
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



def main(runtype="NVT"):

    system_basename = "structure"

    equil_basename = "equil"
    run_basename = "output"

    dt = 5.0 * u.femtosecond
    T = 300 * u.kelvin
    energyFreq = 1 * u.picosecond
    trajectoryFreq = 1 * u.picosecond
    # restartFreq = trajectoryFreq
    equilibrationTime = 10 * u.picosecond
    runTime = 1 * u.nanosecond

    nonbondedCutoff = 1.0*u.nanometers
    frictionCoefficient = 1 / u.picosecond
    hmr = 4 * u.amu

    # forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    """
    system = forcefield.createSystem(pdb.topology,
                                     nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0*u.nanometers,
                                     constraints=app.AllBonds,
                                     rigidWater=True,
                                     ewaldErrorTolerance=0.0005)
    """

    if os.path.exists(f"{system_basename}.prmtop"):
        prmtop = app.AmberPrmtopFile(f"{system_basename}.prmtop")
        system = prmtop.createSystem(nonbondedMethod=app.PME,
                                     nonbondedCutoff=nonbondedCutoff,
                                     constraints=app.AllBonds,
                                     hydrogenMass=hmr)
        topology = prmtop.topology
    elif os.path.exists(f"{system_basename}.psf"):
        psf = app.CharmmPsfFile('{system_basename}.psf')
        params = app.CharmmParameterSet("parameters")
        system = psf.createSystem(params,
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=nonbondedCutoff,
                                  constraints=app.AllBonds,
                                  hydrogenMass=hmr)
    else:
        raise FileNotFoundError("No AMBER nor CHARMM run files found")

    integrator = mm.LangevinIntegrator(T,
                                       frictionCoefficient,
                                       dt)
    integrator.setConstraintTolerance(1e-5)

    """
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
    """

    simulation = app.Simulation(topology, system, integrator)
    ctx = simulation.context
    # ... , platform, properties

    cont=False

    try:
        with open(checkpoint_file, "rb") as cf:
            ctx.loadCheckpoint(cf.read())
        print(f"Successfully loaded {checkpoint_file}")
        cont=True
    except FileNotFoundError:
        print(f"Starting from initial positions")
        pdb = app.PDBFile(f'{system_basename}.pdb')
        ctx.setPositions(pdb.positions)


    log_every = energyFreq/dt
    assert log_every.is_integer() == True
    log_every = int(log_every)

    save_every = trajectoryFreq/dt
    assert save_every.is_integer() == True
    save_every = int(save_every)
    
        
    print('Minimizing...')
    simulation.minimizeEnergy()

    
    print('Equilibrating...')
    ctx.setVelocitiesToTemperature(T)
    
    neq = equilibrationTime/dt
    assert neq.is_integer() == True

    add_reporters(simulation, equil_basename, log_every, save_every, neq, cont)
    simulation.step(int(neq))

    simulation.saveState(f"{equil_basename}.xml")

    


    print('Running Production...')
    # simulation.step(runTime/dt)

    print('Done!')
    return


if __name__ == "__main__":
    main()
