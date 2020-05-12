

import sys
import os.path
from optparse import OptionParser

from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as u

import util
import NAMDBin


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


def run_omm():

    system_basename = "structure"

    minimize_basename = "minimize"
    equil_basename = "equil"
    run_basename = "output"

    dt = options.timestep * u.femtosecond
    temperature = 300 * u.kelvin
    pressure = 1 * u.atmospheres
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

    req_platform = Null
    if options.platform is not Null:
        req_platform = mm.Platform.getPlatformByName(options.platform)

    req_properties = {}
    if options.device is not None and platform.getName() in ('CUDA', 'OpenCL'):
        req_properties['DeviceIndex'] = options.device
    if options.precision is not None and platform.getName() in ('CUDA', 'OpenCL'):
        req_properties['Precision'] = options.precision

    system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

    integrator = mm.LangevinIntegrator(temperature,
                                       frictionCoefficient,
                                       dt)
    integrator.setConstraintTolerance(1e-5)

    simulation = app.Simulation(topology, system, integrator, req_platform, req_properties)
    ctx = simulation.context
    print(f"Using platform {platform.getName()} with properties:")
    for prop in platform.getPropertyNames():
        print(f"    {prop}\t\t{platform.getPropertyValue(ctx,prop)}")

    resuming = False
    try:
        with open(checkpoint_file, "rb") as cf:
            ctx.loadCheckpoint(cf.read())
        print(f"Successfully loaded {checkpoint_file}")
        resuming = True
    except FileNotFoundError:
        print(f"Starting from initial positions")
        pdb = app.PDBFile(f'{system_basename}.pdb')
        ctx.setPositions(pdb.positions)

    log_every = util.every(energyFreq,dt)
    save_every = util.every(trajectoryFreq,dt)

    # -------------------------------------------------------
    if not os.path.exists(f"{minimize_basename}.xml"):
        print('Minimizing...')
        simulation.minimizeEnergy()
        simulation.saveState(f"{minimize_basename}.xml")
    else:
        print(f"{minimize_basename}.xml exists, skipping minimization.")

    # -------------------------------------------------------
    if not os.path.exists(f"{equil_basename}.xml"):
        print('Equilibrating...')
        ctx.setVelocitiesToTemperature(temperature)
        neq = util.every(equilibrationTime, dt)
        util.add_reporters(simulation, equil_basename,
                           log_every, save_every, neq, resuming)
        simulation.step(int(neq))
        simulation.saveState(f"{equil_basename}.xml")
    else:
        print(f"{equil_basename}.xml exists, skipping equilibration.")

    # -------------------------------------------------------
    if not os.path.exists(f"{run_basename}.xml"):
        print('Running Production...')
        nrun = util.every(runTime, dt)
        simulation.reporters = []
        util.remove_barostat(system)
        util.add_reporters(simulation, run_basename,
                           log_every, save_every, nrun, resuming)
        simulation.step(int(nrun))
        simulation.saveState(f"{run_basename}.xml")
    else:
        print(f"{run_basename}.xml exists, skipping run.")

    print('Done!')
    return


def main():
    run_omm()


if __name__ == "__main__":
    global options

    parser = OptionParser()
    platformNames = [mm.Platform.getPlatform(
        i).getName() for i in range(mm.Platform.getNumPlatforms())]
    parser.add_option('--platform', dest='platform',
                      choices=platformNames, help='name of the platform to benchmark')
    parser.add_option('--pme-cutoff', default='0.9', dest='cutoff',
                      type='float', help='direct space cutoff for PME in nm [default: 0.9]')
    parser.add_option('--timestep', default='5', dest='timestep',
                      type='float', help='integration timestep in fs [default: 4.0]')
    parser.add_option('--hours', default='11.5', dest='run_hours', type='float',
                      help='target simulation length in hours [default: 11.5]')
    parser.add_option('--heavy-hydrogens', action='store_true', default=True,
                      dest='heavy', help='repartition mass to allow a larger time step')
    parser.add_option('--device', default=None, dest='device',
                      help='device index for CUDA or OpenCL')
    parser.add_option('--precision', default='single', dest='precision', choices=('single', 'mixed',
                                                                                  'double'),
                      help='precision mode for CUDA or OpenCL: single, mixed, or double [default: single]')
    (options, args) = parser.parse_args()

    main()
