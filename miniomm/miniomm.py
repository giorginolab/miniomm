

import sys
import os.path

from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as u

import miniomm.util as util
from miniomm.namdbin import NAMDBin
from miniomm.config import Config


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

checkpoint_file = "miniomm_restart.chk"


def run_omm(options):

    inp = Config(options.input)

    dt = float(inp.timestep) * u.femtosecond
    temperature = float(inp.temperature) * u.kelvin
    energyFreq = 1 * u.picosecond
    trajectoryFreq = int(inp.trajectoryperiod) * dt
    # restartFreq = trajectoryFreq
    # equilibrationTime = 10 * u.picosecond
    nrun = int(inp.run)
    basename = inp.trajectoryfile.replace(".xtc", "")

    nonbondedCutoff = float(inp.cutoff) * u.angstrom
    switchDistance = 7.5 * u.angstrom
    frictionCoefficient = float(inp.thermostatdamping) / u.picosecond


    util.check_openmm()

    if options.platform is None:
        print("Selecting best platform:")
        req_platform_name = util.get_best_platform()
    else:
        print(f"Requesting platform {options.platform}")
        req_platform_name = options.platform
    req_platform = mm.Platform.getPlatformByName(req_platform_name)

    req_properties = {}
    if options.device is not None and 'DeviceIndex' in req_platform.getPropertyNames():
        print("    Setting DeviceIndex = "+options.device)
        req_properties['DeviceIndex'] = options.device
    if options.precision is not None and 'Precision' in req_platform.getPropertyNames():
        print("    Setting Precision = "+options.precision)
        req_properties['Precision'] = options.precision


    if dt >= 4 * u.femtosecond:
        hmr = 4 * u.amu
        print(f"Enabling hydrogen mass repartitioning at mass(H) = {hmr}")
    else:
        hmr = 1 * u.amu

    if 'parmfile' in inp:
        prmtop = app.AmberPrmtopFile(inp.parmfile)
        system = prmtop.createSystem(nonbondedMethod=app.PME,
                                     nonbondedCutoff=nonbondedCutoff,
                                     switchDistance=switchDistance,
                                     constraints=app.AllBonds,
                                     hydrogenMass=hmr)
        topology = prmtop.topology
    else:
        psf = app.CharmmPsfFile(inp.structure)
        params = app.CharmmParameterSet(inp.parameters, permissive=True)
        system = psf.createSystem(params,
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=nonbondedCutoff,
                                  switchDistance=switchDistance,
                                  constraints=app.AllBonds,
                                  hydrogenMass=hmr)
    if 'barostat' in inp and inp.getboolean('barostat'):
        pressure = float(inp.barostatpressure) * u.bar
        print(f"Enabling barostat at {pressure}")
        system.addForce(mm.MonteCarloBarostat(pressure, temperature))

    # plumed

    integrator = mm.LangevinIntegrator(temperature,
                                       frictionCoefficient,
                                       dt)
    integrator.setConstraintTolerance(1e-5)

    simulation = app.Simulation(topology, system, integrator,
                                req_platform, req_properties)
    ctx = simulation.context
    platform = ctx.getPlatform()
    print(f"Got platform {platform.getName()} with properties:")
    for prop in platform.getPropertyNames():
        print(f"    {prop}\t\t{platform.getPropertyValue(ctx,prop)}")
    print("\n")

    resuming = False
    if os.path.exists(checkpoint_file):
        print(f"Attempting to reload {checkpoint_file}")
        with open(checkpoint_file, 'rb') as cf:
            ctx.loadCheckpoint(cf.read())
        # ctx.loadCheckpoint(str(checkpoint_file))
        print(f"Successfully loaded {checkpoint_file}, resuming simulation")
        # ctx.setTime(100000)
        resuming = True

    else:
        if 'bincoordinates' in inp:
            print(f"Reading positions from NAMDBin: "+inp.bincoordinates)
            coords = NAMDBin(inp.bincoordinates).getPositions()
        else:
            import warnings
            print(f"Reading positions from PDB: ")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pdb = app.PDBFile(inp.coordinates)
            coords = pdb.positions
        ctx.setPositions(coords)

    if not resuming:
        if 'extendedsystem' in inp:
            print("Reading box size from "+inp.extendedsystem)
            (boxa, boxb, boxc) = util.parse_xsc_units(inp.extendedsystem)
        elif 'boxsize' in inp:
            print("Using boxsize from input string "+inp.boxsize)
            (boxa, boxb, boxc) = util.parse_boxsize_units(inp.boxsize)
        else:
            print("Last resort: PDB CRYST1...")
            try:
                (boxa, boxb, boxc) = pdb.topology.getPeriodicBoxVectors()
            except:
                raise ValueError("Failed to load CRYST1 information")

        print("Using this cell:\n   " + str(boxa) +
              "\n   " + str(boxb) + "\n   " + str(boxc))
        ctx.setPeriodicBoxVectors(boxa, boxb, boxc)

        if 'minimize' in inp:
            print(f'Minimizing for max {inp.minimize} iterations...')
            simulation.minimizeEnergy(maxIterations=int(inp.minimize))
            simulation.saveState(f"minimized.xml")
        else:
            if 'binvelocities' in inp:
                print(f"Reading velocities from NAMDBin: "+inp.binvelocities)
                vels = NAMDBin(inp.binvelocities).getVelocities()
                ctx.setVelocities(vels)
            else:
                print(f"Resetting thermal velocities at {temperature}")
                ctx.setVelocitiesToTemperature(temperature)

    # -------------------------------------------------------
    print("\n\n")
    print(
        f'Running for {nrun} timesteps = {nrun * dt.in_units_of(u.nanosecond)}...')

    # nrun = util.every(equilibrationTime, dt)
    log_every = util.every(energyFreq, dt)
    save_every = util.every(trajectoryFreq, dt)

    util.add_reporters(simulation, basename,
                       log_every, save_every,
                       nrun, resuming, checkpoint_file)
    simulation.step(nrun)
    simulation.saveState(f"output.xml")

    print('Done!')
    return


def main():
    from optparse import OptionParser
    parser = OptionParser()
    platformNames = [mm.Platform.getPlatform(
        i).getName() for i in range(mm.Platform.getNumPlatforms())]
    parser.add_option('--input', dest='input',
                      default="input", help='name of the input file')
    parser.add_option('--platform', dest='platform',
                      choices=platformNames, help='name of the platform to benchmark')
    parser.add_option('--device', default=None, dest='device',
                      help='device index for CUDA or OpenCL')
    parser.add_option('--precision', dest='precision', choices=('single', 'mixed', 'double'),
                      help='precision mode for CUDA or OpenCL: single, mixed, or double')

    parser.add_option('--hours', default='11.5', dest='run_hours', type='float',
                      help='target simulation length in hours [default: 11.5]')

    parser.add_option('--pme-cutoff', default='0.9', dest='cutoff',
                      type='float', help='direct space cutoff for PME in nm [default: 0.9]')
    parser.add_option('--timestep', default='5', dest='timestep',
                      type='float', help='integration timestep in fs [default: 4.0]')
    parser.add_option('--heavy-hydrogens', action='store_true', default=True,
                      dest='heavy', help='repartition mass to allow a larger time step')
    (options, args) = parser.parse_args()

    if len(args) > 0:
        print("Remaining args: "+" ".join(args))

    print(util.getBanner())
    run_omm(options)


if __name__ == "__main__":
    main()
    
