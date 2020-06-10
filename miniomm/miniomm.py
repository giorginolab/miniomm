

import sys
import os.path

from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as u

from miniomm.config import Config
import miniomm.util as util
from miniomm.namdbin import NAMDBin
from miniomm.namdxsc import write_xsc


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
    run_steps = int(inp.run)
    basename = inp.trajectoryfile.replace(".xtc", "")

    nonbondedCutoff = float(inp.cutoff) * u.angstrom
    switchDistance = float(inp.switchdistance) * u.angstrom
    frictionCoefficient = float(inp.thermostatdamping) / u.picosecond

    endTime = run_steps * dt


    util.check_openmm()

    if options.platform is None:
        print("Selecting best platform:")
        req_platform_name = util.get_best_platform()
    else:
        print(f"Requesting platform {options.platform}")
        req_platform_name = options.platform
    req_platform = mm.Platform.getPlatformByName(req_platform_name)

    req_properties = {}    # {'UseBlockingSync':'true'}
    if options.device is not None and 'DeviceIndex' in req_platform.getPropertyNames():
        print("    Setting DeviceIndex = "+options.device)
        req_properties['DeviceIndex'] = options.device
    if options.precision is not None and 'Precision' in req_platform.getPropertyNames():
        print("    Setting Precision = "+options.precision)
        req_properties['Precision'] = options.precision


    if dt >= 2.5 * u.femtosecond:
        hmr = 4 * u.amu
        print(f"Enabling hydrogen mass repartitioning at mass(H) = {hmr}")
    else:
        hmr = 1 * u.amu

    if 'parmfile' in inp: 
        print(f"Creating an AMBER system...")
        if 'structure' in inp:
            print("Warning: 'structure' given but irrelevant for AMBER")
        prmtop = app.AmberPrmtopFile(inp.parmfile)
        system = prmtop.createSystem(nonbondedMethod=app.PME,
                                     nonbondedCutoff=nonbondedCutoff,
                                     switchDistance=switchDistance,
                                     constraints=app.AllBonds,
                                     hydrogenMass=hmr)
        topology = prmtop.topology
    else:
        print(f"Creating a CHARMM system...")
        psf = app.CharmmPsfFile(inp.structure)
        params = app.CharmmParameterSet(inp.parameters, permissive=True)
        psf.setBox( 50.*u.angstrom, 50.*u.angstrom, 50.*u.angstrom) # otherwise
                                                                    # refuses
                                                                    # PME
        system = psf.createSystem(params,
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=nonbondedCutoff,
                                  switchDistance=switchDistance,
                                  constraints=app.AllBonds,
                                  hydrogenMass=hmr)
        topology = psf.topology

    if 'barostat' in inp and inp.getboolean('barostat'):
        pressure = float(inp.barostatpressure) * u.bar
        print(f"Enabling barostat at {pressure}...")
        system.addForce(mm.MonteCarloBarostat(pressure, temperature))

    if 'plumedfile' in inp:
        print("Attempting to load PLUMED plugin...")
        from openmmplumed import PlumedForce
        plines = util.plumed_parser(inp.plumedfile)
        system.addForce(PlumedForce(plines))

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
        with open(checkpoint_file, 'rb') as cf:
            ctx.loadCheckpoint(cf.read())
        # ctx.loadCheckpoint(str(checkpoint_file))
        print(f"Successfully loaded {checkpoint_file}, resuming simulation...")
        resuming = True

    else:
        print(f"File {checkpoint_file} absent, starting simulation from the beginning...")
        coords = util.get_coords(inp)
        ctx.setPositions(coords)

    if not resuming:
        (boxa, boxb, boxc) = util.get_box_size(inp)
        ctx.setPeriodicBoxVectors(boxa, boxb, boxc)

        if 'minimize' in inp:
            print(f'Minimizing for max {inp.minimize} iterations...')
            simulation.minimizeEnergy(maxIterations=int(inp.minimize))
            simulation.saveState(f"miniomm_minimized.xml")
        else:
            if 'binvelocities' in inp:
                print(f"Reading velocities from NAMDBin: "+inp.binvelocities)
                vels = NAMDBin(inp.binvelocities).getVelocities()
                ctx.setVelocities(vels)
            else:
                print(f"Resetting thermal velocities at {temperature}")
                ctx.setVelocitiesToTemperature(temperature)

    # -------------------------------------------------------
    print("\n")

    if basename != "output":
        print(f"Warning: basename is '{basename}' instead of 'output'")

    startTime = ctx.getState().getTime()
    startTime_f = startTime.in_units_of(u.nanoseconds).format("%.3f")
    endTime_f = endTime.in_units_of(u.nanoseconds).format("%.3f")
    remaining_steps = round((endTime-startTime)/dt)
    print(f"Current simulation time is {startTime_f},"+
          f"running up to {endTime_f}")
    print(f"Running for {remaining_steps} timesteps = {remaining_steps * dt.in_units_of(u.nanosecond)}...")
    print("\n")


    log_every = util.every(energyFreq, dt)
    save_every = util.every(trajectoryFreq, dt)

    util.add_reporters(simulation, basename,
                       log_every, save_every,
                       remaining_steps, resuming, checkpoint_file)

    # ----------------------------------------
    simulation.saveState(f"miniomm_pre.xml")

    # ----------------------------------------
    simulation.step(remaining_steps)

    # ----------------------------------------
    simulation.saveState(f"miniomm_post.xml")
    final_state = simulation.context.getState(getPositions=True,
                                              getVelocities=True)
    final_coor = final_state.getPositions(asNumpy=True)
    NAMDBin(final_coor).write_file(f"{basename}.coor")

    final_box = final_state.getPeriodicBoxVectors(asNumpy=True)
    write_xsc(f"{basename}.xsc", remaining_steps, final_state.getTime(), final_box)
    # ----------------------------------------
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

    (options, args) = parser.parse_args()

    if len(args) > 0:
        print("Remaining args: "+" ".join(args))

    print(util.getBanner())
    run_omm(options)


if __name__ == "__main__":
    main()
    
