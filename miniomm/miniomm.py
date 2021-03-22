

import sys
import os.path
import socket
import datetime

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


def _printPluginInfo():
    lp = mm.version.openmm_library_path
    print(f"""
           $OPENMM_CUDA_COMPILER: {os.environ.get('OPENMM_CUDA_COMPILER','(Undefined)')}
       OpenMM Library Path (...): {lp}
                  Loaded Plugins: """)
    for p in  mm.pluginLoadedLibNames:
        print("                                  "+p.replace(lp,"..."))
    print("            Loaded Plugin errors:")
    for e in mm.Platform.getPluginLoadFailures():
        print("                                  "+e)
    print("\n")



def run_omm(options):

    inp = Config(options.input)

    dt = float(inp.getWithDefault("timestep", 4)) * u.femtosecond
    temperature = float(inp.getWithDefault("temperature", 300)) * u.kelvin
    thermostattemperature = float(inp.getWithDefault("thermostattemperature", 300)) * u.kelvin
    logPeriod = 1 * u.picosecond
    trajectoryPeriod = int(inp.getWithDefault("trajectoryperiod", 25000)) * dt
    run_steps = int(inp.run)
    basename = "output"
    trajectory_file = basename + ".dcd"

    if 'PME' in inp and not inp.getboolean('PME'):
        nonbondedMethod = app.NoCutoff
    else:
        nonbondedMethod = app.PME
    nonbondedCutoff = float(inp.getWithDefault("cutoff", 9.0)) * u.angstrom
    switchDistance = float(inp.getWithDefault("switchdistance", 7.5)) * u.angstrom
    frictionCoefficient = float(inp.getWithDefault("thermostatdamping", 0.1)) / u.picosecond

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

    # Same logic as https://software.acellera.com/docs/latest/acemd3/reference.html
    if dt > 2 * u.femtosecond:
        hydrogenMass = 4 * u.amu
        constraints = app.AllBonds
        rigidWater = True
    elif dt > 0.5 * u.femtosecond:
        hydrogenMass = None
        constraints = app.HBonds
        rigidWater = True
    else:
        hydrogenMass = None
        constraints = None
        rigidWater = False


    print(f"""
                            Host: {socket.gethostname()}
                            Date: {datetime.datetime.now().ctime()}
                        Timestep: {dt}
                     Constraints: {constraints}
                     Rigid water: {rigidWater}
                       Nonbonded: {nonbondedMethod}
    Hydrogen mass repartitioning: {hydrogenMass}
    """)
    _printPluginInfo()


    if 'parmfile' in inp: 
        print(f"Creating an AMBER system...")
        if 'structure' in inp:
            print("Warning: 'structure' given but ignored for AMBER")
        prmtop = app.AmberPrmtopFile(inp.parmfile)
        system = prmtop.createSystem(nonbondedMethod=nonbondedMethod,
                                     nonbondedCutoff=nonbondedCutoff,
                                     switchDistance=switchDistance,
                                     constraints=constraints,
                                     hydrogenMass=hydrogenMass,
                                     rigidWater=rigidWater)
        topology = prmtop.topology
    else:
        print(f"Creating a CHARMM system...")
        psf = app.CharmmPsfFile(inp.structure)
        params = app.CharmmParameterSet(inp.parameters, permissive=False)
        psf.setBox( 50.*u.angstrom, 50.*u.angstrom, 50.*u.angstrom) # otherwise
                                                                    # refuses
                                                                    # PME
        system = psf.createSystem(params,
                                  nonbondedMethod=nonbondedMethod,
                                  nonbondedCutoff=nonbondedCutoff,
                                  switchDistance=switchDistance,
                                  constraints=constraints,
                                  hydrogenMass=hydrogenMass,
                                  rigidWater=rigidWater)
        topology = psf.topology

    if 'barostat' in inp and inp.getboolean('barostat'):
        pressure = float(inp.barostatpressure) * u.bar
        print(f"Enabling barostat at {pressure}...")
        system.addForce(mm.MonteCarloBarostat(pressure, thermostattemperature))

    if 'plumedfile' in inp:
        print("Attempting to load PLUMED plugin...")
        from openmmplumed import PlumedForce
        plines = util.plumed_parser(inp.plumedfile)
        system.addForce(PlumedForce(plines))

    integrator = mm.LangevinIntegrator(thermostattemperature,
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
    print("")

    resuming = False
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'rb') as cf:
            ctx.loadCheckpoint(cf.read())
        # ctx.loadCheckpoint(str(checkpoint_file))
        util.round_state_time(ctx, 10*dt)
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
    print("")
    inp.printWarnings()

    # -------------------------------------------------------
    print("")

    startTime = ctx.getState().getTime()
    startTime_f = startTime.in_units_of(u.nanoseconds).format("%.3f")
    endTime_f = endTime.in_units_of(u.nanoseconds).format("%.3f")
    remaining_steps = round((endTime-startTime)/dt)
    remaining_ns = remaining_steps * dt.value_in_unit(u.nanosecond)
    print(f"Current simulation time is {startTime_f}, running up to {endTime_f}.")
    print(f"Will run for {remaining_steps} timesteps = {remaining_ns:.3f} ns...")
    print("")

    log_every = util.every(logPeriod, dt)
    save_every = util.every(trajectoryPeriod, dt)
    if remaining_steps % save_every != 0:
        raise ValueError("Remaining steps is not a multiple of trajectoryperiod")

    util.add_reporters(simulation, trajectory_file,
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



