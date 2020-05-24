

import sys
import os.path
from optparse import OptionParser

from simtk.openmm import app
import simtk.openmm as mm
import simtk.unit as u

import util
from namdbin import NAMDBin
from config import Config


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



def run_omm(options, inp):

    dt = float(inp.timestep) * u.femtosecond
    temperature = float(inp.temperature) * u.kelvin
    energyFreq = 1 * u.picosecond
    trajectoryFreq = int(inp.trajectoryperiod) * dt
    # restartFreq = trajectoryFreq
    # equilibrationTime = 10 * u.picosecond
    nrun = int(inp.run)

    nonbondedCutoff = float(inp.cutoff) * u.angstrom
    frictionCoefficient = float(inp.thermostatdamping) / u.picosecond
    hmr = 4 * u.amu

    log_every = util.every(energyFreq,dt)
    save_every = util.every(trajectoryFreq,dt)



    if 'parmfile' in inp:
        prmtop = app.AmberPrmtopFile(inp.parmfile)
        system = prmtop.createSystem(nonbondedMethod=app.PME,
                                     nonbondedCutoff=nonbondedCutoff,
                                     constraints=app.AllBonds,
                                     hydrogenMass=hmr)
        topology = prmtop.topology
    else:
        psf = app.CharmmPsfFile(inp.structure)
        params = app.CharmmParameterSet(inp.parameters, permissive = True)
        system = psf.createSystem(params,
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=nonbondedCutoff,
                                  constraints=app.AllBonds,
                                  hydrogenMass=hmr)
    

    req_platform = None
    if options.platform is not None:
        req_platform = mm.Platform.getPlatformByName(options.platform)

    req_properties = {}
    if options.device is not None:
        req_properties['DeviceIndex'] = options.device
    if options.precision is not None:
        req_properties['Precision'] = options.precision

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
    print(f"Using platform {platform.getName()} with properties:")
    for prop in platform.getPropertyNames():
        print(f"    {prop}\t\t{platform.getPropertyValue(ctx,prop)}")

    
    resuming = False
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, "rb") as cf:
            ctx.loadCheckpoint(cf.read())
        print(f"Successfully loaded {checkpoint_file}")
        resuming = True
        
    else:
        if 'bincoordinates' in inp:
            print(f"Reading NAMDBin positions from "+inp.bincoordinates)
            coords = NAMDBin(inp.bincoordinates).getPositions()
        else:
            print(f"Reading PDB positions")
            pdb = app.PDBFile(inp.coordinates)
            coords = pdb.positions
        ctx.setPositions(coords)

    if not resuming:
        if 'extendedsystem' in inp:
            print("Reading box size from "+inp.extendedsystem)
            box = util.parse_xsc(inp.extendedsystem)
        elif 'boxsize' in inp:
            print("Using boxsize from input")
            box = [float(x) * u.angstrom for x in inp.boxsize.split(" ")]
        else:
            try:
                box = pdb.topology.getPeriodicBoxVectors()
            except:
                print("Last resort PDB CRYST failed")
                raise
            
        boxa = mm.Vec3( box[0], 0., 0.  ) * u.angstrom
        boxb = mm.Vec3( 0., box[1],  0. ) * u.angstrom
        boxc = mm.Vec3( 0., 0., box[2]  ) * u.angstrom

        print("Using this cell: "+ str(boxa) + " " + str(boxb) + " " + str(boxc))
        ctx.setPeriodicBoxVectors(boxa, boxb, boxc)
        

    
    # -------------------------------------------------------
    if 'minimize' in inp and not resuming:
        print('Minimizing...')
        simulation.minimizeEnergy()
        simulation.saveState(f"minimized.xml")

    # -------------------------------------------------------
    if 'binvelocities' in inp:
        print("binvelocities not supported, randomizing")
    
    print(f'Running for {nrun} timesteps = {nrun * dt.in_units_of(u.nanosecond)}...')
    ctx.setVelocitiesToTemperature(temperature)
    
    # nrun = util.every(equilibrationTime, dt)
    util.add_reporters(simulation, inp.trajectoryfile,
                       log_every, save_every, nrun, resuming, checkpoint_file)
    simulation.step(nrun)
    simulation.saveState(f"output.xml")

    print('Done!')
    return


def main(options):
    inp = Config(options.input)
    run_omm(options, inp)


    

if __name__ == "__main__":
#    global options

    parser = OptionParser()
    platformNames = [mm.Platform.getPlatform(
        i).getName() for i in range(mm.Platform.getNumPlatforms())]
    parser.add_option('--input', dest='input',
                      default="input", help='name of the input file')
    parser.add_option('--platform', dest='platform',
                      choices=platformNames, help='name of the platform to benchmark')
    parser.add_option('--device', default=None, dest='device',
                      help='device index for CUDA or OpenCL')
    parser.add_option('--precision', default='single', dest='precision', choices=('single', 'mixed',
                                                                                  'double'),
                      help='precision mode for CUDA or OpenCL: single, mixed, or double [default: single]')

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

    main(options)
