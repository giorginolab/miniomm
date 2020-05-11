from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

prmtop = AmberPrmtopFile('structure.prmtop')
pdb = PDBFile('structure.pdb')

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)

platform = simulation.context.getPlatform()

print(f"Using platform {platform.getName()} with properties...")
for prop in platform.getPropertyNames():
    print(f"    {prop}\t{platform.getPropertyValue(simulation.context,prop)}")

simulation.context.setPositions(pdb.positions)

simulation.saveState("0-initial.xml")
    
simulation.minimizeEnergy()

simulation.saveState("1-minimized.xml")

simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(5000)

simulation.saveState("2-postrun.xml")

