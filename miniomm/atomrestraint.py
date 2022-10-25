import math
import os

import openmm
from openmm import app
import openmm.unit as u


def atrest_parser(txt, run_steps, dt):
    TIME_UNITS = { "us" : u.microsecond,
                   "ns" : u.nanosecond,
                   "ps" : u.picosecond,
                   "fs" : u.femtosecond }

    atrest_dict = {
        "selection" : None,
        "axes" : None,
        "setpoints" : None,
        "n_atoms" : 0
    }

    txt = txt.lower()
    opening_quotes_idx = txt.find('"')
    closing_quotes_idx = txt.find('"', opening_quotes_idx + 1)
    atrest_dict["selection"] = txt[opening_quotes_idx + 1 : closing_quotes_idx]

    tokens = txt[closing_quotes_idx+1:].split()
    for token_idx, token in enumerate(tokens):

        if token == "width":
            print("WARNING: Parameter 'width' in restraint selection is not available. Ignoring.")
            continue

        if token == "axes":
            atrest_dict["axes"] = tokens[token_idx+1]

        if token == "setpoints":
            setpoints_list = []
            for setpoint in tokens[token_idx+1:]:
                if "@" in setpoint:
                    setpoint_data = {}
                    force, time = setpoint.split("@")
                    setpoint_data["force"] = int(force)

                    if time[-2:] in TIME_UNITS.keys():
                        time_size = int(time[:-2])
                        time_unit = TIME_UNITS[time[-2:]]
                        time_in_fs = (time_size * time_unit).value_in_unit(u.femtosecond)
                        setpoint_data["step"] = math.floor(time_in_fs / dt.value_in_unit(u.femtosecond))
                    elif time[-2:].isalpha():
                        raise Exception("ERROR: Restraint time unit not recognised.")
                    else:
                        setpoint_data["step"] = int(time)

                    if setpoint_data["step"] > run_steps:
                        raise Exception("ERROR: Restraint time larger than simulation time")

                    setpoint_data["percent"] = math.floor(setpoint_data["step"] / run_steps * 1000)
                    setpoints_list.append(setpoint_data) 
                else: # if no "@" in token
                    break  
            atrest_dict["setpoints"] = setpoints_list

    # Checks
    if atrest_dict["axes"] == None:
        atrest_dict["axes"] = "xyz"
        print("WARNING: Using default restraint axes = 'xyz'.")
    if atrest_dict["setpoints"] == None:
        raise Exception("ERROR: No restraint setpoints declared.")

    return atrest_dict


def add_restraints(simulation, atrest_dict):
    is_x = is_y = is_z = False
    if "x" in atrest_dict["axes"]: is_x = True
    if "y" in atrest_dict["axes"]: is_y = True
    if "z" in atrest_dict["axes"]: is_z = True

    restraint_formula = f'k*periodicdistance({"x, " if is_x else ""}{"y, " if is_y else ""}{"z, " if is_z else ""}' \
                                         + f'{"x0, " if is_x else ""}{"y0, " if is_y else ""}{"z0, " if is_z else ""}'
    restraint_formula = restraint_formula[:-2] + ")^2" # Remove trailing comma-and-space and end the equation

    restraint = openmm.CustomExternalForce(restraint_formula)
    simulation.system.addForce(restraint)
    restraint.addGlobalParameter('k', 1*u.kilocalories_per_mole/u.angstrom**2) # Dummy value
    if is_x: restraint.addPerParticleParameter('x0')
    if is_y: restraint.addPerParticleParameter('y0')
    if is_z: restraint.addPerParticleParameter('z0')

    atom_positions = simulation.context.getState(getPositions=True).getPositions()
    
    with open("tmp.pdb", "w") as fh:
        app.PDBFile.writeFile(simulation.topology, atom_positions, file=fh)
    import mdtraj as md
    tmp_pdb = md.load("tmp.pdb")

    try:
        selection_list = tmp_pdb.topology.select(atrest_dict["selection"])
    except:
        raise Exception("ERROR: restraint selection not recognised. "
        + "Please, read MDtraj selection reference: https://mdtraj.org/1.9.4/atom_selection.html")

    os.remove("tmp.pdb")

    for atom in simulation.topology.atoms():
        if atom.index in selection_list:
            restraint.addParticle(atom.index, atom_positions[atom.index])
            atrest_dict["n_atoms"] += 1

    simulation.context.reinitialize(preserveState=True)


def get_starting_force_and_gradient(current_percent, atrest_dict):
    min_distance = 1001
    previous_setpoint = None

    # Get closest previous setpoint from current percent
    for setpoint in atrest_dict["setpoints"]:
        distance_to_setpoint = current_percent - setpoint["percent"]
        if setpoint["percent"] <= current_percent and distance_to_setpoint <= min_distance:
            min_distance = distance_to_setpoint
            previous_setpoint = setpoint

    force = previous_setpoint["force"]
    gradient = get_force_gradient(previous_setpoint, atrest_dict)
    starting_force = force + gradient * min_distance

    return starting_force, gradient


def get_force_gradient(current_setpoint, atrest_dict):
    min_distance = 1001
    next_setpoint = None

    # Get closest next setpoint from current setpoint
    for setpoint in atrest_dict["setpoints"]:
        distance_to_setpoint = setpoint["percent"] - current_setpoint["percent"]
        if setpoint["percent"] > current_setpoint["percent"] and distance_to_setpoint < min_distance:
            min_distance = distance_to_setpoint
            next_setpoint = setpoint
    
    if next_setpoint == None: # If no more setpoints, gradient is 0
        next_setpoint = current_setpoint

    gradient = (next_setpoint["force"] - current_setpoint["force"]) / min_distance
    return gradient
    