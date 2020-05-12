import simtk.openmm.openmm


def remove_barostat(system):
    """Remove MonteCarloBarostat if present"""
    fs = system.getForces()
    for i,f in enumerate(fs):
        if type(f) == simtk.openmm.openmm.MonteCarloBarostat or \
           type(f) == simtk.openmm.openmm.MonteCarloAnisotropicBarostat:
            system.removeForce(i)
            return

def every(T,t):
    f = T/t
    assert f.is_integer() is True
    return int(f)


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

