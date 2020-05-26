from simtk.unit import *
from simtk.openmm import *
from simtk.unit import nanometer
from simtk.unit import angstrom
from simtk.unit import picosecond
from simtk.unit import femtosecond


import numpy as np
import time
import math
import os
import sys

# Taken from https://github.com/compsciencelab/MDy/blob/master/python/mdy/reporters.py


class StdoutLogReporter:
    def __init__(self, reportInterval, totalSteps):
        self._reportInterval = reportInterval
        self._inited = False
        self._totalSteps = totalSteps
        self._lastvol = None

    def print_headers(self):
        print("  %10s %10s %11s %11s %11s %8s %8s %8s %8s %11s" % (
            "Step", "Time", "PE", "KE", "Total E", "Temp", "Volume", "Fluct.", "ISpeed", "Completion"))
        print("  %10s %10s %11s %11s %11s %8s %8s %8s %8s %11s" % (
            "", "ns", "kJ/mol", "kJ/mol", "kJ/mol", "K", "nm^3", "%", "ns/day", "dd:hh:mm:ss"))

    def _init(self, simulation, system, state):
        # Compute the number of degrees of freedom.
        dof = 0
        for i in range(system.getNumParticles()):
            if system.getParticleMass(i) > 0 * unit.dalton:
                dof += 3
        dof -= system.getNumConstraints()
        if any(type(system.getForce(i)) == CMMotionRemover for i in range(system.getNumForces())):
            dof -= 3
        self._dof = dof
        self._initialSteps = simulation.currentStep
        self._initialClockTime = time.time()
        self._initialSimulationTime = state.getTime()
        self._lastClockTime = self._initialClockTime
        self._lastSimulationTime = self._initialSimulationTime
        self._inited = True

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, True)

    def report(self, simulation, state):
        if not self._inited:
            self._init(simulation, simulation.system, state)
            self.print_headers()
        clockTime = time.time()
        step = simulation.currentStep
        timex = state.getTime().value_in_unit(nanosecond)
        pe = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        te = pe + ke
        temp = (2 * state.getKineticEnergy() / (self._dof *
                                                unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)

        speed = -1.0
        elapsedSeconds = clockTime - self._initialClockTime
        if elapsedSeconds > 0:
            elapsedDays = (elapsedSeconds) / 86400.0
            elapsedSteps = simulation.currentStep - self._initialSteps
            elapsedNs = (
                state.getTime() - self._initialSimulationTime).value_in_unit(unit.nanosecond)
            speed = elapsedNs / elapsedDays

        instaSpeed = -1.0
        instaSeconds = clockTime - self._lastClockTime
        if instaSeconds > 0:
            instaDays = instaSeconds / 86400.0
            instaNs = (state.getTime() -
                       self._lastSimulationTime).value_in_unit(unit.nanosecond)
            instaSpeed = instaNs / instaDays
            self._lastClockTime = clockTime
            self._lastSimulationTime = state.getTime()

        if elapsedSteps:
            estimatedTotalSeconds = (
                self._totalSteps - self._initialSteps) * elapsedSeconds / elapsedSteps
        else:
            estimatedTotalSeconds = 0
        remainingSeconds = int(estimatedTotalSeconds - elapsedSeconds)
        remainingDays = remainingSeconds // 86400
        remainingSeconds -= remainingDays * 86400
        remainingHours = remainingSeconds // 3600
        remainingSeconds -= remainingHours * 3600
        remainingMinutes = remainingSeconds // 60
        remainingSeconds -= remainingMinutes * 60
        remianingString = "--"
        if remainingDays > 0:
            remainingString = "%d:%d:%02d:%02d" % (
                remainingDays, remainingHours, remainingMinutes, remainingSeconds)
        elif remainingHours > 0:
            remainingString = "%d:%02d:%02d" % (
                remainingHours, remainingMinutes, remainingSeconds)
        elif remainingMinutes > 0:
            remainingString = "%d:%02d" % (remainingMinutes, remainingSeconds)
        else:
            remainingString = "0:%02d" % remainingSeconds

        box = state.getPeriodicBoxVectors()
        v = box[0][0] * box[1][1] * box[2][2]
        volume = v.value_in_unit(nanometer ** 3)
        if self._lastvol:
            fluctuation = 100. * (volume / self._lastvol - 1.)
        else:
            fluctuation = 0.
        self._lastvol = volume

        print("# %10ld %10.3f %11.2f %11.2f %11.2f %8.2f %8.2f %8.2f %8.2f %11s" % (
            step, timex, pe, ke, te, temp, volume, fluctuation, instaSpeed, remainingString))

        if (math.isnan(pe) or math.isnan(ke) or math.isnan(temp)):
            raise ValueError("Simulation has become unstable. Aborted")
