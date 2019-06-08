
import math
import matplotlib.pyplot as plt
import sys
from Core import RK4FirstOrder
from Core import Euler 
from Core import PhysicalQuantity

metresInParsec = 3.08567758149137e16
parsecsInMpc = 1e6
distanceMpc = 1e2
speedOfLight = 299792458

def metresToMpc(metres):
    return metres / (metresInParsec * parsecsInMpc)

def hubbleToSi(hubble):
    return hubble / (metresInParsec * parsecsInMpc)

def MpcToMetres(distanceMpc):
    return distanceMpc * metresInParsec * parsecsInMpc

def yearsToSeconds(years):
    return years * 365.25 * 60 * 60 * 24

def secondsToYears(seconds):
    return seconds / (60 * 60 * 24 * 365.25)

class Universe():
    def __init__(self, name):
        self.name = name
        
    def Hubble(self, Time):
        return 1

class StaticUniverse(Universe):
    def __init__(self):
        Universe.__init__(self, "Static Hubble parameter")
        HubbleConstant = 70e3
        self.HubbleConstantSiUnits = hubbleToSi(HubbleConstant)

    def Hubble(self, Time):
        return self.HubbleConstantSiUnits
        
class LinearIncreaseWithTimeUniverse(Universe):
    def __init__(self):
        Universe.__init__(self, "Hubble parameter increasing linearly with time")
        HubbleConstantNow = 70e3
        HubbleConstantSiUnitsNow = hubbleToSi(HubbleConstantNow)
        self.HubbleConstantTimeZero = hubbleToSi(1e3)
        TimeNow = 13.8e9
        TimeNowSeconds = yearsToSeconds(TimeNow)
        self.HubbleConstantRateOfIncrease = (HubbleConstantSiUnitsNow - self.HubbleConstantTimeZero) / TimeNowSeconds

    def Hubble(self, Time):
        return self.HubbleConstantRateOfIncrease * Time + self.HubbleConstantTimeZero

class MovingObject():
    def __init__(self, initialDisplacement, initialVelocity, universe):
        self.Displacement = initialDisplacement
        self.Velocity = initialVelocity 
        self.Universe = universe

    def RateOfChangeOfDisplacement(self, Time, Displacement):
        return self.Universe.Hubble(Time) * Displacement - self.Velocity

    def Update(self, TimeStep, Time):
        self.Displacement = RK4FirstOrder(self.Displacement, TimeStep, lambda time, Displacement: self.RateOfChangeOfDisplacement(time, Displacement), Time) 

universes = {'static' : StaticUniverse(), 'linear' : LinearIncreaseWithTimeUniverse()}

def SimulateConstantVelocityTravel(InitialDisplacement, hubbleModel, testNumber):
    speed = speedOfLight
    Displacements = []
    DistancesFromStart = []
    StartAndEndSeparation = []
    Timestamps = []

    # Try to use a sensible duration to stop the computer running out of memory.
    TimeStep = (2 * InitialDisplacement / speedOfLight) / 1e6

    time = 0
    universe = universes[hubbleModel]
    movingObject = MovingObject(InitialDisplacement, speedOfLight, universe)
    startingPoint = MovingObject(0, -speedOfLight, universe)
    terminationCondition = False
    while terminationCondition == False:
        movingObject.Update(TimeStep, time)
        startingPoint.Update(TimeStep, time)
        Timestamps.append(secondsToYears(time))    
        time += TimeStep
        Displacements.append(metresToMpc(movingObject.Displacement))
        DistancesFromStart.append(metresToMpc(startingPoint.Displacement))
        StartAndEndSeparation.append(metresToMpc(movingObject.Displacement) + metresToMpc(startingPoint.Displacement))
        if len(Displacements) > 1:
            terminationCondition = Displacements[-2] < Displacements[-1]
        if movingObject.Displacement < 0:
            terminationCondition = True

    plt.figure(testNumber)
    plt.plot(Timestamps, Displacements, label='distance to destination')
    # plt.plot(Timestamps, DistancesFromStart, label='distance from start')
    plt.plot(Timestamps, StartAndEndSeparation, label='distance from start to end')
    plt.xlabel('time (years)')
    plt.ylabel('displacement (Mpc)')
    plt.title(universe.name)
    plt.legend()
    plt.show()

def main(argv):
    distance = float(argv[1])
    InitialDisplacement = MpcToMetres(distance) 
    
    testNumber = 1
    SimulateConstantVelocityTravel(InitialDisplacement, 'static', testNumber)
    testNumber += 1
    SimulateConstantVelocityTravel(InitialDisplacement, 'linear', testNumber)

if __name__ == "__main__":
    main(sys.argv)
