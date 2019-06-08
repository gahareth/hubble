import math
import matplotlib.pyplot as plt
import sys

AccelerationDueToGravity = 9.8

class PhysicalQuantity:
    def __init__(self, InitialValue, InitialDerivative):
        self.value = InitialValue
        self.derivative = InitialDerivative

def Euler(value, step, function, time = 0):
    initialValue = value.value
    value.value += step * value.derivative
    value.derivative += step * function(time, value.derivative, initialValue)

def EulerFirstOrder(value, step, function, time = 0):
    value += step * function(time, value)

def RK4SecondOrder(value, step, function, time = 0):
    derivativeCoefficients = []
    valueCoefficients = []

    derivativeCoefficients.append(step * function(time, value.derivative, value.value))
    valueCoefficients.append(step * value.derivative)

    for i in range(2):
        derivativeCoefficients.append(step * function(time + step / 2, value.derivative + derivativeCoefficients[-1] / 2, value.value + valueCoefficients[-1] / 2)) 
        valueCoefficients.append(step * (value.derivative + derivativeCoefficients[-1] / 2))

    derivativeCoefficients.append(step * function(time + step, value.derivative + derivativeCoefficients[-1], value.value + valueCoefficients[-1]))
    valueCoefficients.append(step * (value.derivative + derivativeCoefficients[-1]))

    value.derivative += (derivativeCoefficients[0] + 2 * derivativeCoefficients[1] + 2 * derivativeCoefficients[2] + derivativeCoefficients[3]) / 6
    value.value += (valueCoefficients[0] + 2 * valueCoefficients[1] + 2 * valueCoefficients[2] + valueCoefficients[3]) / 6
 
def RK4FirstOrder(value, step, function, time = 0):
    valueCoefficients = []
    valueCoefficients.append(step * function(time, value))
    for i in range(2):
        valueCoefficients.append(step * function(time + step / 2, value + valueCoefficients[-1] / 2))
    valueCoefficients.append(step * function(time + step, value + valueCoefficients[-1]))
    value += (valueCoefficients[0] + 2 * valueCoefficients[1] + 2 * valueCoefficients[2] + valueCoefficients[3]) / 6
    return value

