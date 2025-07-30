import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
import csv
import pandas as pd

# All in meters
def heatCooling(busbarWidth, busbarLength, busbarThickness, CFM, heat):
    carAir = 30 # C
    busBarSurfaceTemp = 55 # C 
    velocityAir = 10 # m/s
    airDensity = 1.145 # kg/s
    prandtlNumber = 0.7268 # For 35C air
    airThermalConductivity = 0.027 # W/mK
    airViscosity = 1.65 * 10**-5 # m^2/s
    reynolds = (busbarLength * velocityAir) / airViscosity
    nusseltsNumber = 0.0308*(reynolds**(0.8)) * (prandtlNumber**(1/3))
    convectionHeatTC = (airThermalConductivity * nusseltsNumber) / busbarLength
    surfaceArea = busbarLength * busbarWidth + 2 * busbarLength * busbarThickness + 2* busbarWidth * busbarThickness
    heatTransfer = convectionHeatTC * surfaceArea * (busBarSurfaceTemp - carAir)
    return heatTransfer


df = pd.read_csv('AllData.csv')
dfRegen = pd.read_csv('RegenData.csv')

currentMultiplier = 1.02 # New pack is lower voltage

timeRegen = dfRegen.iloc[:, 0].tolist()

averageRPM = []
for j in range(len(timeRegen)):
    averageRPM.append(((-1 * dfRegen.iloc[j, 1]) + dfRegen.iloc[j, 2] + (-1 * dfRegen.iloc[j, 3]) + dfRegen.iloc[j, 4])/4)

brakePressure = dfRegen.iloc[:, 5].tolist()
regenVoltage = dfRegen.iloc[:, 7].tolist()

# Convert individual columns to lists
time = df.iloc[:, 0].tolist()  # First column
current = np.multiply(currentMultiplier, df.iloc[:, 1].tolist())  # Second column
voltage = df.iloc[:, 2].tolist()  # Third column
cells = 128
cellIR = 0.0016
packResistance = cells * cellIR * 1.33 # ohms
packWeight = 32 # kg
cellCp = 717 
dt = 0.005 # seconds
startTemp = 35 # C

airCp = 1007 
airDensity = 1.145
allowableTDelta = 20 # C

coolingHeat = heatCooling(0.017, 0.0300, 0.002, 0, 0) # W
print(coolingHeat)

tabTempFactor = 1.35 # Tabs get hotter than body

brakeThreshold = 25
rpmHigh = 0
rpmLow = 0
regenCurrentProfile = []
b = 0

timeRegening = 0
timeRegeningList = []
powerLossRegening = []
powerRegening = []

regenResistnaceIncrease = 1.5 # Multiplier of pack ressistnace increase when regening

while b < len(brakePressure):
    while brakePressure[b] > brakeThreshold:

        rpmHigh = averageRPM[b]
        currentRegen = -5*10**-12*rpmHigh**3 + 3*10**-7*rpmHigh**2 - 0.0045 * rpmHigh + 7.2698
        if currentRegen > 0:
            regenCurrentProfile.append(current[b])
        else:
            regenCurrentProfile.append(currentRegen * 2)
            timeRegening = dt + timeRegening
            timeRegeningList.append(timeRegening)
            powerLossRegening.append(((currentRegen * 2)**2 * (packResistance * regenResistnaceIncrease)/1000) * dt)
            powerRegening.append(regenVoltage[b] * abs(currentRegen * 2)/1000) # kW

        b = b + 1
    regenCurrentProfile.append(current[b])
    b = b + 1

#plt.plot(timeRegen, regenCurrentProfile)
#plt.plot(timeRegen, current)
#plt.show()

averageCurrent = 60 # np.average(current)
packHeat = averageCurrent**2 * packResistance

tempSpikesRegen = []
tempSpikes = (np.multiply(np.multiply(current, current), packResistance))/(packWeight * cellCp) * tabTempFactor
tempCooled = (coolingHeat * cells)/(packWeight * cellCp)

airVolume = 0.0254 * 0.620 * 0.385 #m^3
airMass = airVolume * airDensity

tempRise = []
cooledPack = []
cooledPackRegen = []
tempRiseRegen = []
prevNoCooling = 0
prevCooling = 0
prevRegen = 0
prevRegenCooling = 0
addingHeat= 0
cfmDesire = []
heatInAir = 0
batteryWattsLoss = []
airTemp = []

for i in range(len(time)):
    currentNow = regenCurrentProfile[i]
    heatInAir = (coolingHeat + heatInAir)
    airTemperature = ((heatInAir / (airMass * airCp)) + startTemp) * dt
    airTemp.append(airTemperature)

    if currentNow > 0:
        tempSpikesRegen.append((currentNow * currentNow * packResistance)/(packWeight * cellCp) * tabTempFactor)
        batteryWattsLoss.append(regenCurrentProfile[i] * regenCurrentProfile[i] * packResistance)
        addingHeat = (regenCurrentProfile[i] * regenCurrentProfile[i] * packResistance)
        
    else:
        tempSpikesRegen.append((currentNow * currentNow * packResistance * regenResistnaceIncrease)/(packWeight * cellCp) * tabTempFactor)
        batteryWattsLoss.append(0 * packResistance)
        addingHeat = (currentNow * currentNow * packResistance * regenResistnaceIncrease)


    prevRegen = tempSpikesRegen[i] + prevRegen
    prevRegenCooling = tempSpikesRegen[i] - tempCooled + prevRegenCooling
    prevNoCooling = tempSpikes[i] + prevNoCooling
    prevCooling = tempSpikes[i] - tempCooled + prevCooling
    massFlowCFM = (voltage[i] * currentNow * 3.41)/(1.08 * (airTemperature - startTemp))
    cfmDesire.append(massFlowCFM)
    tempRiseRegen.append(startTemp + prevRegen * dt)
    cooledPackRegen.append(startTemp + prevRegenCooling * dt)
    tempRise.append(startTemp + prevNoCooling * dt)
    cooledPack.append(startTemp + (prevCooling) * dt)

plt.plot(time, tempRise, label="No Regen Without Cooling")
plt.plot(time, cooledPack, label="No Regen With Cooling")
plt.plot(time, tempRiseRegen, label="Regen withtout Cooling")
plt.plot(time, cooledPackRegen, label="Regen With Cooling")

plt.title("Temperature Rise Endurance w/o Heat Fins")
plt.xlabel("Time (Seconds)")
plt.ylabel("Temperature (C)")

plt.legend()

print(addingHeat)

print(sum(powerLossRegening))

print(trapezoid(powerRegening, timeRegeningList)/3600)

print(packHeat)

massFlow = packHeat/(airCp * allowableTDelta) # kg/s
massFlowCFM = massFlow/airDensity * 35.3147 * 60

print(massFlowCFM)

averageHeatGen = integral = trapezoid((np.multiply(np.multiply(current, current), packResistance)), time)/time[-1]

print(averageHeatGen)

plt.show()
