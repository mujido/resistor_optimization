import numpy as np
import random
import sys
import traceback
from cuckoo.filter import BCuckooFilter
import struct

e12Bases = np.array([1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2])

e12 = np.array([10**e * e12Bases for e in range(6)]).flatten()
e12.flags.writeable = False

sortedE12 = np.sort(e12)
sortedE12.flags.writeable = False

targetPotentials = np.array([
    0.0722,
    0.2126,
    0.2848,
    0.7874,
    0.9278,
    1,
    1.0722,
    1.2126,
    1.2848,
    1.7152,
    1.7874,
    1.9278,
    2
])
targetPotentials.flags.writeable = False

# Vo = Vi * Rb / (Rb + Rt)
# Rb / (Rb + Rt) = Vo / Vi
# (Rb + Rt) / Rb = Vi / Vo
# 1 + (Rt / Rb) = Vi / Vsh
# Rt / Rb = Vi / Vo - 1

sourceVoltage = 5.0
targetRatios = sourceVoltage / targetPotentials - 1.0

def parallel(resistors):
    return 1.0 / np.sum(resistors)

def bestRatios(desiredRatios, resistors):
    count = len(resistors)
    allOn = 2**count - 1
    assignedRatios = np.zeros(len(desiredRatios))

    for i in range(allOn + 1):
        if i == 0:
            ratio = 0
        elif i == allOn:
            ratio = 1
        else:
            top = 0
            bottom = 0

            for j in range(count):
                if i & (1 << j) != 0:
                    top += 1.0 / resistors[j]
                else:
                    bottom += 1.0 / resistors[j]

            ratio = top / bottom

        for j, dr in enumerate(desiredRatios):
            if np.abs(dr - ratio) < np.abs(assignedRatios[j] - dr):
                assignedRatios[j] = ratio

    return assignedRatios

def checkAccuracy(desiredRatios, resistors):
    assignedRatios = bestRatios(desiredRatios, resistors)
    # return targetRatios - assignedRatios
    return np.sum(np.abs(100.0 * (targetRatios - assignedRatios) / targetRatios))
    # return assignedRatios

searchRange = np.append(np.arange(-6, 0), np.arange(1, 7))
# searchRange = np.arange(1, 4)
# searchRange = [-1, 1]
# searchRange = np.arange(-73, 74)

def findBiggestChange(desiredRatios, resistors, currentScore):
    bestIndex = None
    bestResistor = None
    bestScore = currentScore

    for (i, r) in enumerate(resistors):
        r = resistors[i]
        c = np.argmin(np.abs(sortedE12 - r))
        rcopy = np.copy(resistors)

        for d in searchRange:
            k = c + d
            if k < 0:
                continue
            elif k >= len(sortedE12):
                break

            adj = sortedE12[k]
            # if (i > 0 and adj < resistors[i - 1]) or (i + 1 < len(resistors) and adj > resistors[i + 1]):
            #     break

            rcopy[i] = adj

            crBytes = struct.pack(structFormat, *rcopy)
            if not seenFilter.contains(crBytes):
                seenFilter.insert(crBytes)
            else:
                # print("collision")
                continue

            acc = checkAccuracy(desiredRatios, rcopy)
            if acc < bestScore:
                bestIndex = i
                bestResistor = adj
                bestScore = acc

    if bestResistor is not None:
        rcopy = np.copy(resistors)
        rcopy[bestIndex] = bestResistor
        return (np.sort(rcopy), bestScore)
    else:
        return (resistors, currentScore)

seenFilter = BCuckooFilter(capacity=int(100e6), error_rate=1e-5)

# bestResistors = np.array([330000, 68000, 1800, 12000, 22000, 6800, 4700])
# startingResistors = np.array([150000,  68000,   1800,  27000,  18000,   6800,   4700])
# startingResistors = np.array([470, 6800, 560, 2200, 5600, 150, 56000])
# startingResistors = np.array([220, 2200, 180, 3.9, 56, 18, 12])
# startingResistors = np.array([3.3, 10.0, 15.0, 47.0, 150.0, 180.0, 1500.0])
startingResistors = random.sample(e12.tolist(), 6)
# startingResistors = np.array([4.7, 68.0, 15.0, 180.0, 18.0])
# startingResistors = [27.0, 82.0, 120.0, 470.0, 1200.0, 2200.0]

startingResistors = np.array(startingResistors, dtype=np.double)
startingResistors.sort()
startingScore = checkAccuracy(targetRatios, startingResistors)

bestResistors = np.copy(startingResistors)
bestResistors.flags.writeable = False
bestScore = startingScore

currentResistors = bestResistors
currentScore = bestScore
randomizations = 0

structFormat = f"{len(startingResistors)}d"

def searchSortedRandom():
    global currentResistors
    global currentScore
    global bestResistors
    global bestScore

    (newResistors, newScore) = findBiggestChange(targetRatios, currentResistors, currentScore)
    if newScore < currentScore:
        currentResistors = newResistors
        currentScore = newScore
        # print((currentScore, currentResistors))
        print(f"[{randomizations}] Score: {currentScore:7.2f} ({currentScore - bestScore:+7.2f}) best: {bestScore:7.2f}")
    else:
        if currentScore < bestScore:
            bestResistors = currentResistors
            bestScore = currentScore

        currentResistors = np.copy(bestResistors)

        for i in range(2):
            idx = random.randint(0, len(currentResistors) - 1)
            # ma = np.ma.masked_values(sortedE12, currentResistors[idx])
            ma = np.ma.masked_array(sortedE12)
            ma = np.ma.masked_less_equal(ma, currentResistors[idx - 1]) if idx > 0 else ma
            ma = np.ma.masked_greater_equal(ma, currentResistors[idx + 1]) if idx + 1 < len(currentResistors) else ma
            domain = ma.compressed()
            if domain.any():
                nr = random.choice(domain)
                currentResistors[idx] = nr

try:
    while True:
        (newResistors, newScore) = findBiggestChange(targetRatios, currentResistors, currentScore)
        if newScore < currentScore:
            currentResistors = newResistors
            currentScore = newScore
            # print((currentScore, currentResistors))

            print(f"[{randomizations}] Score: {currentScore:7.2f} ({currentScore - bestScore:+7.2f}) best: {bestScore:7.2f}")
        else:
            if currentScore < bestScore:
                bestResistors = currentResistors
                bestScore = currentScore

            while True:
                currentResistors = np.copy(bestResistors)

                randCount = random.randrange(1, len(currentResistors))
                randIndices = random.sample(range(len(currentResistors)), randCount)
                for idx in randIndices:
                    nr = random.choice(e12)
                    currentResistors[idx] = nr

                currentResistors.sort()
                crBytes = struct.pack(structFormat, *currentResistors)
                if not seenFilter.contains(crBytes):
                    seenFilter.insert(crBytes)
                    break
                else:
                    print(f"[{randomizations}] (collision) best: {bestScore:7.2f}")

            currentResistors.flags.writeable = False

            # currentResistors.sort()
                    
                # diffs = sortedE12 - currentResistors[idx] lowerBound = currentResistors[idx - 1] if idx > 0 else 0
                # upperBound = currentResistors[idx + 1] if idx + 1 < len(currentResistors) else 
                # if idx > 0:
                #     lowerBound = np.ma.masked_greater_equal(, 0.0).argmin()
                # upperBound = np.ma.masked_less_equal(diffs, 0.0).argmax()
                # if diffs[lowerBound] is np.ma.masked:
                #     lowerBound = idx
                # if diffs[upperBound] is np.ma.masked:
                #     upperBound = idx - 1
                # if upperBound - lowerBound > 0:
                    # nr = random.choice(sortedE12[lowerBound:upperBound + 1])
                    # currentResistors[idx] = nr

            randomizations += 1
            currentScore = 1000000
except Exception as e:
    traceback.print_exc()
except KeyboardInterrupt:
    pass

print(targetRatios)

br = bestRatios(targetRatios, bestResistors)
print(f"Score: {bestScore} ({bestScore - startingScore})")
print(br)
print(100 * (targetRatios - br) / targetRatios)
print(bestResistors.tolist())
