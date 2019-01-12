import re
import csv

p = re.compile("[-+]?\d*\.\d*e*[-+]*\d*|\d+")

def readFromFileChFig1():
    data = []
    with open("TremorData/ChFig1.txt", "r") as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            data.append([float(i) for i in row])
    data = list(zip(*data))

    return data

def readFromFileChFig3():
    # Bi-axial hand tremor, RH Sampling rate: 100/sec
    # x: RH measuring P+/S deg/sec
    # y: RH measuring F+/E deg/sec
    # z: RL measuring ankle F+/E, deg/sec
    # a: LL measuring ankle F/E+, deg/sec

    data = []
    with open("TremorData/ChFig3.txt", "r") as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            data.append([float(i) for i in row])
    data = list(zip(*data))

    return data[1:]

def readFromFileChFig4():
    # Bi-axial hand tremor, RH Sampling rate: 100/sec
    # x: RH measuring P+/S deg/sec
    # y: RH measuring F+/E deg/sec

    data = []
    with open("TremorData/ChFig4.txt", "r") as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            data.append([float(i) for i in row])
    data = list(zip(*data))

    return data

def readFromFileChFig5():
    # Concurrent rest tremors, RH, LH.   Sampling rate: 100/sec
    # x: LH measuring P/S+   deg/sec
    # y: RH    measuring P+/S  deg/sec
    data = []
    with open("TremorData/ChFig5.txt", "r") as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            data.append([float(i) for i in row])
    data = list(zip(*data))

    return data



def readFromFileChFig6():
    # Bi-axial, bilateral records, RH, LH. LH rest tremor, RH voluntary P/S
    # x: RH measuring P+/S    deg/sec
    # y: LH measuring P+/S    deg/sec
    # z: RH measuring F+/E    deg/sec
    # a: LH measuring F/E+    deg/sec

    data = []
    with open("TremorData/ChFig6.txt", "r") as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            data.append([float(i) for i in row])
    data = list(zip(*data))

    return data[1:]


def readFromFileChFig7():
    # Bilateral postural tremor
    # x: RH   measuring F+/E   deg/sec
    # y: LH  measuring F/E+  deg/sec

    data = []
    with open("TremorData/ChFig7.txt", "r") as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            data.append([float(i) for i in row])
    data = list(zip(*data))

    return data[1:]

def readFromFileChFig8():
    # Action tremor LH  Accompanies voluntary LH P/S
    # x: RH measuring P+/S  deg/sec
    # y: LH measuring P+/S  deg/sec

    data = []
    with open("TremorData/ChFig8.txt", "r") as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            data.append([float(i) for i in row])
    data = list(zip(*data))

    return data[1:]

def readFromAllFiles():

    allFiles = []
    allMeasurments = []

    allFiles.append(readFromFileChFig1())
    allFiles.append(readFromFileChFig3())
    allFiles.append(readFromFileChFig4())
    allFiles.append(readFromFileChFig5())
    allFiles.append(readFromFileChFig6())
    allFiles.append(readFromFileChFig7())
    allFiles.append(readFromFileChFig8())

    for f in allFiles:
        for m in f:
            allMeasurments.append(m)

    return allMeasurments