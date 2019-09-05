# %% File to read in the output from the STEAM code
import numpy as np
import matplotlib as mpl
import os
import dataclasses 
from dataclasses import dataclass
# import os.path
from os import path
import csv
import pandas as pd

# Allocating an object to store all of the data
@dataclass
class allData:
    type: any
    name: any
    data: any
    File: any

allData.type = []
allData.name = []
allData.data = []
allData.File = []

# Reading through every file in the output file director
outputFileLocation = "OutputFiles"
for filename in os.listdir(outputFileLocation):
    # Trying to open the files
    cFileString = outputFileLocation + "/" + filename
    if not path.exists(cFileString):
        print("file: ", filename, "didn't open")
        break
    elif not (".csv" in cFileString):
        print("Not a csv")
        break

    # Opening the file
    count_line = 1
    allData.File.append(csv.reader(open(cFileString)))
    # print(data)
    tst = []
    # next(data)
    # for row in data:
    #     tst.append(row)
    #
    #     count_line = count_line + 1
    #     if (count_line > 10):
    #         break
    #
    # print(tst)
    # exit()

    # # Reading the first line of the file to determine if it is the right type
    # count_line = 1
    for spltLine in allData.File[-1]:

        # This checks to see if the file is a csv that is expected. (Skips files not formatted correctly)
        # try:
        #     dumCheck = cLine.split(",")[1]
        # except:
        #
        #     break

        # Pulling the entire line split into sections
        # spltLine  = cLine.split(",")

        # Pulling the first word of the csv to decide how to handle the file
        if (count_line == 1):
            if (spltLine[0] == "MATRIX"):
                allData.name.append(spltLine[1])
                allData.type.append("matrix")
            else:
                for i in range(len(spltLine)):
                    if (spltLine[i] != ""):
                        allData.name.append(spltLine[i])
                        allData.type.append("value")
        else:
            cData = []

            for i in range(len(spltLine)):
                # print(np.array(float(spltLine[i].strip())))
                # print(spltLine[i])
                if (spltLine[i] != ""):
                    cData.append(np.array(float(spltLine[i].strip())))
            if (allData.type[-1] == "matrix"):
                if (count_line == 2):
                    allData.data.append(np.array(cData))
                else:
                    allData.data[-1] = np.vstack((allData.data[-1],cData))

        count_line = count_line + 1
        if (count_line > 100):
            break

    break


# print(allData.data[0][:][:].shape)
print(allData.data[0][:][:].shape)
