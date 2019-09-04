# %% File to read in the output from the STEAM code
import numpy as np
import matplotlib as mpl
import os
from dataclasses import dataclass

# Allocating an object to store all of the data
@dataclass
class allData:
    type: any
    name: any
    data: any

allData.type = []
allData.name = []
allData.data = []

# Reading through every file in the output file director
outputFileLocation = "OutputFiles"
for filename in os.listdir(outputFileLocation):
    # Trying to open the files
    try:
        cFile = open(outputFileLocation + "/" + filename)
        # with open(outputFileLocation + "/" + filename) as csv_file:
    except:
        print("file: ", filename, "didn't open")

    # Reading the first line of the file to determine if it is the right type
    count_line = 1
    for cLine in cFile:
        # This checks to see if the file is a csv that is expected. (Skips files not formatted correctly)
        try:
            dumCheck = cLine.split(",")[1]
        except:

            break

        # Pulling the entire line split into sections
        spltLine  = cLine.split(",")

        # Pulling the first word of the csv to decide how to handle the file
        if (count_line == 1):
            if (spltLine[0] == "MATRIX"):
                allData.name.append(spltLine[1])
                allData.type.append("matrix")
            else:
                for i in range(len(spltLine)):
                    if (spltLine[i] != "\n"):
                        allData.name.append(spltLine[i])
                        allData.type.append("value")
        else:
            cData = []
            for i in range(len(spltLine)):
                if (spltLine[i] != "\n"):
                    cData.append(float(spltLine[i]))
            if (allData.type[-1] == "matrix"):
                if (count_line == 2):
                    allData.data.append(np.array(cData))
                else:
                    allData.data[-1] = np.vstack((allData.data[-1],cData))

        count_line = count_line + 1
        if (count_line > 1e6):
            break
    break

print(allData.data[0][:][:].shape)
