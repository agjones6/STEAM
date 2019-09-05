# %% File to read in the output from the STEAM code
import numpy as np
import matplotlib.pyplot as plt
import os
import dataclasses
from dataclasses import dataclass
# import os.path
from os import path
import csv
import pandas as pd

# Functions that are to be used
def find_index(name_matrix,desired_string):
    # This function is made to search the matrix of strings and find the index of it
    for i in range(len(name_matrix)):
        if name_matrix[i] == desired_string:
            break

    return i

def find_data(data_matrix,name_matrix,desired_string):
    return data_matrix[find_index(name_matrix, desired_string)]

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
allData.file = []

# User Input
myStr = "rho_CR"

# This is the number of lines to skip when pulling in data
skip_lines_num = 30

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
    allData.file.append(csv.reader(open(cFileString)))

    # This variable is used to pull only as much data as is desired
    sk_count = skip_lines_num
    for spltLine in allData.file[-1]:



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
        elif (sk_count >= skip_lines_num): # This only allows the data to be pulled in if
            # Current Line's data
            cData = []

            # Going through each of the comma seperated values to fill the current data matrix
            for i in range(len(spltLine)):
                if (spltLine[i] != ""):
                    cData.append(np.array(float(spltLine[i].strip())))

            # If the data is a matrix, it is handled differently
            if (allData.type[-1] == "matrix"):
                # If this is the first time the matrix has been touched
                #   (first time step pulled in) the the based numpy array is created
                # Else the numpy array is built
                if (count_line == 2):
                    allData.data.append(np.array(cData))
                else:
                    allData.data[-1] = np.vstack((allData.data[-1],cData))

            # Else the data is handled as a single value
            else:
                # If this is the first time the matrix has been touched
                #   (first time step pulled in) the the based numpy array is created
                # Else the numpy array is built
                if (count_line == 2):
                    for i in range(len(cData)):
                        allData.data.append(np.array(cData[i]))

                else:
                    i2 = len(cData)
                    for i in range(len(cData)):
                        allData.data[i-i2] = np.vstack((allData.data[i-i2],cData[i]))

            # Resetting the skip counter
            sk_count = 1

        else:
            sk_count = sk_count + 1

        count_line = count_line + 1
        # if (count_line > 1e10):
        #     break


# %% Plotting

# %% Finding the desired names
print(allData.name)

# %%
time = find_data(allData.data,allData.name,"el_time")


desired_data = find_data(allData.data,allData.name,myStr)
plt.plot(time, desired_data)
plt.title(myStr)
plt.show()
