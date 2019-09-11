import numpy as np
import matplotlib.pyplot as plt
import os
import dataclasses
from dataclasses import dataclass
# import os.path
from os import path
import csv

# Functions that are to be used
def get_data(file_location, skip_lines_num):
    # Allocating an object to store all of the data
    @dataclass
    class allData:
        type: any
        name: any
        data: any
        file: any
    allData.type = []
    allData.name = []
    allData.data = []
    allData.file = []

    try:
        os.listdir(file_location)
        print("Reading Files from: ", file_location)
    except:
        # print("The output directory is not correct, trying to back out a directory")
        try:
            file_location = "../" + file_location
            os.listdir(file_location)
            print("Reading Files from: ", file_location)
        except:
            print("Output location cannot be found")
            exit()

    for filename in os.listdir(file_location):
        # Trying to open the files
        cFileString = file_location + "/" + filename
        if not path.exists(cFileString):
            print("file: ", filename, "didn't open trying to back out a directory ")
            if not path.exists("../" + cFileString):
                "The Folder can not be found "
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

    return allData

def get_data2(file_location, skip_lines_num, desired_string):
    # Allocating an object to store all of the data
    @dataclass
    class allData:
        type: any
        name: any
        data: any
        file: any
    allData.type = []
    allData.name = []
    allData.data = []
    allData.file = []

    try:
        os.listdir(file_location)
        print("Reading Files from: ", file_location)
    except:
        # print("The output directory is not correct, trying to back out a directory")
        try:
            file_location = "../" + file_location
            os.listdir(file_location)
            print("Reading Files from: ", file_location)
        except:
            print("Output location cannot be found")
            exit()

    for filename in os.listdir(file_location):
        # Trying to open the files
        cFileString = file_location + "/" + filename
        if not path.exists(cFileString):
            print("file: ", filename, "didn't open trying to back out a directory ")
            if not path.exists("../" + cFileString):
                "The Folder can not be found "
                break
        elif not (".csv" in cFileString):
            print("Not a csv")
            break

        # Opening the file
        count_line = 1
        allData.file.append(csv.reader(open(cFileString)))

        # print(count_line)
        # This variable is used to pull only as much data as is desired
        sk_count = skip_lines_num
        for spltLine in allData.file[-1]:

            # Pulling the first word of the csv to decide how to handle the file
            if (count_line == 1):
                incl_mat = np.zeros(len(spltLine))

                if (spltLine[0] == "MATRIX"):
                    for i2 in range(len(desired_string)):
                        if (desired_string[i2] in spltLine[1]):
                            allData.name.append(spltLine[1])
                            allData.type.append("matrix")
                            incl_mat[0] = 1
                    if incl_mat[0] != 1:
                        break
                else:
                    for i in range(len(spltLine)):
                        if (spltLine[i] != ""):
                            for i2 in range(len(desired_string)):
                                if (desired_string[i2] in spltLine[i]):
                                    allData.name.append(spltLine[i])
                                    allData.type.append("value")
                                    incl_mat[i] = 1

                print(incl_mat)
                    # exit()
            elif (sk_count >= skip_lines_num): # This only allows the data to be pulled in if
                # Current Line's data
                cData = []

                # Going through each of the comma seperated values to fill the current data matrix
                # if ()
                for i in range(len(spltLine)):
                    # if not desired_string in allData.name[-1]:
                    #     print(allData.name[-1])
                    #     continue

                    if (spltLine[i] != "") and (incl_mat[i] != 0):
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

    return allData

def find_index(name_matrix,desired_string):
    # This function is made to search the matrix of strings and find the index of it
    for i in range(len(name_matrix)):
        if name_matrix[i] == desired_string:
            break

    return i

def find_data(my_structure,desired_string):
    # This function is made to pull in data to a matrix
    data_array = []
    for t in range(len(desired_string)):
        data_array.append(my_structure.data[find_index(my_structure.name, desired_string[t])])

    return data_array

# %% Finding the desired names
def plot_data(my_structure,desired_string,*args):
    # This function is designed to plot the data in a desired manor
    #   if there is an optional word at the end of the function "matrix"

    # Pulling the time
    time = find_data(my_structure, ["el_time"])[0]

    # Getting the data that matches the input string
    desired_data = find_data(my_structure, desired_string)

    # Default values for the plot
    plot_type   = 0
    plot_title  = ""
    plot_xlabel = "time (hr)"
    plot_ylabel = ""
    plot_yscale = "linear"

    # Resolving the optional input arguments
    for t in range(len(args)):
        # Determing the plot type
        if (("type" in args[t]) and (len(args) > (t + 1))):
            if (args[t + 1] == "matrix"):
                plot_type = 1

        # Determing the plot title
        elif (("title" in args[t]) and (len(args) > (t + 1))):
            plot_title = args[t + 1]

        # Determing the plot xlabel
        elif (("xlabel" in args[t]) and (len(args) > (t + 1))):
            plot_xlabel = args[t + 1]

        # Determing the plot ylabel
        elif (("ylabel" in args[t]) and (len(args) > (t + 1))):
            plot_ylabel = args[t + 1]

        # Determing the plot yscale
        elif (("yscale" in args[t]) and (len(args) > (t + 1))):
            plot_yscale = args[t + 1]

    # The legend is populated differently if nodal/junction values are desired
    if (plot_type == 1):
        leg_str = range(1,len(desired_data[0][0,:])+1)
        if (plot_title == ""): plot_title = desired_string[0]
    else:
        leg_str = desired_string

    # Plotting the values
    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1)
    for i in range(len(desired_string)):
        plt.plot(time, desired_data[i])

    # Plot formatting
    plt.legend(leg_str)
    plt.title(plot_title)
    plt.xlabel(plot_xlabel)
    plt.ylabel(plot_ylabel)
    ax.set_yscale(plot_yscale)
