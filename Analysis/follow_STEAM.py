# %% File to read in the output from the STEAM code

# --------------------------- PACKAGES -----------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import os
import dataclasses
from dataclasses import dataclass
# import os.path
from os import path
import csv
import sys

c = 0
for i in range(len(sys.path)):
    if not ("Analysis/" in sys.path[i]): c = c + 1
if c == len(sys.path): sys.path.append("./Analysis/")

from functions_STEAM import *
# ------------------------------------------------------------------------------


# --------------------------- USER INPUT----------------------------------------
# Location of the output files
#  NOTE: This assumes the pythong script is being ran from the main STEAM directory
outputFileLocation = "OutputFiles"

# This is the number of time steps to pass over when pulling in data
#  NOTE: The lower the number the slower the script but the finer the resolution
#        The higher the number the faster the script but the coarser the resolution
skip_lines_num = 100

# %% Getting the data for plotting
#  NOTE: This is the most time intensive part of the script
myStr = ["rho_CR", "rho_fuel","rho_Xe"] # "nodal_density", IT ISNT WORKING FOR MATRICES YET
my_data = get_data2(outputFileLocation, skip_lines_num,myStr)

# print(my_data.data)
# # %% Plotting
#
# myStr = ["junction_velocity"]
# plot_data(my_data,myStr,"type","matrix","yscale","log")
# plt.show()
