### General Use of STEAM (Scaled Time-dependent ESBWR Analysis Model)
1) Running the code (using unix commands)
  - run without compiling fortran code './runSTEAM'
  - the options listed below can be used together or seperately by putting a space
    between the options (ie: './runSTEAM -m -s')
  - run with compiling fortran code './runSTEAM -m'
    * do this if the source code is changed
    * do this if changing operating systems
  - run without executing the code './runSTEAM -s'
    * this can be used to check if there are any compiling errors without running the actual code
2) Changing the input deck (in InputFiles)
  - the file is named 'input_deck.txt'
  - the standard setup starts the reactor at full power
  - all of the input files needed are included in the InputFiles directory
  - a description of all of the options will be include in the wiki eventually
  - do not change the names in the lefthand columns. These are identifiers used by the code
3) Output
  - The outputs are currently very spread out into different files
  - The matlab code in the analysis folder can read all of the files.
  - In the future, a python script is going to be made and the output files will be consolidated.

### Analysis
#### Python File
- The python script reads all of the .csv files in the 'OutputFiles' directory
- All of the functions used in the plot_STEAM.py script are in functions_STEAM.py
- plot_STEAM.py
  * Pulling in Data
    + function: get_data(output_directory, skip_lines_num) -> This gives a structure
      with fields of 'type' (type of data), 'name' (name of the data),
      'data' (the actual numbers), and 'file' (file identifier).
    + 'output_directory' -> The location of the csv output files from running
      the STEAM code. Should be "OutputFiles"
    + 'skip_lines_num' -> This dictates how many state points are read in as data.
       For example: if the value is 3, only every 3rd state point will be read in
       as data. This is used to speed the code up. The lower the number, the higher
       the resolution of the plot but the code will run slower. The higher the number
       the lower the resolution of the plot but the code will run faster.
  * Basic Plotting
    + function: plot_data(data_structure, string_array)
    + 'plot_data' plots values based on the inputted string matrix vs elapsed time.
    + 'data_structure' -> The strucure of data from the 'get_data' function
    + 'string_array' -> An array of strings used to pull in 1 or more sets of data.
      All of the strings in the array will be plotted. The options for plotting can
      be found in data_structure.name. By default, this will be printed to the window
      to show what can be plotted.
    + example plotting 3 values on the same plot:
      plot_data(data_structure, "rho_CR", "rho_fuel", "rho_Xe"])
      plt.show()
    + NOTE: plt.show() NEEDS to be used after calling the function. The plot
      WILL NOT show unless this is used. 
  * Plot Options
    + When calling the plot function: plot_data(data,string_array,option,value)
      listed below are the options and result from using the value. All of the
      options are treated as optional input variables meaning they aren't needed
      for the plotting to work. Additionally, the user can use as many options as desired.
    + "type":
      - "matrix" -> This is for viewing nodal/junction values. The legend will
        be integers up to the number of values. ie plot_data(data, string_array, "type", "matrix")
      - else -> This is for plotting values that only have one value per time step.
        This does not require the us
    + "title":
      - string that is the title of the plot. ie plot_data(data, string_array, "title", "Temperature") will have a title of Temperature
    + "xlabel":
      - string that is the xlabel of the plot
    + "ylabel":
      - string that is the ylabel of the plot
    + "yscale":
      - "linear" -> This provides a linear y-axis
      - "log" -> This provides a logarithmic y-axis
    + DEFAULT: type = else, title = "", xlabel = "time (hr)", ylabel = "", yscale = "linear"


### Specific Case Runs
#### USING HEAT INPUT STARTING AT DESIRED CONDITIONS
1) Change StartUpCond = .TRUE.
2) Change SysAtPwr    = .TRUE.
3) Change the PwrFile to be the desired Power profile (time left column, power right column)
4) If a start location different from full power is desired, change the StartFile to the desired file that contains starting conditions

#### STARTING THE REACTOR FROM 0 POWER
1) Change StartUpCond = .TRUE.
2) Change SysAtPwr    = .False.
3) Change the StartFile = ""

#### USE STEAM DEMAND TO CONTROL REACTOR POWER
1) Change StartUpCond = .FALSE.
2) Change SysAtPwr    = .TRUE.
3) Change SteamFile to the desired profile (time left column, relative demand in the right column)
4) Set time_mult to desired multiplier (1 should be used if no time adjustments are wanted)
5) Set stFolloG to a desired level to change how fast the reactor will follow load. (~8e-3 is relatively fast where ~1e-4 is fairly slow)
6) If a start location different from full power is desired, change the StartFile to the desired file that contains starting conditions
