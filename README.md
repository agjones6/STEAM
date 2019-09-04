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
