# SysCoX300 (STEAM)
Repository for NCSU Nuclear Engineering senior design group 007's BWRX-300 systems code.
STEAM: Scaled Time-dependent Esbwr Analysis Model

USING HEAT INPUT STARTING AT DESIRED CONDITIONS
    1) Change StartUpCond = .TRUE.
    2) Change SysAtPwr    = .TRUE.
    3) Change the PwrFile to be the desired Power profile (time left column, power right column)
    4) If a start location different from full power is desired, change the StartFile to the desired file that contains starting conditions 

STARTING THE REACTOR FROM 0 POWER
    1) Change StartUpCond = .TRUE.
    2) Change SysAtPwr    = .False.
    3) Change the StartFile = ""

USE STEAM DEMAND TO CONTROL REACTOR POWER
    1) Change StartUpCond = .FALSE.
    2) Change SysAtPwr    = .TRUE.
    3) Change SteamFile to the desired profile (time left column, relative demand in the right column)
    4) Set time_mult to desired multiplier (1 should be used if no time adjustments are wanted)
    5) Set stFolloG to a desired level to change how fast the reactor will follow load. (~8e-3 is relatively fast where ~1e-4 is fairly slow)
    6) If a start location different from full power is desired, change the StartFile to the desired file that contains starting conditions
