! Globals Module
! Includes system geometry, core distributions, etc.
MODULE globals

    USE state
    USE files
    IMPLICIT NONE

    REAL(8),PARAMETER:: pi    = acos(-1.0_8)
    REAL(8),PARAMETER:: grav  = 32.17405*3600**2 ![ft/hr^2]
    REAL(8),PARAMETER:: gc    = 32.17405*3600**2 ![lbm*ft / lbf*hr^2]
    REAL(8),PARAMETER:: ff    = 144/778.169  ![(in^2*Btu)/(lbf*ft^3)]
    REAL(8),PARAMETER:: trunc = 1e-11
    REAL(8),PARAMETER:: W2Btu = 3.41214 ![Btu/(hr*W)]

    !Input Variables
    REAL(8):: T_in      = 530, &  ![°F]
              H_G       = 1000, & ![Btu/(hr*ft^2*°F)]
              k_clad    = 9.6, &  ![Btu/(hr*ft*°F]
              lambda    = 0.90252, & !0.217922729884828 , & !Test Value pulled from HWK 7 ![ft]
              gf        = 0.974, &
              X_sepe    = 0.99, &
              Pcond     = 1., & ![psia]
              Fz        = 1.4  , &
              Fq        = 2.24 , &
              RampLimit = 100 ! Limit on the reactor ramping [%/min]

    ! Variables read in through the input files
    REAL(8):: D_o, t_clad, s, Hch, D_ch, D_vessel, n_frods, n_rods, Hsep,H_lowerPlen, &
              H_SteamDome, Hfuel, NomSteamM, nomXe, nomI, Pref, init_Press, &
              Nom_Power, SD_AREA_MULT, Init_Void_SD, CR_INS_LIMIT
    ! Variables read in through the gain file
    REAL(8):: TCV_Gain_SU, TCV_Special_SU, TCV_Gain_FP, TCV_Special_FP, FCV_Gain1_SU, &
              FCV_Gain2_SU, FCV_Gain3_SU, FCV_Gain1_FP, FCV_Gain2_FP, FCV_Gain3_FP
              ! D_o       = 0.041083333333333, & ![ft]
              ! t_clad    = 0.025/12., & ![ft]
              ! s         = 0.053333333333333, & ![ft]
              ! Hch       = 35, &    !28.5945 , & !35,&! ![ft]
              ! D_ch      = 12.6, &  !11.6652, & ![ft]
              ! D_vessel  = 13.9, &  !12.94108 ![ft]
              ! n_frods   = 22080, & ! Number of Fuel Rods
              ! n_rods    = 24000  ! Total number of Rod Locations

    !Geometry Matrices
    REAL(8):: zn(6), zj(7), Aj(12), An(11), Vn(11), Ln(11), Kn(11), De_n(11)
    INTEGER:: looplim = 5

    !-------Variable Definitions-------!
    !**INTEGER**   ==>  looplim     == Newton iteration loop limit
    !**REAL**      ==>  pi
    !                   grav        == Gravitational acceleration on Earth
    !                   gc          == Correction factor for graviational constant
    !                                  in imperial units
    !                   trunc       == Value for truncating decimals
    !                   W2Btu       == Watt to Btu multiplicative conversion factor
    !                   Hfuel       == BWRX-300 fuel height
    !                   q0          == Maximum heat flux [Btu/hr*ft^2]
    !                   T_in        == Inlet Temperature [°F]
    !                   Pref        == System Pressure [psia]
    !                   D_o         == Outer cladding diameter [ft]
    !                   s           == Pitch of the lattice
    !                   gf          == Percentage of heat deposited in the fuel
    !                   lambda      == Egtrapolated distance [ft]
    !                   G           == Initial mass flux [lbm/ft^2*hr]
    !                   X_sepe      == Quality leaving the seperator
    !                   Hch         == Height of the chimney [ft]
    !                   zn()        == Core axial node positions [ft]
    !                   zj()        == Core axial junction positions [ft]
    !                   Aj()        == Junction cross-sectional areas
    !                   An()        == Nodal cross-sectional areas
    !                   Vn()        == Nodal volumes
    !                   Ln()        == Nodal Lengths
    !                   Kn()        == Local loss coefficients
    !                   Pcond       == Pressure of the condensor
    !**CHARACTER** ==>
    !**LOGICAL**   ==>
    !----------------------------------!


CONTAINS
    ! Function for the general cosine shape
    REAL(8) FUNCTION shapez(z)
        IMPLICIT NONE

        REAL(8):: z

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==> z == Axial position [ft]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        shapez  = SIN( pi * ( (z + lambda)/(Hfuel+2*lambda))); !General sine shape function

    END FUNCTION

    !Function for the integral of the general cosine shape function
    REAL(8) FUNCTION shape_int(z) ![ft]
        IMPLICIT NONE

        REAL(8),INTENT(IN):: z

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==> z == Axial position [ft]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        shape_int = -((cos((pi*(lambda + z))/(Hfuel + 2*lambda)) - &
                 cos((lambda*pi)/(Hfuel + 2*lambda))) &
                 *(Hfuel + 2*lambda))/pi    !Shape integrated in MATLAB

    END FUNCTION

    SUBROUTINE AssignGeometry(Aj,An,Vn,zj,zn,Ln,Kn,De_n)
        IMPLICIT NONE

        INTEGER:: i, n
        REAL(8),INTENT(INOUT):: Aj(:), An(:), Vn(:), zj(:), zn(:), Ln(:), Kn(:), De_n(:)
        REAL(8):: D_LP, D_core, D_steam

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  i,n == Counter variables
        !**REAL**      ==>  zn  == [ARRAY] Core axial node positions
        !                          [ft]
        !                   zj  == [ARRAY] Core axial junction
        !                          positions [ft]
        !                   Aj  == [ARRAY] Junction cross-sectional
        !                          areas [ft^2]
        !                   An  == [ARRAY] Nodal cross-sectional
        !                          areas [ft^2]
        !                   Vn  == [ARRAY] Nodal volumes
        !`                  Ln  == [ARRAY] Length array for nodes [ft]
        !                   Hsep== Height of the separator [ft]
        !                   Kn  == Loss coefficient of each node
        !                   D_LP== Lower Plenum Diamter
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        ! bn shoulkd have 144/778 * pack 0.0410833333333
        ! DelP is in psi
        ! Values pulled from ESBWR DCD
        ! Hsep    = 26.7762  ![ft] Seperator Length
        D_LP    = D_ch !D_vessel ! Lower Plenum Diameter
        D_core  = D_ch !10.7222  ! [ft] Effective Core Diameter
        D_steam = 2.5      ! [ft] Steam line diameter

        !Junction cross-sectional areas [ft^2]
        !Aj(1)     = S**2 - (pi/4)*D_o**2 !Area of one channel, not full core

        !Declaring Length of Nodes 2-7
        Ln(2:7) = Hfuel / 6. !zj(n+1) - zj(n)
        !Declaring Length of Nodes 1,9,10,11
        Ln(1)  = H_lowerPlen        !Lower Plenum
        Ln(8)  = Hch         !Chimney
        Ln(9)  = Hsep        !Seperator
        Ln(10) = Hfuel + Hch !Downcomer
        Ln(11) = H_SteamDome !Steam Dome

        !Node cross-sectional areas [ft^2]
        An(1)  = (pi/4)*D_LP**2.0  !Lower plenum-- core area minus rod area
                           !Core area: 9.522
        DO i = 2,7
            An(i) = (pi/4)*D_core**2.0 - ((pi/4.0)*D_o**2.0)*n_frods
            ! (s**2 - (pi/4)*D_o**2) * n_frods
            ! An(i) = An(1)
        END DO

        An(8)     = pi * D_ch**2.0 / 4.0 !Chimney
        An(9)     = An(8) !Separator !Assumed half of the Chimney
        An(10)    = pi * D_vessel**2.0 / 4.0 - An(8) !Downcomer
        An(11)    = An(1)*SD_AREA_MULT       !Steam Dome, assumed same as lower plenum

        !Junction cross-sectional areas [ft^2]
        DO i = 1,11
            Aj(i) = An(i)
        END DO
        Aj(1)  = An(2) ! Outlet of the lower plenum is equal to the core inlet
        Aj(12) = pi * D_steam**2 / 4

        !Node volumes ^[ft]
        !ESBWR Volumes [ft^3]:
            !Lower Plenum       : 3567
            !Core               : 3390
            !Chimney + Separator: 9923
            !Steam Dome         : 7946
            !Downcomer          : 9040

        Vn(1)     = An(1) * Ln(1) !Lower Plenum !ASSUMED

        Vn(2)     = An(2) * Hfuel/6. !Core nodes
        DO i = 3,7
            Vn(i) = Vn(2)
        END DO

        Vn(8)     = An(8)  * Hch           !Chimney !ASSUMED
        Vn(9)     = An(9)  * Hsep          !Separator !ASSUMED
        Vn(10)    = An(10) * (Hfuel + Hch) !Downcomer
        Vn(11)    = An(11) * Ln(11)        !Steam Dome

        !Declare axial position of core junctions [ft]
        i = SIZE(zj)
        DO n = 1,i

            IF(n .EQ. 1) THEN
                zj(n) = 0
            ELSE
                zj(n) = zj(n-1) + Hfuel/6.
            END IF

        END DO

        !Declare axial position of core nodes [ft]
        i = SIZE(zj)
        DO n = 1,i

            IF(n .GT. 1) THEN
                zn(n-1) = (zj(n) + zj(n-1))/2.
            END IF

        END DO

        !Assigning Loss Coefficients
        Kn(1)  = 5.  !Inlet Plenum
        Kn(2)  = 0.25 !Core
        Kn(3)  = 0.25 !Core
        Kn(4)  = 0.25 !Core
        Kn(5)  = 0.25 !Core
        Kn(6)  = 0.25 !Core
        Kn(7)  = 0.25 !Core
        Kn(8)  = 1.5 !Chimney
        Kn(9)  = 2.  !Seperator
        Kn(10) = 10.  !Downcomer
        Kn(11) = 8 !15 !SteamDome

        !Assigning Equivalent Diameters
        De_n(1)  = 4*An(1) / (pi*D_vessel)           !Inlet plenum
        De_n(2)  = 4*(S**2-D_o**2*(pi/4)) / (pi*D_o) !Core nodes
        DO n = 3,7
            De_n(n) = De_n(2)
        END DO
        De_n(8)  = 4*An(8)  / (pi*D_ch)            !Chimney
        De_n(9)  = 4*An(9)  / (pi*D_vessel)        !Separator
        De_n(10) = 4*An(10) / (pi*(D_ch+D_vessel)) !Downcomer
        De_n(11) = 4*An(11) / (pi*D_vessel)        !Steam dome

    END SUBROUTINE


    !Function to evaluate point-wise heat flux
    REAL(8) FUNCTION q(z,q0)
        IMPLICIT NONE

        REAL(8):: z, q0

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==> z  == Axial position [ft]
        !                  q0 == Max heat flux at given time
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        q = q0 * shapez(z)

    END FUNCTION

    !Function to evaluate the (0-z) Integral of Heat Flux
    REAL(8) FUNCTION q_int(z , q0)
        IMPLICIT NONE

        REAL(8):: z, q0

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==> z  == Axial position [ft]
        !                  q0 == Max heat flux at given time
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        q_int = q0 * shape_int(z)

    END FUNCTION

    !Function evaluate point-wise quality (equilibrium model)
    REAL(8) FUNCTION X(P,alpha,G)
        IMPLICIT NONE

        !Inputs
        REAL(8):: P, alpha, G
        !Calculations
        REAL(8):: Vgj, Co

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  P      == Pressure [psia]
        !                   alpha == Vapor volume void fraction
        !                   G      == Mass flux
        !                   Vgj    == Correlation Value
        !                   Co     == Correlation value
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        Co    = 1.13
        Vgj   = 2.9*(sigma(P)*grav*gc*(rhof(P) - rhog(P)) / rhof(P)**2)**0.25

        IF((alpha .LE. 0.) .OR. (G .EQ. 0.)) THEN !Prevent dividing by 0
            X = 0
        ELSE
            X = alpha*(Co * rhog(P)/rhof(P) + (rhog(P)*Vgj)/G)/ &
                (1-alpha*Co*(1-rhog(P)/rhof(P)))
        END IF
        ! IF(X .GT. 1) THEN
        !     X = 0.97
        ! END IF
        IF(X .LT. 0) X = 0

    END FUNCTION


    !Function to evaluate point-wise void fraction (equilibrium model)
    !Uses the Zuber-Findlay and Dix correlations
    REAL(8) FUNCTION alphag(P,rho1)
        IMPLICIT NONE

        REAL(8):: P, rho1
        ! REAL(8):: C_o, Vgj, bet, b

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  P     == Pressure [psia]
        !                   rho1  == Density [lbm/ft^3]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        alphag = (rhof(P) - rho1)/(rhof(P)-rhog(P))

        IF(alphag .LE. 0.) THEN !Two phase condition
            alphag = 0
        END IF

    END FUNCTION


    !Function to evaluate point-wise bulk fluid temperature [°F]
    REAL(8) FUNCTION bulktemp(P,u)
        IMPLICIT NONE

        REAL(8):: P, u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  P  == Pressure [psia]
        !                   u  == Specific Internal Energy [Btu/lbm]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        bulktemp = Temp(u)

        !Bulk fluid temperature cannot be greater than saturation
        IF(bulktemp .GT. Tsat(P)) THEN
            bulktemp = Tsat(P)
        END IF

    END FUNCTION


    !Subroutine to test convergence
    SUBROUTINE converge(v_k,n1,v_k1,n2,tolerance,crit)
        IMPLICIT NONE

        INTEGER,INTENT(IN):: n1, n2
        REAL(8),INTENT(IN):: v_k(n1,2), v_k1(n2,1), tolerance
        LOGICAL,INTENT(OUT):: crit

        INTEGER:: n
        REAL(8):: test

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n1,n2       == Matrix dimension sizes
        !                   n           == Counter variable
        !**REAL**      ==>  v_k()       == Velocity (k) matrix
        !                   v_k1()      == Velocity (k+1) matrix
        !                   tolerance   == Convergence tolerance limit
        !                   test        == Convergence test
        !**CHARACTER** ==>
        !**LOGICAL**   ==>  crit        == Convergence boolean
        !----------------------------------!

        DO n = 1,SIZE(v_k,1)

            test = ABS( v_k(n,1) - v_k1(n,1) )
            IF( test .GT. tolerance ) THEN
                crit = .FALSE.
                EXIT
            ELSE IF((n .EQ. SIZE(v_k,1)) .AND. (test .LE. tolerance)) THEN
                crit = .TRUE.
            END IF

        END DO

    END SUBROUTINE

    !Finds the maximum of a set of REAL(8) values
    REAL(8) FUNCTION Maximum(vals)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: vals(:,:)

        INTEGER             :: i

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  size    == Number of values
        !                   i       == Counter variable
        !**REAL**      ==>  vals    == Values to be compared
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        Maximum = 0
        DO i = 1, SIZE(vals,1)
            IF( vals(i,1) .GT. Maximum ) Maximum = vals(i,1)
        END DO

    END FUNCTION

    SUBROUTINE ReadGeometry(FILENAME)
        ! This subroutine will read in the geometry of the system
        IMPLICIT NONE

        ! Inputs
        CHARACTER(LEN=*),INTENT(IN):: FILENAME

        ! Local
        REAL(8):: dum2
        CHARACTER(LEN=50):: dum1
        INTEGER:: i, io, num_lines, UNIT

        !-------Variable Definitions-------!
        !**INTEGER**   ==> UNIT      == The unit of the file to be opened
        !                  i         == Counting for the loop
        !                  io        == for error handling
        !                  num_lines == number of lines in the file
        !**REAL**      ==> dum2      == The current value of the line
        !**CHARACTER** ==> FILENAME  == The name of the input file
        !                  dum1      == Just a space for handling inputs
        !**LOGICAL**   ==>
        !----------------------------------!

        ! Opening File
        UNIT = 43
        OPEN(UNIT=UNIT,FILE=FILENAME,ACTION="READ")

        ! Counting the number of lines in the file
        num_lines = 0
        DO
            READ(UNIT,*,iostat=io) dum1
            IF (io/=0) EXIT
            IF(dum1 .NE. "") THEN
                num_lines = num_lines + 1
            END IF
        END DO

        ! Reading in each line of the file
        REWIND(UNIT = UNIT)
        DO i = 1,num_lines-1
            ! Skipping the first line in the file
            IF(i .EQ. 1) THEN
                READ(UNIT,*)
            END IF

            ! Reading in the string and the value next to the string
            READ(UNIT,*,iostat=io) dum1, dum2

            !--> Assigning all of the variables to the appropriate value
            IF(TRIM(dum1) .EQ. "Clad_Diameter")       D_o          = dum2
            IF(TRIM(dum1) .EQ. "Clad_Thickness")      t_clad       = dum2
            IF(TRIM(dum1) .EQ. "Rod_Pitch")           s            = dum2
            IF(TRIM(dum1) .EQ. "Chimney_Height")      Hch          = dum2
            IF(TRIM(dum1) .EQ. "Chimney_Diameter")    D_ch         = dum2
            IF(TRIM(dum1) .EQ. "Vessel_Diameter")     D_vessel     = dum2
            IF(TRIM(dum1) .EQ. "Seperator_Length")    Hsep         = dum2
            IF(TRIM(dum1) .EQ. "Number_Fuel_Rods")    n_frods      = dum2
            IF(TRIM(dum1) .EQ. "Number_Total_Rods")   n_rods       = dum2
            IF(TRIM(dum1) .EQ. "Lower_Plenum_Height") H_lowerPlen  = dum2
            IF(TRIM(dum1) .EQ. "SteamDome_Height")    H_SteamDome  = dum2
            IF(TRIM(dum1) .EQ. "Fuel_Height")         Hfuel        = dum2
            IF(TRIM(dum1) .EQ. "Nom_Xenon")           nomXe        = dum2
            IF(TRIM(dum1) .EQ. "Nom_Iodine")          nomI         = dum2
            IF(TRIM(dum1) .EQ. "Nom_mSteam")          NomSteamM    = dum2
            IF(TRIM(dum1) .EQ. "ref_Press")           Pref         = dum2
            IF(TRIM(dum1) .EQ. "initial_Press")       init_Press   = dum2
            IF(TRIM(dum1) .EQ. "Nom_Rx_Power")        Nom_Power    = dum2
            IF(TRIM(dum1) .EQ. "SD_Area_mult")        SD_AREA_MULT = dum2
        END DO

        ! Closing the file
        CLOSE(UNIT = UNIT, STATUS = 'KEEP')
    END SUBROUTINE

    SUBROUTINE ReadGains(FILENAME)
        ! This subroutine will read in the geometry of the system
        IMPLICIT NONE

        ! Inputs
        CHARACTER(LEN=*),INTENT(IN):: FILENAME

        ! Local
        REAL(8):: dum2
        CHARACTER(LEN=50):: dum1
        INTEGER:: i, io, num_lines, UNIT

        !-------Variable Definitions-------!
        !**INTEGER**   ==> UNIT      == The unit of the file to be opened
        !                  i         == Counting for the loop
        !                  io        == for error handling
        !                  num_lines == number of lines in the file
        !**REAL**      ==> dum2      == The current value of the line
        !**CHARACTER** ==> FILENAME  == The name of the input file
        !                  dum1      == Just a space for handling inputs
        !**LOGICAL**   ==>
        !----------------------------------!

        ! Opening File
        UNIT = 43
        OPEN(UNIT=UNIT,FILE=FILENAME,ACTION="READ")

        ! Counting the number of lines in the file
        num_lines = 0
        DO
            READ(UNIT,*,iostat=io) dum1
            IF (io/=0) EXIT
            IF(dum1 .NE. "") THEN
                num_lines = num_lines + 1
            END IF
        END DO

        ! Reading in each line of the file
        REWIND(UNIT = UNIT)
        DO i = 1,num_lines-1
            ! Skipping the first line in the file
            IF(i .EQ. 1) THEN
                READ(UNIT,*)
            END IF

            ! Reading in the string and the value next to the string
            READ(UNIT,*,iostat=io) dum1, dum2

            !--> Assigning all of the variables to the appropriate value
            IF(TRIM(dum1) .EQ. "TCV_Gain_SU")    TCV_GAIN_SU     = dum2
            IF(TRIM(dum1) .EQ. "TCV_Gain_FP")    TCV_GAIN_FP     = dum2
            IF(TRIM(dum1) .EQ. "TCV_Special_SU") TCV_SPECIAL_SU  = dum2
            IF(TRIM(dum1) .EQ. "TCV_Special_FP") TCV_SPECIAL_FP  = dum2

            IF(TRIM(dum1) .EQ. "FCV_Gain1_SU")   FCV_GAIN1_SU    = dum2
            IF(TRIM(dum1) .EQ. "FCV_Gain1_FP")   FCV_GAIN1_FP    = dum2

            IF(TRIM(dum1) .EQ. "FCV_Gain2_SU")   FCV_GAIN2_SU    = dum2
            IF(TRIM(dum1) .EQ. "FCV_Gain2_FP")   FCV_GAIN2_FP    = dum2

            IF(TRIM(dum1) .EQ. "FCV_Gain3_SU")   FCV_GAIN3_SU    = dum2
            IF(TRIM(dum1) .EQ. "FCV_Gain3_FP")   FCV_GAIN3_FP    = dum2
        END DO

        ! Closing the file
        CLOSE(UNIT = UNIT, STATUS = 'KEEP')
    END SUBROUTINE

    SUBROUTINE ReadCond(FILENAME, total_time, psuedoTime, &
                        StartUpCond, switch_Prof_time, SysAtPwr,switch_Gain_time,&
                        StartFile, TakeSnap, snapInt, SnapFile, SteamFile, stFollowG,&
                        readPoisonFile, time_mult, cont_meth, G_feed, &
                        PoisonEquil, writePoisonFile, PwrFile, PumpFile, CPRtol, &
                        tStep_tol_SU,tStep_tol_FP)
        ! This subroutine will read in the geometry of the system
        IMPLICIT NONE

        ! --> Inputs
        CHARACTER(LEN=*),INTENT(IN):: FILENAME

        ! --> Outputs
        CHARACTER(LEN=*),INTENT(INOUT):: StartFile, &
                                          SnapFile, SteamFile, readPoisonFile, &
                                          writePoisonFile, PwrFile, PumpFile

        REAL(8),INTENT(INOUT):: total_time, psuedoTime, snapInt, CPRtol,stFollowG, &
                                G_feed, switch_Prof_time, switch_Gain_time, &
                                tStep_tol_SU, tStep_tol_FP
        REAL(4),INTENT(INOUT):: time_mult

        LOGICAL,INTENT(INOUT):: StartUpCond, SysAtPwr, TakeSnap, PoisonEquil

        INTEGER,INTENT(INOUT):: cont_meth

        ! Local
        REAL(8):: dum_real
        CHARACTER(LEN=75):: dum1, dum_str, test
        INTEGER:: i, io, num_lines, UNIT
        LOGICAL:: Test_Log

        !-------Variable Definitions-------!
        !**INTEGER**   ==> UNIT      == The unit of the file to be opened
        !                  i         == Counting for the loop
        !                  io        == for error handling
        !                  num_lines == number of lines in the file
        !**REAL**      ==> dum2      == The current value of the line
        !**CHARACTER** ==> FILENAME  == The name of the input file
        !                  dum1      == Just a space for handling inputs
        !**LOGICAL**   ==>
        !----------------------------------!

        ! Opening File
        UNIT = 43
        OPEN(UNIT=UNIT,FILE=FILENAME,ACTION="READ")

        ! Counting the number of lines in the file
        num_lines = 0
        DO
            READ(UNIT,*,iostat=io) dum1
            IF (io/=0) EXIT
            IF(dum1 .NE. "") THEN
                num_lines = num_lines + 1
            END IF
        END DO

        ! Reading in each line of the file
        REWIND(UNIT = UNIT)
        DO i = 1,num_lines
            ! Reading in the line identifier
            READ(UNIT,*,iostat=io) dum1

            !--> Assigning all of the variables to the appropriate value

            ! ---------------------- CONDITIONS --------------------------------
            IF(TRIM(dum1) .EQ. "total_time") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, total_time
            END IF
            IF(TRIM(dum1) .EQ. "start_time") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, psuedoTime
            END IF
            IF(TRIM(dum1) .EQ. "Use_Power_Prof") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, StartUpCond
            END IF
            IF(TRIM(dum1) .EQ. "Switch_Profile_Time") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, switch_Prof_time
            END IF
            IF(TRIM(dum1) .EQ. "Full_Power_Gains") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, SysAtPwr
            END IF
            IF(TRIM(dum1) .EQ. "Switch_Gains_Time") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, switch_Gain_time
            END IF
            IF(TRIM(dum1) .EQ. "Start_File") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, StartFile
            END IF
            IF(TRIM(dum1) .EQ. "Take_Snapshot") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, TakeSnap
            END IF
            IF(TRIM(dum1) .EQ. "Snapshot_Interval") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, snapInt
            END IF
            IF(TRIM(dum1) .EQ. "Snapshot_Write_File") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, SnapFile
            END IF
            IF(TRIM(dum1) .EQ. "Steam_Demand_File") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, SteamFile
            END IF
            IF(TRIM(dum1) .EQ. "Steam_Profile_mult") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, time_mult
            END IF
            IF(TRIM(dum1) .EQ. "Control_Method") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, cont_meth
            END IF
            IF(TRIM(dum1) .EQ. "Xe_and_Sm_Equilibrium") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, PoisonEquil
            END IF
            IF(TRIM(dum1) .EQ. "Read_Poison_File") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, readPoisonFile
            END IF
            IF(TRIM(dum1) .EQ. "Write_Poison_File") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, writePoisonFile
            END IF
            IF(TRIM(dum1) .EQ. "Power_Profile") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, PwrFile
            END IF
            IF(TRIM(dum1) .EQ. "Pump_Profile") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, PumpFile
            END IF
            IF(TRIM(dum1) .EQ. "CPR_Tolerance") THEN
                BACKSPACE(UNIT = UNIT) ! Going to the line if there is a match
                READ(UNIT,*) dum1, CPRtol
            END IF

            ! ---------------------- GEOMETRY ----------------------------------
            IF(TRIM(dum1) .EQ. "Clad_Diameter") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, D_o
            END IF
            IF(TRIM(dum1) .EQ. "Clad_Thickness") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, t_clad
            END IF
            IF(TRIM(dum1) .EQ. "Rod_Pitch") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, s
            END IF
            IF(TRIM(dum1) .EQ. "Chimney_Height") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, Hch
            END IF
            IF(TRIM(dum1) .EQ. "Chimney_Diameter") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, D_ch
            END IF
            IF(TRIM(dum1) .EQ. "Vessel_Diameter") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, D_vessel
            END IF
            IF(TRIM(dum1) .EQ. "Seperator_Length") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, Hsep
            END IF
            IF(TRIM(dum1) .EQ. "Number_Fuel_Rods") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, n_frods
            END IF
            IF(TRIM(dum1) .EQ. "Number_Total_Rods") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, n_rods
            END IF
            IF(TRIM(dum1) .EQ. "Lower_Plenum_Height") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, H_lowerPlen
            END IF
            IF(TRIM(dum1) .EQ. "SteamDome_Height") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, H_SteamDome
            END IF
            IF(TRIM(dum1) .EQ. "Fuel_Height") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, Hfuel
            END IF
            IF(TRIM(dum1) .EQ. "Nom_Xenon") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, nomXe
            END IF
            IF(TRIM(dum1) .EQ. "Nom_Iodine") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, nomI
            END IF
            IF(TRIM(dum1) .EQ. "Nom_mSteam") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, NomSteamM
            END IF
            IF(TRIM(dum1) .EQ. "ref_Press") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, Pref
            END IF
            IF(TRIM(dum1) .EQ. "initial_Press") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, init_Press
            END IF
            IF(TRIM(dum1) .EQ. "Nom_Rx_Power") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, Nom_Power
            END IF
            IF(TRIM(dum1) .EQ. "SD_Area_mult") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, SD_AREA_MULT
            END IF
            IF(TRIM(dum1) .EQ. "Initial_SD_Void") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, Init_Void_SD
            END IF

            ! ----------------------- GAINS ------------------------------------
            IF(TRIM(dum1) .EQ. "TCV_Gain_SU") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, TCV_GAIN_SU
            END IF
            IF(TRIM(dum1) .EQ. "TCV_Gain_FP") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, TCV_GAIN_FP
            END IF
            IF(TRIM(dum1) .EQ. "TCV_Special_SU") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, TCV_SPECIAL_SU
            END IF
            IF(TRIM(dum1) .EQ. "TCV_Special_FP") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, TCV_SPECIAL_FP
            END IF

            IF(TRIM(dum1) .EQ. "FCV_Gain1_SU") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, FCV_GAIN1_SU
            END IF
            IF(TRIM(dum1) .EQ. "FCV_Gain1_FP") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, FCV_GAIN1_FP
            END IF

            IF(TRIM(dum1) .EQ. "FCV_Gain2_SU") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, FCV_GAIN2_SU
            END IF
            IF(TRIM(dum1) .EQ. "FCV_Gain2_FP") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, FCV_GAIN2_FP
            END IF

            IF(TRIM(dum1) .EQ. "FCV_Gain3_SU") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, FCV_GAIN3_SU
            END IF
            IF(TRIM(dum1) .EQ. "FCV_Gain3_FP") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, FCV_GAIN3_FP
            END IF

            IF(TRIM(dum1) .EQ. "Steam_Follow") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, stFollowG
            END IF
            IF(TRIM(dum1) .EQ. "Feedwater_Follow") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, G_feed
            END IF
            IF(TRIM(dum1) .EQ. "Time_Step_tol_SU") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, tStep_tol_SU
            END IF
            IF(TRIM(dum1) .EQ. "Time_Step_tol_FP") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, tStep_tol_FP
            END IF
            IF(TRIM(dum1) .EQ. "Control_Rod_Speed") THEN
                BACKSPACE(UNIT = UNIT)
                READ(UNIT,*) dum1, CR_INS_LIMIT
            END IF

            IF(TRIM(dum1) .EQ. "END_FILE") RETURN
        END DO

        ! Closing the file
        CLOSE(UNIT = UNIT, STATUS = 'KEEP')

    END SUBROUTINE

    ! !Numerical Integrator using the Trapezoid Rule
    ! REAL(8) FUNCTION NumInt(eqbins)
    !     IMPLICIT NONE
    !
    !     REAL(8),INTENT(IN)  :: eq(:)
    !     INTEGER,INTENT(IN)  :: bins
    !     REAL(8)             :: min, max,interval, val(bins-1)
    !     INTEGER             :: n
    !
    !     !-------Variable Definitions-------!
    !     !**INTEGER**   ==>  n        == Counter variable
    !     !**REAL**      ==>  eq      == Function to be integrated
    !     !                   min      == Lower bound of integration
    !     !                   max      == Upper bound of integration
    !     !                   bins     == Number of bins (AKA mesh/grid size)
    !     !                   interval == Size of each bin
    !     !                   val()    == Value of integration value at end of each bin
    !     !**CHARACTER** ==>
    !     !**LOGICAL**   ==>
    !     !----------------------------------!
    !
    !     min = 1
    !     max = SIZE(eq)
    !
    !     interval = (max - min) / bins
    !     DO n = 1,bins-1
    !         val(n) = (n * interval) + min
    !     END DO
    !
    !     !Integration
    !     NumInt = 0
    !     DO n = 1,bins
    !         IF(n .EQ. 1) THEN
    !             NumInt = NumInt + (eq(min) + eq(val(n))) / (2*interval)
    !         ELSE IF((n .GT. 1) .AND. (n .LT. bins)) THEN
    !             NumInt = NumInt + (eq(val(n-1)) + eq(val(n))) / (2*interval)
    !         ELSE
    !             NumInt = NumInt + (eq(val(n-1)) + eq(max)) / (2*interval)
    !         END IF
    !     END DO
    !
    ! END FUNCTION


END MODULE
