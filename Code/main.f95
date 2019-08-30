! This is the main body of the GE-Hitachi BWRX-300 systems code
! Created by NCSU Senior Design Team 007, 2018-2019
! Coded by Benjamin Laramee and Andy Jones

!! This program assumes an equilibrium model where, when two phases are...
!! present, both are saturated

PROGRAM Systems_Code_X300

    USE state
    USE globals
    USE controllers
    USE conditional
    USE files
    USE steamdome
    USE pressuredrop
    USE matrices
    USE safety
    USE heat
    IMPLICIT NONE

    !===============================================================================
    !                   VARIABLE DECLARATIONS and DEFINITIONS
    !===============================================================================

    INTEGER:: IPIV(1,13), stat, i, j, n, o, loopde, Err
    !Loop/Convergence Variables
    INTEGER:: tloop = 0, kloop = 0, count = 0, loop = 0
    LOGICAL:: conv = .FALSE., dt_change = .FALSE.
    REAL(8):: tol = 0.001
    !Thermodynamic Matrices
    REAL(8):: Tn(11,2), Tj(12,2), agn(11,2), agj(12,2), rhon(11,2), &
    rhoj(12,2), uln(11,2), ulj(12,2), rho_un(11,2), rho_uj(12,2),   &
    Xn(11,2), Xj(12,2), Gj(12,2), Xi(11), v_j(12,2), vg_j(12,2)
    !Equation Matrices
    REAL(8):: Amat(13,13), Bmat(13,1), alpha(11,11), beta(11), gam(11), &
    psi(3), a(11), aprev(11), b(11), bprev(11), M(6), E(6), LM(6),      &
    SM(11), SE(11), qn(11,1), dPsi_n(11)

    !Allocatable Arrays
    REAL(8),ALLOCATABLE:: Power(:), t_q(:), Pump(:), t_p(:), Steam(:), t_st(:)

    !Other
    REAL(8):: dt, C12, dP, y, w, r, y2, w2, level, mfeed, msteam,      &
    mfeednew, hFD, S_Ml11, Ksys, Ktcv(2), Xi_cond, SME(2,2), vg1, vg2, &
    vl1, vl2, Psys(1,2), tcvPos, Aj12, dqdt, total_time, qi, qf, d_ul11,   &
    al_rhol11(1,2), cPower, pPower, cPump, cSteamD, snapCur = 0, snapInt, tcvVel, Pprev,&
    Vgj, dPrev, TFD, CPR, Gain, Gainf, psuedoTime, CPRtol, G1, G2, G3, stFollowG, &
    cC, react, react_CR, avg_void, rho_ex, rho_fuel, rho_void, G_feed, Q0_Rx, Q0_f, &
    N_I, N_Xe, RelXe, avg_void1, avg_void2, avg_void3, avg_void4,&
    avg_relPwr, avg_relPwr1, avg_relPwr2, avg_relPwr3, avg_relPwr4, &
    switch_Prof_time, switch_Gain_time, tStep_tol_SU, tStep_tol_FP, Linear_Heat_Rate

    LOGICAL:: WantDebug = .FALSE., StartUpCond = .TRUE., SysAtPwr = .FALSE., &
    TakeSnap = .TRUE., PoisonEquil = .TRUE., has_Prof_Switched = .FALSE., has_Gain_Switched = .FALSE.

    CHARACTER(LEN=75):: StartFile, SnapFile, SteamFile, PwrFile, PumpFile, &
                          readPoisonFile, writePoisonFile, INPUT_FILE
    REAL(4):: time_mult

    !Time Control Variables
    REAL(8):: dtmin, dtmax, el_time
    !Temporary
    REAL(8):: a12, a11, a9, b12, b11, b9, al9, DP1(7), K, Vjn, &
                test(13,1), dum1(11,1), q0hot
    INTEGER:: testcounter, tcvCount, cont_meth


    !-------Variable Definitions-------!
    !**INTEGER**   ==>  IPIV()    ==  Pivot array, necessary to use DGESV
    !                   stat      ==  Info (error) variable
    !                   i,j,n,o   ==  Counter/placeholder variables
    !                   tloop     ==  Time loop counter
    !                   kloop     ==  Convergence loop counter
    !                   count     ==  Convergence counter
    !**REAL**      ==>  tol       ==  Convergence tolerance
    !                   Tn()      ==  Temperatures:      node     [°F]
    !                   Tj()      ==  Temperatures:      junction [°F]
    !                   agn()     ==  Void Fractions:    node     []
    !                   agj()     ==  Void Fractions:    junction []
    !                   rhon()    ==  Densities:         node     [lbm/ft^3]
    !                   rhoj()    ==  Densities:         junction [lbm/ft^3]
    !                   uln()     ==  Internal Energies: node     [Btu/lbm]
    !                   ulj()     ==  Internal Energies: junction [Btu/lbm]
    !                   rho_un()  ==  rho*u:             node     [Btu/ft^3]
    !                   rho_uj()  ==  rho*u:             junction [Btu/ft^3]
    !                   Xj()      ==  Quality:           junction
    !                   Gj()      ==  Mass Flux:         junction [lbm/ft^2*hr]
    !                   Xi()      ==  Xi matrix value
    !                   Xn()      ==  Quality:           node
    !                   v_j()     ==  velocity (t,k)   node     [ft/hr]
    !                   qn()      ==  Heat Generation rate: node  [btu/hr]
    !                   Amat()    ==  A matrix
    !                   Bmat()    ==  B matrix, after DGESV is solution matrix
    !                   alpha()   ==  Alpha values
    !                   beta()    ==  Beta values
    !                   gam()     ==  Gamma values
    !                   psi()     ==  Psi values
    !                   a()       ==  a(n) values
    !                   aprev()   ==  a(n-1) values
    !                   b()       ==  b(n) values
    !                   bprev()   ==  b(n-1) values
    !                   M(),E()   ==  Matrices for Steam dome
    !                   LM        ==  Matrices for Steam dome
    !                   SM()      ==  Mass conservation for Xi values
    !                   SE()      ==  Energy conservation for Xi values
    !                   dPsi_n()  ==  Change in internal Energy or void fraction
    !                   dt        ==  Current time step
    !                   C12       ==  C value for junction 12
    !                   dP        ==  delta_P change in pressure
    !                   y,w,r     ==  Subroutine output variables
    !                   y2, w2    ==  Subroutine output placeholder
    !                   DP1       ==  Pressure drop testing variable
    !                   K         ==  Local Loss Coefficient
    !                   level     ==  Level of water in the downcomer [ft]
    !                   mfeed     ==  Mass flow rate of the feed water [lbm/hr]
    !                   msteam    ==  Mass flow rate of steam [lbm/hr]
    !                   hFD       ==  Feed water enthalpy [Btu/lbm]
    !                   Ksys      ==  Total loss coefficient around the loop
    !                   Ktcv()    ==  Loss coefficient for the TCV (old,new)
    !                   Xi_cond   ==  Xi in the condensor
    !                   SME()     ==  Contains SME1 & SME2 (1|t,k ; 2|t,k)
    !                   Psys()    ==  System pressure (t,k) [psia]
    !                   tcvPos    ==  Current Turbine Control Valve Position
    !                   Aj12      ==  Constant cross section area of the turbine control valve
    !                   dqdt      ==  Change in heat over time
    !                   total_time    ==  Total Change in time
    !                   qi        ==  initial q0
    !                   qf        ==  Final q0
    !                   el_time   ==  Elapsed Time
    !                   dtmin     ==  Minimum time step size
    !                   dtmax     ==  Maximum time step size
    !                   d_ul11    ==  Change in specific internal energy for the steam dome
    !                   Pump()    ==  Allocatable Array for the change in pressure given by the pump [psia]
    !                   t_p()     == Time when the pump changes power levels
    !                   cPump     == Current Pump output
    !                   tcvVel    == The instantaneous change in valve position over time
    !                   al_rhol11 == void * density in the steam dome (t,k+1)
    !                   Power()   == Allocatable array for the Power Profile
    !                   t_q()     == Allocatable array for the time of the power manuever
    !                   cPower    == Current Power Level
    !                   pPower    == Previous Power Level
    !                   Steam()   == Allocatable array for steam demand profile
    !                   t_st()    == Allocatable arrau for the time of the steam demand
    !                   cSteam    == Current Steam Demand
    !                   Pprev     == Previous time step pressure
    !                   TFD       == Feed water temperature
    !                   q0hot     == q0 (heat flux) in the hot channel
    !**CHARACTER** ==>
    !**LOGICAL**   ==>  conv      ==  Convergence boolean
    !----------------------------------!

    !---> Debugging Files
    ! CALL openfile(11,'Output.csv')
    ! CALL openfile(13,'Vel.csv')
    ! CALL openfile(14,'Vel2.csv')
    ! CALL openfile(15,'VelMat.csv')

    CALL openfile(16,'../OutputFiles/velocityMat.csv') ! /SysCoX300_v2

    CALL openfile(17,'../OutputFiles/uMat.csv')
    CALL openfile(18,'../OutputFiles/densMat.csv')
    CALL openFile(19,'../OutputFiles/voidMat.csv')
    CALL openfile(20,'../OutputFiles/convMat.csv')
    CALL openfile(21,'../OutputFiles/relVMat.csv')
    CALL openfile(22,'../OutputFiles/reactivity.csv')
    CALL openfile(23,'../OutputFiles/random.csv')



    ! CALL openFile(21,'TestyTest.csv')


    ! CALL openFile(31,'IntervalSnap.csv')

    ! Xenon and Iodine File
    ! CALL openfile(62,'myPoison.txt')

    !CALL testsub()

    !===============================================================================
    !                   INITIALIZING THE SYSTEM and INPUTS
    !===============================================================================
    !----------------------> Initializing Values <------------------------------
    v_j     = 0.                  ![ft/hr]]
    el_time = 0.                  ![hr]
    dtmin   = 1e-4                ![hr]
    dt      = 1e-2                ![hr]
    dtmax   = 0.1                 ![hr]
    Q0_Rx   = 0.00                ![Btu/lbm*ft^2]

    !-------------------------------> INPUTS <----------------------------------
    INPUT_FILE = "../InputFiles/input_deck.txt" ! <--------------------- THE ONLY THING THAT NEEDS TO BE ALTERED PER RUN
    ! BWRX_10.txt is the newest file

    CALL ReadCond(INPUT_FILE, total_time, psuedoTime, &
                  StartUpCond, switch_Prof_time, SysAtPwr,switch_Gain_time,&
                  StartFile, TakeSnap, snapInt, SnapFile, SteamFile, stFollowG,&
                  readPoisonFile, time_mult, cont_meth, G_feed, &
                  PoisonEquil, writePoisonFile, PwrFile, PumpFile, CPRtol, &
                  tStep_tol_SU, tStep_tol_FP)
    !-------------------------> USING INPUTS <----------------------------------
    !Assign Geometry
    CALL AssignGeometry(Aj,An,Vn,zj,zn,Ln,Kn,De_n)

    ! Handling the switch of profiles and gains
    IF(ABS(switch_Prof_time) .LE. 1e-4) has_Prof_Switched = .TRUE.
    IF(ABS(switch_Gain_time) .LE. 1e-4) has_Gain_Switched = .TRUE.

    !--> Initializing Steam Demand
    OPEN(UNIT = 41, FILE = SteamFile, ACTION = 'READ')
    CALL ReadProf(MATRIX = Steam, TIME = t_st, FILE = 41, TIME_MULT = time_mult)
    CLOSE(UNIT = 41, STATUS = 'KEEP')

    !--> Initializing Power Profile (only used when StartUpCond = TRUE)
    OPEN(UNIT = 41, FILE = PwrFile, ACTION = 'READ')
    CALL ReadProf(MATRIX = Power, TIME = t_q, FILE = 41)
    CLOSE(UNIT = 41, STATUS = 'KEEP')

    !--> Initializing Pump Profile
    OPEN(UNIT = 41, FILE = PumpFile, ACTION = 'READ')
    CALL ReadProf(MATRIX = Pump, TIME = t_p, FILE = 41) !TIME_MULT is a REAL(4) input
    CLOSE(UNIT = 41, STATUS = 'KEEP')
    ! Pump is always 0 if there is no start up
    IF(StartUpCond .EQV. .FALSE.) Pump = 0

    !--> Initializing Xenon and Samarium
    IF(PoisonEquil) THEN

    ELSE
        OPEN(UNIT = 62, FILE = readPoisonFile, ACTION = 'READ')
        READ(62,*) N_Xe, N_I
        CLOSE(UNIT = 62, STATUS = 'KEEP')
    END IF

    !---> INITIAL STATE VALUES
    N_Xe    = nomXe
    N_I     = nomI
    RelXe   = N_Xe/nomXe

    ! Initial System Pressure
    Psys  = init_Press             ![psia]
    dPrev = 0.001

    ! Counting Variable for determining the tcv change
    tcvCount = 0

    ! Initial Turbine Control Valve values
    tcvPos  = 0.     !Full-closed
    Aj12    = Aj(12) !Original cross-sectional area of the steam line
    Ktcv(2) = 1e6

    ! Initial Feed Flow Mass Flow Rate and Feed Enthalpy
    mfeed  = 0.      !Testing value pulled from project
    TFD    = 420     ![degF] Feed Temperature

    ! Iniitializing the State Point Variables with a snapshot file or from 0 conditions
    ! Opens Desired Start File
    IF((StartFile .EQ. "")) THEN
        CALL InitialProps(v_j,T_in,Psys,al_rhol11,Tn,Tj,agn,agj,rhon,rhoj,uln,ulj, &
        rho_un,rho_uj,Xn,Xj,Gj,Q0_Rx,dt,el_time,tcvPos,mfeed)
        msteam  = v_j(12,1) * rhoj(12,1) * Aj(12)
        cPower  = (msteam/NomSteamM)*Nom_Power
    ELSE
        OPEN(UNIT = 31, FILE = StartFile, ACTION = 'READ')
        CALL InitialProps(v_j,T_in,Psys,al_rhol11,Tn,Tj,agn,agj,rhon,rhoj,uln,ulj, &
        rho_un,rho_uj,Xn,Xj,Gj,Q0_Rx,dt,el_time,tcvPos,mfeed,31)
        CLOSE(UNIT = 31, STATUS = 'KEEP')
        Pprev = Psys(1,1)
        msteam  = v_j(12,1) * rhoj(12,1) * Aj(12)
        cPower  = (msteam/NomSteamM)*Nom_Power

        !---> INPUT starting time of simulation to match profiles
        el_time = psuedoTime
    END IF
    Q0_f = Q0_Rx


    !===============================================================================
    !                         STARTING OUTER LOOP
    !===============================================================================

    DO WHILE(el_time .LE. total_time) !Start outer loop (t)
        !---> This changes the tcv gains depending on the power level of the reactor
        IF((el_time .GT. switch_Prof_time) .AND. (has_Prof_Switched .EQV. .FALSE.)) THEN
            IF(StartUpCond) THEN
                StartUpCond = .FALSE.
            ELSE
                StartUpCond = .TRUE.
            END IF
            has_Prof_Switched = .TRUE.
            tloop = 0
        END IF
        IF((el_time .GT. switch_Gain_time) .AND. (has_Gain_Switched .EQV. .FALSE.)) THEN
            IF(SysAtPwr) THEN
                SysAtPwr = .FALSE.
            ELSE
                SysAtPwr = .TRUE.
            END IF
            has_Gain_Switched = .TRUE.
        END IF

        msteam = v_j(12,1) * rhoj(12,1) * Aj(12) ![lbm/hr]

        !---> Calculating the average over time for void and relative power
        IF(tloop .EQ. 0) THEN
            ! Void
            avg_void = 0
            DO n = 2,7
                avg_void = agn(n,1) + avg_void
            END DO
            avg_void = avg_void/6.0

            ! Power
            avg_relPwr = cPower/Nom_Power
        ELSEIF(tloop .EQ. 1) THEN
            ! Void
            avg_void1 = avg_void
            avg_void  = 0
            DO n = 2,7
                avg_void = agn(n,1) + avg_void
            END DO
            avg_void = avg_void/6.0
            avg_void = (avg_void + avg_void1)/2.0

            ! Power
            avg_relPwr1 = avg_relPwr
            avg_relPwr  = cPower/Nom_Power
            avg_relPwr  = (avg_relPwr + avg_relPwr1)/2.0
        ELSEIF(tloop .EQ. 2) THEN
            ! Void
            avg_void2 = avg_void
            avg_void  = 0
            DO n = 2,7
                avg_void = agn(n,1) + avg_void
            END DO
            avg_void = avg_void/6.0
            avg_void = (avg_void + avg_void1 + avg_void2)/3.0

            ! Power
            avg_relPwr2 = avg_relPwr
            avg_relPwr  = cPower/Nom_Power
            avg_relPwr  = (avg_relPwr + avg_relPwr1 + avg_relPwr2)/3.0
        ELSEIF(tloop .EQ. 3) THEN
            ! Void
            avg_void3 = avg_void
            avg_void  = 0
            DO n = 2,7
                avg_void = agn(n,1) + avg_void
            END DO
            avg_void = avg_void/6.0
            avg_void = (avg_void + avg_void1 + avg_void2 + avg_void3)/4.0

            ! Power
            avg_relPwr3 = avg_relPwr
            avg_relPwr  = cPower/Nom_Power
            avg_relPwr  = (avg_relPwr + avg_relPwr1 + avg_relPwr2 + avg_relPwr3)/4.0
        ELSE
            ! Void
            avg_void1 = avg_void
            avg_void2 = avg_void1
            avg_void3 = avg_void2
            avg_void  = 0
            DO n = 2,7
                avg_void = agn(n,1) + avg_void
            END DO
            avg_void = avg_void/6.0
            avg_void = (avg_void + avg_void1 + avg_void2 + avg_void3)/4.0

            ! Power
            avg_relPwr1 = avg_relPwr
            avg_relPwr2 = avg_relPwr1
            avg_relPwr3 = avg_relPwr2
            avg_relPwr  = cPower/Nom_Power
            avg_relPwr  = (avg_relPwr + avg_relPwr1 + avg_relPwr2 + avg_relPwr3)/4.0
        END IF

        !Writing data at the desired Interval
        IF((el_time - snapCur .GT. snapInt) .AND. TakeSnap) THEN
            snapCur = snapCur + snapInt
            CALL openFile(30, SnapFile)
            CALL SnapShot(VEL = v_j, DENS_n = rhon, DENS_j = rhoj, DENSu_n = rho_un, &
            DENSu_j = rho_uj, VOID_n = agn, VOID_j = agj, intENERGY_n = uln, &
            intENERGY_j = ulj, QUAL_n = Xn, QUAL_j = Xj, TEM_n = Tn, TEM_j = Tj, &
            PRESS = Psys, VOID_DENS11 = al_rhol11, HEATq = Q0_f, &
            TIME = el_time, VALVE = tcvPos, DELTA_t = dt, FEEDFLOW = mfeed, &
            MY_COLUMN = 1, FILE_NUM = 30)
            CLOSE(UNIT = 30, STATUS = 'KEEP')
        END IF

        !Increasing time loop counter
        tloop = tloop + 1

        ! Calculating Xenon
        CALL XenonTrack2(N_Xe, N_I, avg_relPwr, dt, tloop, PoisonEquil)
        RelXe = N_Xe / nomXe

        !---> Getting Current Power, Pump, and Steam Demand info
        pPower = cPower
        CALL PumpProf(cPump,Pump,t_p,el_time)
        ! During a startup, reactor power is controlled directly
        IF(StartUpCond) THEN
            CALL PowerProf(Q0_f,q0hot,cPower,Power,t_q,el_time)
            Q0_Rx = Q0_f
        ELSE ! Otherwise, the power demand is set by the steam output demand
            CALL SteamProf(cSteamD, Steam, t_st, el_time)
            ! CALL Steam_Pwr(cSteamD, msteam, cPower, pPower, Q0_f, q0hot, dt, &
            !                SteamFollowG = stFollowG)

            ! Testing Reactivity Subroutine
            CALL React_Cont(cSteamD, msteam, react, react_CR, avg_relPwr, cPower, cC, avg_void, &
                            Q0_f, q0hot, el_time, dt, tloop, rho_ex, rho_fuel, rho_void, TFD, &
                            cont_meth, G_feed, RelXe, N_Xe, N_I, SteamFollowG = stFollowG)

            ! Calculating the Heat Going into the fluid
            Q0_Rx = HeatTransfer(Q0_Rx,Q0_f,dt)
        END IF

        !--> Updating the feed enthalpy
        hFD    = hliq(TFD,Psys(1,1)) ![Btu/lbm]
        !Heat Generation rate array
        DO n = 1,SIZE(qn)
            IF( (n .GE. 2) .AND. (n .LE. 7) ) THEN !Core
                qn(n,1) = q(zn(n - 1),Q0_Rx) * Ln(n)
            ELSE !Not Core
                qn(n,1) = 0. ![Btu/hr]
            END IF
        END DO

        !Updating Turbine Valve
        Ktcv(1) = Ktcv(2) !Old time coefficient is equal to the last solved coefficient
        IF(SysAtPwr) THEN ! System is considered at power
            Gain  = TCV_Gain_FP
            Gainf = TCV_Special_FP
        ELSE ! System is still starting up
            Gain  = TCV_Gain_SU
            Gainf = TCV_Special_SU
        END IF
        CALL TCVControl(Psys(1,1),Pprev,Ktcv(1),Ktcv(2),tcvPos,tcvVel,dt,Gain,Gainf)

        !Sets steam line cross-sectional area to zero when there is no void
        IF(tcvPos .LE. 0) THEN
            Aj(12)  = 0.   ![ft^2]
        ELSE
            Aj(12)  = Aj12 ![ft^2]
            tcvCount = tcvCount + 1
        END IF

        !---> Updating Feed
        !   G1: Level Gain
        !   G2: Feed/Steam Mismatch Gain
        !   G3: Total Gain
        IF(SysAtPwr) THEN
            G1  = FCV_Gain1_FP
            G2  = FCV_Gain2_FP
            G3  = FCV_Gain3_FP
        ELSE
            G1  = FCV_Gain1_SU
            G2  = FCV_Gain2_SU
            G3  = FCV_Gain3_SU
        END IF
        CALL FeedControl(agn(11,1),mfeed,msteam,mfeednew,G1,G2,G3)

        ! Defines the new feed mass flow rate
        mfeed = mfeednew ![lbm/hr]

        !---> Finding a system total loss coefficient
        Ksys = 0
        DO n = 1,SIZE(Kn)
            Ksys = Ksys + Kn(n) * Psi_loss(Xn(n,2),Psys(1,2))
        END DO

        testcounter = 0
        kloop       = 0
        loop        = 0
        conv        = .FALSE.
        DO WHILE( conv .EQV. .FALSE. ) !Start convergence (k) loop

            ! kloop keeps track of the number of k loops taken to converge
            loop  = loop + 1
            kloop = kloop + 1


            !Eliminates zero-velocity guess when the TCV first opens
            IF((tcvCount .EQ. 1) .OR. tcvPos .GT. 0) THEN
                IF(v_j(12,2) .EQ. 0) v_j(12,2) = v_j(11,2)*0.01
                IF(v_j(12,1) .EQ. 0) v_j(12,1) = v_j(11,1)*0.01
            END IF

            !===============================================================================
            !                             MATRIX VALUES
            !===============================================================================

            !---> Assign a, aprev, b and bprev matrix values
            a12  = (rhog(Psys(1,1)) * Aj(12) / Vn(11) ) * dt                 ![lbm*hr / ft^4]
            a11  = (rhol(uliq(Tj(11,1)),Psys(1,1)) * Aj(11) / Vn(11)) * dt   ![lbm*hr / ft^4]
            a9   = (rhoj(9,1) * Aj(9) / Vn(11) ) * dt                        ![lbm*hr / ft^4]

            b12  = ((rhog(Psys(1,1))*ug(Psys(1,1)) + Psys(1,1)*ff) * Aj(12) / Vn(11) ) * dt                 ![Btu*hr / ft^4]
            b11  = ((rhol(uliq(Tj(11,1)),Psys(1,1))*uliq(Tj(11,1)) + Psys(1,1)*ff) * Aj(11) / Vn(11) ) * dt ![Btu*hr / ft^4]
            b9   = ((rho_uj(9,1) + Psys(1,1)*ff) * Aj(9) / Vn(11) ) * dt                                    ![Btu*hr / ft^4]

            CALL Assign_a(a,aprev,rhoj,Aj,Vn,dt)             ![lbm*hr / ft^4]
            CALL Assign_b(b,bprev,rho_uj,Aj,Vn,dt,Psys(1,1)) ![Btu*hr / ft^4]

            !---> Assign alpha matrix values
            alpha = 0
            CALL Assign_alpha(alpha,rhon,Tn,Psys(1,2),a,aprev,b,bprev)
            !alpha units (single phase): [Btu*hr / lbm*ft]
            !alpha units (two-phase)   : [hr / ft]

            !---> Calculating Betas for the last column of the A matrix
            CALL Assign_beta(beta,rhon,Tn,Psys(1,2),agn)
            !beta units (single phase): [Btu/(lbm*psia)]]
            !beta units (two-phase)   : [1/psia]

            !---> Calculating Steam Dome Values (M,E,LM)
            !SE, SM values
            CALL Assign_SMSE(SM,SE,rhon,rho_un,qn,dt,Vn,agj,rhoj,Gj,Xj,Tj,Psys,Aj,mfeed,hFD)
            !SM units: [lbm/ft^3]
            !SE units: [Btu/ft^3]

            !M matrix
            CALL Assign_M(M,Tn,Psys,agn,a9,a11,a12,SM)

            !E matrix
            CALL Assign_E(E,Tn,Psys,agn,b9,b11,b12,SE)

            !LM matrix
            al9   = ((1 - agj(9,1))*rhol(uliq(Tj(9,1)),Psys(1,1)) * Aj(9) / Vn(11) ) * dt ![lbm*hr / ft^4]
            S_Ml11 = al_rhol11(1,1) + mfeed * dt / Vn(11) ![lbm / ft^3]
            CALL Assign_LM(LM,Tn,Tj,Psys,agn,agj,Aj,al9,a11,S_Ml11)

            !---> Assign psi matrix, beta(11) values
            CALL Assign_Psi(psi,M,E,LM)

            !---> Calculating Xi
            CALL Assign_Xi(Xi,M,E,LM,tcvpos,rhon,Tn,Psys,agn,SM,SE,uln)
            !Xi units (n=1..10) (single phase): [Btu / lbm]
            !Xi units (n=1..10) (two-phase)   : [dimensionless]


            !---> Calculating C12 and Xi_cond
            C12     = (Ktcv(2)+Ksys)*rhoj(12,1)*v_j(12,2) / gc  ![psf]

            Xi_cond = ((Psys(1,2)*144 + (Ktcv(2)+Ksys)*(rhoj(12,1)*v_j(12,2)**2) / (2*gc) - Pcond*144))  ![psf]

            !---> Calculating Gammas
            Do n = 1,SIZE(gam)
                gam(n) = (Ln(n))*rhoj(n,1)/dt + & ! Changed to new time
                DP_friction(v_j(n,2),rhon(n,:),rhoj(n,1),Xj(n,1),Tj(n,:), &
                Psys(1,2),De_n(n),Ln(n)) + &
                DP_loss(v_j(n,2),rhon(n,1),rhoj(n,1),Xj(n,1),Tj(n,1), &
                Psys(1,2),Kn(n))
            END DO

            SME = 0.

            !---> Calculating SME1, SME2
            CALL Assign_SME(SME,Gj,Xj,Xn,agj,Psys,rhoj,rhon,Ln,v_j,Tj,Tn,De_n,Kn,dt)

            !There seems to be a 'rounding' error that is causing problems
            IF(ABS(SME(1,1)) .LE. trunc) SME(1,1) = 0.
            IF(ABS(SME(2,2)) .LE. trunc) SME(2,2) = 0.
            ! IF(ABS(C12) .LE. trunc) C12 = 0.


            !===============================================================================
            !                            WRITING 'A' MATRIX
            !===============================================================================
            !-> Assign A matrix to all 0s
            Amat = 0.

            !-> Filling A matrix
            ! Diagonal
            i  = SIZE(Amat,1)
            DO o = 1,i !Row
                DO n = 1,i !Column
                    !Setting 'A' matrix equal to the alpha matrix
                    IF((o .LT. i-1) .AND. (n .LT. i-1)) THEN !Fill in Alpha matrix
                        Amat(o,n) = alpha(o,n)
                    ELSE IF((n .EQ. i) .AND. (o .LT. i-1)) THEN !Fill in Betas in the last column
                        Amat(o,n) = beta(o)
                    END IF

                    !Filling Row 12 for gamma dealing with pressure drops
                    IF((o .EQ. 12) .AND. (n .LE. 11)) THEN
                        Amat(o,n) = gam(n)
                    END IF
                END DO

            END DO

            !Irregular Values
            Amat(11,9)  = -psi(3)
            Amat(11,11) = psi(2)
            Amat(11,12) = psi(1)
            Amat(11,13) = beta(11)
            IF(tcvPos .LE. 0) THEN !Condition when the turbine control valve is closed
                Amat(13,12) = 1 !This sets vj12 = 0
            ELSE !This happens whent the turbine control valve is open
                Amat(13,12) = C12
                Amat(13,13) = -1
            END IF

            !Writing the 'A' matrix to a csv file (Output.csv)
            IF(WantDebug)THEN
                WRITE(11,'(A,F10.6, A,I2, A,F10.6, A,F14.8, A,F11.8)') "Elapsed Time = ", el_time,&
                "  ||  k loop = "   , loop , &
                "  ||  q0 = "       , Q0_Rx , &
                "  ||  Psys(1,2) = ", Psys(1,2) , &
                "  ||  dt = ", dt
                WRITE(11,*)
                i  = SIZE(Amat,1)
                DO o = 1,i !Row
                    DO n = 1,i !Column
                        !Writing matrix to csv
                        WRITE(11,'(F20.8,A)',ADVANCE='NO') Amat(o,n), ", "
                    END DO
                    WRITE(11,*)
                END DO
            END IF
            ! WRITE(*,'(A)') "A matrix has been written to 'Output.csv'"


            !===============================================================================
            !                            WRITING 'B' MATRIX
            !===============================================================================
            !Writing it in a matrix
            DO n = 1,SIZE(Bmat,1)
                IF(n .EQ. 12) THEN
                    Bmat(n,1) = (SME(1,1) + SME(2,2)) + cPump * (144*gc)
                ELSE IF(n .EQ. 13) THEN
                    IF(tcvPos .LE. 0) THEN
                        Bmat(n,1) = 0
                    ELSE
                        Bmat(n,1) = Xi_cond
                    END IF

                ELSE
                    Bmat(n,1) = Xi(n)
                END IF
            END DO

            !Writing the B matrix to a csv file (Output.csv)
            IF(WantDebug) THEN
                i  = SIZE(Bmat,1)
                WRITE(11,*)
                DO n = 1,i !Column
                    !Writing matrix to csv
                    WRITE(11,'(F20.8,A)',ADVANCE='NO') Bmat(n,1), ", "
                END DO
                WRITE(11,*)
            END IF

            !===============================================================================
            !                         SOLVE SYSTEM OF EQUATIONS
            !===============================================================================
            !Solve linear system of equations, B matrix overwritten with solution matrix
            CALL DGESV(13,1,Amat,13,IPIV,Bmat,13,stat)
            IF(stat .GT. 0) THEN
                WRITE(*,'(A,I2)') "Singular, no solution. Current iteration: ", loopde
                ! STOP
            ELSE IF(stat .LT. 0) THEN
                WRITE(*,'(A,I1,A)') "Illegal entry in argument ", stat, "."
                STOP
            ENDIF

            DO n = 1,SIZE(Bmat,1)
                IF(ABS(Bmat(n,1)) .LE. trunc) Bmat(n,1) = 0.
            END DO
            !===============================================================================
            !                            DETERMINING CONVERGENCE
            !===============================================================================
            IF(WantDebug)THEN
                WRITE(15,'(A,F10.7,A,I6,A,F6.3,A,F10.4,A,F15.12)') "Elapsed Time = ", el_time,&
                "  ||  k loop = "   , loop, &
                "  ||  q0     = "       , Q0_Rx , &
                "  ||  Psys   = ", Psys(1,2), &
                "  ||  dt     = ", dt
                CALL writeVel2(v_j,1,el_time,Psys(1,1),dt,15)
                CALL writeVel2(v_j,2,el_time,Psys(1,2),dt,15)
            END IF

            IF(tloop .EQ. 1) dPrev = Bmat(13,1)


            IF(SysAtPwr) THEN
                CALL TimeControl(dtmin,dtmax,dt,v_j,Bmat,Psys(1,2),Pprev,dPrev, &
                                 count,conv,test,kloop,dt_change,rhon, tStep_tol_FP)
            ELSE
                CALL TimeControl(dtmin,dtmax,dt,v_j,Bmat,Psys(1,2),Pprev,dPrev, &
                                 count,conv,test,kloop,dt_change,rhon, tStep_tol_SU)
            END IF


            !===============================================================================
            !                UPDATING VARIABLES FOR K vs K+1 CONVERGENCE
            !===============================================================================
            ! Writing to debugging file
            IF(WantDebug) CALL writeVel2(v_j,2,el_time,Psys(1,2),dt,15)


            !Writing Velocities to the Output.csv
            IF(WantDebug) THEN
                WRITE(11,*)
                DO n = 1,SIZE(Bmat,1) !Column
                    !Writing matrix to csv
                    WRITE(11,'(F20.8,A)',ADVANCE='NO') Bmat(n,1), ", "
                END DO
                WRITE(11,*)
                WRITE(11,*)
                WRITE(11,*)
            END IF

            !---> Calculating a dPsi
            !Core to Seperator
            DO n = 2,9
                CALL   drhodP(rhon(n,2),y,Tn(n,2),Psys(1,2),agn(n,2))
                CALL drhodpsi(rhon(n,2),w,Tn(n,2),Psys(1,2))
                dPsi_n(n) = ( -(dP/144)*y - a(n)*Bmat(n,1) &
                + aprev(n)*Bmat(n-1,1) &
                + SM(n) - rho(rhon(n,2), uln(n,2), Psys(1,2), agn(n,2))) * (1/w)

                !Alternate methods of calculation
                ! CALL   drho_udP(rhon(n,2),y,Tn(n,2),Psys(1,2),agn(n,2))
                ! CALL drho_udpsi(rhon(n,2),w,Tn(n,2),Psys(1,2))
                ! dPsi_n(n) = (-(dP/144)*y - b(n)*v_j(n,2) + bprev(n)*v_j(n-1,2) + SE(n) - (rho_un(n,2))) &
                !             *(1/w)
            END DO

            !Steam Dome
            ! IF(tcvPos .LE. 0)THEN
            !     CALL   drhodP(rhon(11,2),y,Tn(11,2),Psys(1,2),agn(11,2))
            !     CALL drhodpsi(rhon(11,2),w,Tn(11,2),Psys(1,2))
            !     dPsi_n(11) = (- (dP/144)*y - a(11)*v_j(11,2) &
            !                + aprev(11)*v_j(9,2) &
            !                + SM(11) - rho(rhon(11,2), uln(11,2), Psys(1,2), agn(11,2))) * (1/w)
            ! ELSE
            CALL   drhodP(rhon(11,2),y,Tn(11,2),Psys(1,2),agn(11,2))
            CALL drhodpsi(rhon(11,2),w,Tn(11,2),Psys(1,2))
            ! d_ul11     = (-(M(2) - LM(2))*dP - (M(3) - LM(3))*v_j(12,2) - (M(4) - LM(4))*v_j(11,2) &
            !             + (M(5) - LM(5))*v_j(9,2) + (M(6) - LM(6))) * (1/(M(1) - LM(1)))
            !
            ! !The change in VOID FRACTION
            ! dPsi_n(11) = -d_ul11*M(1) - dP*M(2) - M(3)*v_j(12,2) - M(4)*v_j(11,2) + M(5)*v_j(9,2) + M(6)

            !Alternate methods of calculation
            d_ul11     = (-(E(2) - LM(2))*dP - (E(3) - LM(3))*v_j(12,2) - (E(4) - LM(4))*v_j(11,2) &
            + (E(5) - LM(5))*v_j(9,2) + (E(6) - LM(6))) * (1/(E(1) - LM(1)))
            !Delta_alpha in steam dome
            dPsi_n(11) = -d_ul11*E(1) - dP*E(2) - E(3)*v_j(12,2) - E(4)*v_j(11,2) + E(5)*v_j(9,2) + E(6)

            ! dPsi_n(11) = -d_ul11*LM(1) - dP*LM(2) - LM(3)*v_j(12,2) - LM(4)*v_j(11,2) + LM(5)*v_j(9,2) + LM(6)
            ! END IF

            !Down Comer
            CALL   drhodP(rhon(10,2),y,Tn(10,2),Psys(1,2),agn(10,2))
            CALL drhodpsi(rhon(10,2),w,Tn(10,2),Psys(1,2))
            dPsi_n(10) = (- (dP/144)*y - a(10)*Bmat(10,1) &
            + aprev(10)*Bmat(11,1) &
            + SM(10) - rho(rhon(10,2), uln(10,2), Psys(1,2), agn(10,2))) * (1/w)

            !Lower Plenum
            CALL   drhodP(rhon(1,2),y,Tn(1,2),Psys(1,2),agn(1,2))
            CALL drhodpsi(rhon(1,2),w,Tn(1,2),Psys(1,2))
            dPsi_n(1) = (- (dP/144)*y - a(1)*Bmat(1,1) &
            + aprev(1)*Bmat(10,1) &
            + SM(1) - rho(rhon(1,2), uln(1,2), Psys(1,2), agn(1,2))) * (1/w)

            !Setting the k value of the junction velocity array to B matrix
            i = SIZE(Bmat,1)
            DO n = 1,i
                IF(n .LT. i) THEN
                    v_j(n,2) = Bmat(n,1)
                ELSE
                    dP = Bmat(n,1)
                END IF
            END DO

            CALL UpdateProps_k(dP,Psys,dPsi_n,d_ul11,al_rhol11,Tn,Tj,agn,agj,rhon,rhoj,uln,ulj, &
            rho_un,rho_uj,Xn,Xj,Gj,v_j,qn,Vn,dt,dt_change,tcvPos)


            ! Writing to Debugging matrix
            IF(WantDebug) THEN
                CALL writeVel2(test,1,el_time,Psys(1,2),dt,15)
                WRITE(15,*) "Convergence = ", conv
                WRITE(15,*)
            END IF


            IF( conv .EQV. .TRUE. ) THEN
                ! conv  = .FALSE.
                count = 0 !Reset counter for next test
                ! EXIT
            END IF
            testcounter = testcounter + 1

            n = 300
            IF(kloop .GT. n) THEN
                WRITE(*,*) "No Convergence After: t = ", tloop
                WRITE(*,*) "No Convergence After: k = ", kloop
                STOP
            END IF

            dPrev = dP !Previoud iterates dP
        END DO !End convergence (k) loop

        ! Convergence Logical for updating values
        conv = .FALSE.

        ! Writing Beta Matrix to a dummy matrix for debugging
        DO n = 1,SIZE(Xi,1)
            dum1(n,1) = Xi(n)
        END DO

        !===============================================================================
        !                          Safety Calculations
        !===============================================================================
        IF(SysAtPwr) THEN ! Setting to NAN during startup
            CALL CHF(G_CORE = Gj(1,2), PRESS = Psys(1,2), q0_HOT = q0hot, &
                        U0 = ulj(1,2), CPR = CPR, TOLERANCE = CPRtol)
        ELSE
            CPR = 0
            CPR = 0/CPR
        END IF



        !===============================================================================
        !                                  Output
        !===============================================================================
        ! ---> Writing values to csv files to be read into MATLAB
        ! Calculting phase velocities in the 9th & 8th node
        CALL phase_vel(w,w2,Gj(9,2),Xj(9,2),agj(9,2),Psys(1,2),9) ! w=vg, w2=vl
        CALL phase_vel(y,y2,Gj(8,2),Xj(8,2),agj(8,2),Psys(1,2),8)

        Linear_Heat_Rate = (cPower * gf)/(n_frods * Hfuel)
        ! Writing to files
        y2 = kloop
        CALL writeVel2(v_j,2,el_time,Psys(1,2),dt,16)         ! Velocity
        CALL writeVel2(Tn,2,agn(11,1),TFD,cPower,17)          ! Internal Energy
        CALL writeVel2(rhon,2,el_time,tcvPos,cSteamD,18)      ! Density
        CALL writeVel2(agn,2,msteam,mfeed,react_CR,19)        ! Void Fraction
        CALL writeVel2(v_j,2,N_I,CPR,N_Xe,21)                  ! Relative Velocity
        CALL writeSingle(23,Linear_Heat_Rate,cPower)

        !---> Updating properties for the t loop
        CALL UpdateProps_t(Psys,al_rhol11,Tn,Tj,agn,agj,rhon,rhoj,uln,ulj, &
                           rho_un,rho_uj,Xn,Xj,Gj,v_j,Pprev)

        el_time = el_time + dt

        ! WRITE(*,*) Gj(1,1), Psys(1,1), q0hot, ulj(1,1)
        WRITE(*,*) "<----------------------------------------------"
        WRITE(*,*) cPower, cSteamD, el_time, tloop, agn(11,1)
        ! WRITE(*,*)
        ! IF(tloop .GE. 3) STOP
        !Writing the velocities to Vels.csv
        ! CALL writeVel(Bmat,SIZE(Bmat,1),el_time,loop,q0,Psys(1,2),SIZE(Tn,1), 2, Tn)

        ! Writing the last poison value to a file
        IF(writePoisonFile .NE. "") THEN
            IF(tloop .EQ. 1) CALL openFile(65, writePoisonFile)
            WRITE(65,*) N_Xe, N_I
            REWIND(65)
        END IF
    END DO  !End time loop

END PROGRAM
