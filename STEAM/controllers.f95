! Feed Water and Turbine Control Valve (TCV) Controllers Module
MODULE controllers

    USE globals
    USE files
    IMPLICIT NONE

CONTAINS

    !Turbine Control Valve (TCV) controller subroutine
    SUBROUTINE TCVControl(P,Pprev,K,Knew,tcvPos,tcvVel,dt,G,Gf)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: P, Pprev, K, dt, Gf
        REAL(8),INTENT(OUT):: Knew
        REAL(8),INTENT(INOUT):: tcvPos, tcvVel, G
        REAL(8):: Error, K0, DK, tcvPos_dum, max_tcvVel, tcvVel_dum, &
                  max_tcvAcc, tcvAcc

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  P          == Instantaneous (current) pressure
        !                   Pprev      == Previoues iteration pressure
        !                   K          == Instantaneous (current) TCV loss
        !                                 coefficient
        !                   dt         == Change current change in time
        !                   G          == Gain on the turbine control valve
        !                   Gf         == Gain fraction when pressure is moving toward the reference
        !                   Knew       == New (t + dt) TCV loss coefficient
        !                   Pref       == Reference (system) pressure
        !                   Error      == TCV error signal
        !                   K0         == Full-open loss coefficient
        !                   Taunew     == New (t + dt) relative valve position [0,1]
        !                   G          == Error gain
        !                   DK         == Change in TCV loss coefficient
        !                   tcvPos     == Turbine control valve position
        !                   tcvVel     == Maximum rate of change for the turbine control valve
        !                   tcvPos_dum == intermediate turbine control valve
        !                                 position for determining the rate
        !                   tcvVel_fum == intermediate change in valve position
        !                                 to calc acceleration of the valve
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

!Second way
        !Full Open Loss Coefficient
        K0     = 350 !390.102 should give full opern to be at full power

        ! Limits on Velocity and Acceleration of the valve
        !   The values are currently large enough to not effect the TCV
        max_tcvVel = 1100!1e-4
        max_tcvAcc = 1e4!30

        ! Normalized Error Value
        Error = ((P - Pref)/Pref)
        G     = G * Error
        ! IF(Error .LT. 0) THEN ! Closing
        !     G  = 15 * Error ! 5
        ! ELSE IF(Error .GT. 0) THEN ! Opening
        !     G  = 5 * Error ! 5
        ! ELSE
        !     G  = 0.
        ! END IF

        !--> Checks for determining if the tcv needs to change
        !IF the valve is alread at 0. It cant go down
        IF((tcvPos .LE. 0) .AND. (G .LE. 0)) RETURN
        ! When the pressure is decreasing the valve shouldnt keep opening
        IF((P-Pprev .LE. 0) .AND. (G .GT. 0)) G = G * Gf !0.6
        ! WHen the pressure is increasing the valve shouldnt continue to close
        IF((P-Pprev .GE. 0) .AND. (G .LT. 0)) G = G * Gf !0.6

        !If the pressure is less than Pref, the valve will close.
        !If the pressure is greater than Pref, the valve will open
        IF((tcvPos .LE. 0) .AND. (G .GT. 0)) THEN
            tcvPos_dum = tcvPos + ABS(Error)*0.0001
        ELSE IF((tcvPos .LE. 0) .AND. (G .LE. 0)) THEN
            tcvPos_dum = 0
        ELSE
            tcvPos_dum = tcvPos * (1 + G)
        END IF

        !Making sure tcvPos isnt too small
        IF(tcvPos_dum .LE. 1e-11)THEN
            tcvPos = 1e-11
            RETURN
        END IF

        !Calculating the change in valve position over time
        tcvVel_dum = (tcvPos_dum - tcvPos)/dt

        !Calculating the acceleration of the valve
        tcvAcc = (tcvVel - tcvVel_dum)/dt

        !Checking the acceleration of the TCV
        !  If the acceration is too high, the dummy tcvPos is adjusted
        IF(( ABS(tcvAcc) .GT. max_tcvAcc ) .AND. (tcvPos .GT. 0 )) THEN
            tcvVel_dum = tcvVel + max_tcvAcc * dt * (tcvVel_dum/ABS(tcvVel_dum))
        END IF


        !Checking to see if the change in time is greater than the allowed max
        IF((ABS(tcvVel_dum) .GT. max_tcvVel) .AND. (tcvPos .GT. 0 )) THEN
            tcvPos = tcvPos + max_tcvVel*dt*(tcvVel_dum/ABS(tcvVel_dum))
        ELSE
            tcvPos = tcvPos + tcvVel_dum * dt
        END IF

        tcvVel = tcvVel_dum

        ! Checking to make sure the tcvPos is not less than 0 or greater than 1
        IF(tcvPos .LT. 0. )THEN !.OR. tcvPos .LT. 1e-20) THEN
            tcvPos = 0.
            ! WRITE(*,'(A)') "TCV is fully closed"
        ELSEIF(tcvPos .GE. 1) THEN
            ! tcvPos = 1.
            WRITE(*,'(A)') "TCV is fully open"
        END IF

        !Preventing Division by 0
        IF(tcvPos .GT. 0.) THEN
            Knew = K0/(tcvPos**2)
        ELSE
            Knew = 0
        END IF

        !Note: K(j) -> as the valve is closed. j = [0=closed,1=open]
    END SUBROUTINE


    !Feed flow controller subroutine
    SUBROUTINE FeedControl(level,mfeed,msteam,mfeednew,G1,G2,G3)
        IMPLICIT NONE

        REAL(8),INTENT(INOUT):: level
        REAL(8),INTENT(IN):: mfeed, msteam, G1, G2, G3
        REAL(8),INTENT(OUT):: mfeednew
        REAL(8):: levelref, Error1, Error2, Error, mfeednom, Dfeed

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  level    == Instantaneous (current) level in the
        !                               downcomer/steam dome
        !                   mfeed    == Instantaneous (current) feed flow rate
        !                   msteam   == Instantaneous (current) steam flow rate
        !                   mfeednew == Adjusted (t + dt) feed flow rate
        !                   levelref == Reference level in the downcomer/steam
        !                               dome
        !                   Error1   == Level error signal
        !                   Error2   == Feed error signal
        !                   Error    == Combined error signal
        !                   mfeednom == Nominal (full power) feed flow rate [lbm/hr]
        !                   G1,G2,G3 == Error gains
        !                   Dfeed    == Change in feed flow rate
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        ! Putting the level in terms of conventional level instead of void
        level    = 1 - level
        levelref = 0.5 !Void in the Steam Dome
        mfeednom = 4e6 !6.099291757605433e+06 ! This is very unknown

        Error1 = ( levelref - level )/levelref



        IF(msteam .GT. 0) THEN
            Error2 = ( msteam - mfeed )/msteam
        ELSEIF(mfeed .GT. 0) THEN
            Error2 = ( msteam - mfeed )/mfeed
        ELSE
            Error2 = 0
        END IF

        WRITE(*,*)
        WRITE(*,*) "Error 1 -->", Error1
        WRITE(*,*) "Error 2 -->", Error2

        Error  = G1*Error1 + G2*Error2

        Dfeed  = G3*Error
        mfeednew = mfeed + Dfeed !mfeed * (1 + Dfeed)
        ! IF((mfeed .LE. 0) .AND. (Error1 .LE. 0)) THEN
        !     mfeednew = msteam
        ! ELSE
        !     mfeednew = mfeed * (1 + Dfeed)
        ! END IF
        WRITE(*,*) "Dfeed     -->", Dfeed
        WRITE(*,*) "Old mfeed -->", mfeed
        WRITE(*,*) "New mfeed -->", mfeednew
        IF(mfeednew .LT. 0) mfeednew = 0
        WRITE(*,*)
    END SUBROUTINE


    SUBROUTINE TimeControl(dtmin,dtmax,dt,vtk,vk1,Press,Pprev,dPrev,count, &
                            converged,test,kloop,dt_change,rhon, tol )
        IMPLICIT NONE

        !I/O Variables
        INTEGER,INTENT(INOUT) :: count, kloop
        LOGICAL,INTENT(INOUT) :: converged, dt_change
        REAL(8),INTENT(IN)    :: dtmin, dtmax, Press, rhon(:,:), Pprev, dPrev, tol
        REAL(8),INTENT(INOUT) :: vtk(:,:), vk1(:,:), dt
        REAL(8),INTENT(OUT)   :: test(:,:)

        !Local Variables
        REAL(8) :: maxeps, tol2, tol3
        REAL(8),ALLOCATABLE:: eps(:,:)
        INTEGER :: i, j, conv_n, s, c2p
        LOGICAL :: zerocase, zerocase2, LowTol = .FALSE.

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  count       == Number of times converged
        !                   kloop       == Current convergence iteration
        !                   i,j         == Counter variables
        !                   conv_n      == Index convergence counter
        !                   s           == Size of velocity array (dim 1)
        !**REAL**      ==>  dtmin       == Minimum time step [hr]
        !                   dtmax       == Maximum time step [hr]
        !                   vtk()       == Velocity matrix (t,k)
        !                   vk1()       == Velocity matrix (k+1)
        !                   dP          == Current system pressure increment
        !                   test()      == Relative difference in iterate velocities
        !                   tol         == Convergence tolerance
        !                   eps()       == Relative difference in time velocities
        !                   maxeps      == Maximum relative difference in time velocities
        !                   tol2        == Relative difference in time velocities tolerance
        !                   tol3        ==
        !                   dPrev       == Previous iteration's change in pressure
        !**CHARACTER** ==>
        !**LOGICAL**   ==>  converged   == .TRUE. if converged enough times
        !                   zerocase    == Initial case where all velocities are zero
        !                   zerocase2   == Initial case where old time velocities (t) are zero
        !----------------------------------!

        s = SIZE(vtk,1)
        ALLOCATE(eps(s,1))

        c2p = 0
        DO i = 1,SIZE(rhon,1)
            IF(rhon(i,1) .LE. rhof(Press)) c2p = c2p + 1
        END DO
        IF(c2p .GT. 1) LowTol = .TRUE.


        tol2    = 1e-3
        tol3    = 1e-6
        conv_n  = 0
        eps     = 0
        test    = 0.
        maxeps  = 0
        j       = 0
        zerocase  = .FALSE.
        zerocase2 = .FALSE.
        dt_change = .FALSE.

        !Adjust for possible rounding errors
        DO i = 1,SIZE(vtk,1)
            IF( ABS(vtk(i,1)) .LT. trunc ) vtk(i,1) = 0
            IF( ABS(vtk(i,2)) .LT. trunc ) vtk(i,2) = 0
            IF( ABS(vk1(i,1)) .LT. trunc ) vk1(i,1) = 0
        END DO

        !Determine if all new time velocities are zero
        DO i = 1, s
            IF( vtk(i,2) .EQ. 0 ) j = j + 1
        END DO
        IF( j .EQ. s ) zerocase = .TRUE.

        !Determining if all old time velocities are zero
        j = 0
        DO i = 1, s
            IF( vtk(i,1) .EQ. 0 ) j = j + 1
        END DO
        IF( j .EQ. s ) zerocase2 = .TRUE.

        !Evaluate relative difference of the iterate velocities (k, k+1)
        DO i = 1,SIZE(vk1,1)
            IF(i .NE. 13) THEN
                IF( (zerocase .EQV. .TRUE.) ) THEN
                    test(i,1) = ABS( vk1(i,1) - vtk(i,2) )
                ELSE IF(ABS(vtk(i,2)) .EQ. 0.) THEN
                    test(i,1) = ABS( (vk1(i,1) - vtk(i,2)) )
                ELSE
                    test(i,1) = ( vk1(i,1) - vtk(i,2))  / vtk(i,2)
                END IF
            ELSE ! Pressure
                IF(kloop .GT. 100) THEN
                    test(i,1) =  (vk1(13,1)/144)*0.1!((dPrev - vk1(13,1))/dPrev)*0.00004!(vk1(13,1)/144)*15
                ELSEIF(kloop .GT. 20) THEN
                    test(i,1) =  (vk1(13,1)/144)*0.1
                ELSE
                    IF(ABS(Press - Pprev) .LT. 0.001) THEN
                        test(i,1) = (vk1(13,1)/144)*1!((dPrev - vk1(13,1))/dPrev)*0.00004(vk1(13,1)/144)*15
                    ELSE
                        test(i,1) = (vk1(13,1)/144)*1!((dPrev - vk1(13,1))/dPrev)*0.00004!ABS((vk1(13,1)/144 + Press - Press))/Press
                    END IF
                END IF
            END IF
            !IF(test(i) .LE. tol) conv_n = conv_n + 1
        END DO
        ! test(13,1) = 0

        !If all relative differences converge, increment convergence counter
        IF( Maximum(ABS(test)) .LE. tol ) THEN
            !count = count + 1
            !If all relative differences have converged 3 times, fully converged
            ! IF( count .EQ. 3 ) THEN
                converged = .TRUE.
            ! ELSE
                ! converged = .FALSE.
            ! END IF
        ! ELSE
            ! count = 0 !Reset counter if another pass does not converge
        END IF


        ! IF( (count .GT. 1) .AND. (count .LT. 3) ) RETURN

        !Case to exit the subroutine when k and k+1 have converged and the t velocities are 0
        IF((converged .EQV. .TRUE.) .AND. (zerocase2 .EQV. .TRUE.)) RETURN

        !Evaluate relative differences in time step velocities (t, t+dt)
        DO i = 1, s
            IF( (zerocase .EQV. .TRUE.) .OR. (zerocase2 .EQV. .TRUE.) ) THEN
                eps(i,1) = ABS( vtk(i,2) - vtk(i,1) )
            ELSE
                eps(i,1) = ABS(( vk1(i,1) - vtk(i,1) ) / vtk(i,1))
            END IF
        END DO

        maxeps = Maximum(eps)
        ! WRITE(*,*)
        ! WRITE(*,*) maxeps

        !If fully converged
        IF( converged .EQV. .TRUE. ) THEN
            !For the non-zero case
            IF( ( zerocase .EQV. .FALSE. ) ) THEN
                !Difference is greater than the allowed tolerance
                IF( maxeps .GT. tol2 ) THEN
                    IF(dt/2. .GT. dtmin) THEN !When the loop actually changes time step
                        WRITE(*,'(A)')"Converged went back down due to epsilon "
                        dt_change = .TRUE.
                        converged = .FALSE.
                        kloop     = 0

                    ELSE !When the loop time step is at it's minimum
                        dt_change = .FALSE.
                        converged = .TRUE.
                    END IF
                    dt = MAX(dtmin , dt/2.)
                    ! WRITE(*,'(A)')"Lower"
                    RETURN

                !Difference is acceptable and needs to converge
                ELSE IF( (maxeps .GT. tol3) .AND. (maxeps .LE. tol2) ) THEN
                    dt = dt
                    converged = .TRUE.
                    ! WRITE(*,'(A)')"Middle"
                    RETURN

                !Difference is more than it needs to be in order to converge
                ELSE !(maxeps .LE. tol3)
                    dt = MIN(dtmax, 1.2*dt)
                    dt_change = .FALSE.
                    WRITE(*,*) "Converged and Increased Time Step <<----"
                    ! WRITE(*,'(A)')"Over"
                    RETURN
                END IF

            !For the zero case
            ELSE
                IF( maxeps .EQ. 0. ) THEN
                    dt = MIN(dtmax, 1.2*dt)
                    RETURN
                ELSE IF( (maxeps .GT. tol3) .AND. (maxeps .LE. 0.) ) THEN
                    dt = dt
                    converged = .FALSE.
                    RETURN
                ELSE
                    IF(dt/2. .GT. dtmin) dt_change = .TRUE.
                    dt = MAX(dtmin , dt/2.)
                    converged = .FALSE. !Repeat time step
                    RETURN
                END IF

            END IF

        !If not fully converged
        ELSE IF( (converged .EQV. .FALSE.) .AND. (kloop .EQ. looplim)) THEN
            !Difference is acceptable and needs to converge
            IF( (maxeps .GT. tol3) .AND. (maxeps .LE. tol2) ) THEN
                dt = dt
                converged = .FALSE.
                RETURN
            ELSE
                IF(dt/2. .GT. dtmin) THEN !When the loop actually changes time step
                    dt_change = .TRUE.
                    converged = .FALSE.
                    kloop     = 0
                ELSE !When the loop time step is at it's minimum
                    dt_change = .FALSE.
                END IF
                dt = MAX(dtmin , dt/2.)
            END IF

        !If loop limit not reached
        ELSE IF( (converged .EQV. .FALSE.) .AND. (kloop .LT. looplim)) THEN
            RETURN
        END IF

    END SUBROUTINE


    SUBROUTINE PowerProf(cq0,cq0hot,cPower,Power,t_q,t,Err)
        IMPLICIT NONE

        !Inputs
        REAL(8), INTENT(IN)  :: Power(:), t_q(:), t
        REAL(8), INTENT(OUT) :: cq0, cq0hot, cPower
        INTEGER, OPTIONAL    :: Err

        !Calculations
        INTEGER :: s, n
        REAL(8) :: y, qflux

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  s           == number of manuevers
        !                   n           == Counting variable
        !                   Err         == Optional Error Handling Variable
        !                                  1 -> The times entered are not sequential
        !                                  2 -> The time enterend is not in range
        !                                  3 -> Both Errors 1 & 2
        !
        !**REAL**      ==>  Power(:)    == Power profile desired [W].
        !                                  The power given 'Power(n)' corresponds to the 't(n)'
        !                                  (ie, if t(2) = 1hr and Power(2) = 1000W, the reactor will be at
        !                                  1000W at an el_time of 1hr)
        !                                  | Q1 , Q2, Q3, Q4, Q5, ... |
        !                   t_q(:,:)    == Time [hr] corresponding to the power level desired from Power()
        !                                  The time given 't(n)' corresponds to the desired 'Power(n)'
        !                                  | t1 , t2 , t3 , t4 , t5 , ... |
        !                   t           == Current time (el_time)
        !                   y           == Dummy Variable
        !                   cPower      == Current Power Level
        !                   cq0         == Current q0 associated with the power level
        !                   cq0hot      == Current q0 (heat flux) for the hot channel
        !                   qflux       == Average heat flux
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        s = SIZE(Power,1)
        IF(PRESENT(Err)) Err = 0
        y = 0

        ! Checking to ensure each time was imported sequentially
        DO n = 1,s
            !If the next start or end time is less than the previous start or end times
            IF (n .LT. s) THEN
                IF ((t_q(n+1) .LT. t_q(n)) .AND. PRESENT(Err)) THEN
                    Err = 1
                END IF
            END IF
        END DO

        !Finding which time span is used
        DO n = 1,s
            IF( (t .GE. t_q(n)) .AND. (t .LE. t_q(n+1)) .AND. (n .NE. s)) THEN
                cPower = (Power(n+1) - Power(n)) / (t_q(n+1) - t_q(n))*(t - t_q(n)) + Power(n)
                cq0    = (cPower/(shape_int(Hfuel))) * (3.41214)
                qflux  = cPower*3.41214*gf/(n_frods*pi*Hfuel*D_o)
                cq0hot = Fq*qflux
                RETURN
            ELSE !The time entered is not in range (used if elapsed time goes past defined time range)
                cPower = Power(n)
                cq0    = (cPower/(shape_int(Hfuel))) * (3.41214)
                qflux  = cPower*3.41214*gf/(n_frods*pi*Hfuel*D_o)
                cq0hot = Fq*qflux

                IF(PRESENT(Err)) Err = Err + 2
            END IF
        END DO

    END SUBROUTINE

    SUBROUTINE PumpProf(cPump,Pump,t_p,t,Err)
        IMPLICIT NONE

        !Inputs
        REAL(8), INTENT(IN)  :: Pump(:), t_p(:), t
        REAL(8), INTENT(OUT) :: cPump
        INTEGER, OPTIONAL    :: Err

        !Calculations
        INTEGER :: s, n
        REAL(8) :: y

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  s           == number of manuevers
        !                   n           == Counting variable
        !                   Err         == Optional Error Handling Variable
        !                                  1 -> The times entered are not sequential
        !                                  2 -> The time enterend is not in range
        !                                  3 -> Both Errors 1 & 2
        !
        !**REAL**      ==>  Pump(:)     == Pump profile desired [psia].
        !                                  The power given 'Pump(n)' corresponds to the 't(n)'
        !                                  (ie, if t(2) = 1hr and Pump(2) = 10psia, the pump will be at
        !                                  10psia at an el_time of 1hr)
        !                                  | cPump1 , cPump2, cPump3, ... |
        !                   t_p(:)      == Time [hr] corresponding to the pump level desired from Pump()
        !                                  The time given 't(n)' corresponds to the desired 'Pump(n)'
        !                                  | t1 , t2 , t3 , ... |
        !                   t           == Current time (el_time)
        !                   y           == Dummy Variable
        !                   cPump       == Current pump Level
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        s = SIZE(Pump,1)
        IF(PRESENT(Err)) Err = 0
        y = 0

        ! Checking to ensure each time was imported sequentially
        DO n = 1,s
            !If the next start or end time is less than the previous start or end times
            IF (n .LT. s) THEN
                IF ((t_p(n+1) .LT. t_p(n)) .AND. PRESENT(Err)) THEN
                    Err = 1
                END IF
            END IF
        END DO

        !Finding which time span is used
        DO n = 1,s
            IF( (t .GE. t_p(n)) .AND. (t .LE. t_p(n+1)) ) THEN
                cPump = (Pump(n+1) - Pump(n)) / (t_p(n+1) - t_p(n))*(t - t_p(n)) + Pump(n)
                RETURN
            ELSE !The time entered is not in range (used if elapsed time goes past defined time range)
                cPump = Pump(n)
                IF(PRESENT(Err)) Err = Err + 2
            END IF
        END DO

    END SUBROUTINE

    SUBROUTINE SteamProf(cSteam,Steam,t_st,t,Err)
        IMPLICIT NONE

        !Inputs
        REAL(8), INTENT(IN)  :: Steam(:), t_st(:), t
        REAL(8), INTENT(OUT) :: cSteam
        INTEGER, OPTIONAL    :: Err

        !Calculations
        INTEGER :: s, n
        REAL(8) :: y

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  s           == number of manuevers
        !                   n           == Counting variable
        !                   Err         == Optional Error Handling Variable
        !                                  1 -> The times entered are not sequential
        !                                  2 -> The time enterend is not in range
        !                                  3 -> Both Errors 1 & 2
        !
        !**REAL**      ==>  Steam(:)     == Pump profile desired [psia].
        !                                  The steam demand given 'Steam(n)' corresponds to the 't(n)'
        !                                  (ie, if t(2) = 1hr and Steam(2) = 100lbm/hr, the
        !                                  steam demand will be at 100lbm/hr at an el_time of 1hr)
        !                                  | Steam1 , Steam2, Steam3, ... |
        !                   t_st(:)     == Time [hr] corresponding to the Steam Demand desired from Steam()
        !                                  The time given 't(n)' corresponds to the desired 'Pump(n)'
        !                                  | t1 , t2 , t3 , ... |
        !                   t           == Current time (el_time)
        !                   y           == Dummy Variable
        !                   cSteam      == Current Steam Demand
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        s = SIZE(Steam,1)
        IF(PRESENT(Err)) Err = 0
        y = 0

        ! Checking to ensure each time was imported sequentially
        DO n = 1,s
            !If the next start or end time is less than the previous start or end times
            IF (n .LT. s) THEN
                IF ((t_st(n+1) .LT. t_st(n)) .AND. PRESENT(Err)) THEN
                    Err = 1
                END IF
            END IF
        END DO

        !Finding which time span is used
        DO n = 1,s
            IF( (t .GE. t_st(n)) .AND. (t .LE. t_st(n+1)) .AND. (n .NE. s) ) THEN
                cSteam = (Steam(n+1) - Steam(n)) / (t_st(n+1) - t_st(n))*(t - t_st(n)) + Steam(n)
                RETURN
            ELSE !The time entered is not in range (used if elapsed time goes past defined time range)
                cSteam = Steam(n)

                IF(PRESENT(Err)) Err = Err + 2
            END IF
        END DO

    END SUBROUTINE

    SUBROUTINE React_Cont(cSteamD, msteam, react, react_CR, avg_relPwr, Pnew, Cnew, avg_void,  &
                          Q0, q0hot, el_time, dt, tloop, rho_ex, rho_fuel, rho_void, TFD, &
                          cont_meth, G_feed, RelXe, N_Xe, N_I, SteamFollowG, nomSteam )
        ! This controls reactor power based on reactivity
        IMPLICIT NONE

        !Inputs
        REAL(8), INTENT(IN):: msteam, cSteamD, dt, avg_relPwr, el_time
        INTEGER, INTENT(IN):: tloop
        INTEGER, INTENT(INOUT):: cont_meth
        REAL(8), OPTIONAL,INTENT(IN):: nomSteam, SteamFollowG, avg_void
        REAL(8), INTENT(INOUT):: Pnew, Cnew, react, react_CR, rho_ex, rho_fuel, &
                                rho_void, TFD, G_feed, RelXe, N_Xe, N_I

        !Outpus
        REAL(8), INTENT(OUT):: Q0, q0hot

        !Calculations
        ! INTEGER ::
        REAL(8) :: Error1, G1, Error2, G2, nSteam, qflux, cmSteamD, cRamp
        REAl(8):: l, betaR, Tou1, Tou2, bR, cR, lambdaR, alpha1, alpha2, C0, P0, &
                  reactD, RelPwr

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  cSteamD   == Current steam Demand as a fraction of nominal (input)
        !                   msteam    == Current Steam mass flow rate (input)   [lbm/hr]
        !                   RxPwr     == Current Reactor Power (input/output) [W]
        !                   Q0        == Heat Constant applied to the shaping function (output) [btu]
        !                   q0hot     == Hot channel heat flux constant [btu/ft^2]
        !                  [nomSteam] == Optional input for the full power steam flow [lbm/hr]
        !
        !                   Error1    == Error on the steam demand versus steam output
        !                   G1        == Gain on the Error1
        !                   nSteam    == Local nominal steam flow rate [lbm/hr]
        !                   qflux     == max average heat flux [btu/ft^2*hr]
        !                   cmSteam   == Current steam mass flow rate [lbm/hr]
        !
        !                   reactD    == The Current Reactivity Demand
        !                   react_CR  == Control Rod Reactivity
        !                   TFD       == Feedwater temperature [°F]
        !                   cont_meth == the method desired to follow load
        !
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        !Checking if there is an input for the nominal full power steam flow rate
        IF(PRESENT(nomSteam))THEN
            nSteam = nomSteam
        ELSE
            nSteam = NomSteamM
        END IF
        IF(PRESENT(SteamFollowG))THEN
            G1 = SteamFollowG
        ELSE
            G1 = 0.008
        END IF

        G2 = 0. ! Error on reactor power

        RelPwr   = (msteam / nSteam  + Pnew / Nom_Power)/2 !cSteamD !Pnew / Nom_Power

        ! Calculating a steam flow rate demand based on the relative steam flow desired
        cmSteamD = nSteam*cSteamD

        ! Error on the Desired steam flow rate versus the actual steam flow rate
        Error1 = (cmSteamD - msteam) / nSteam
        ! If using feed temp to control Rx Pwr, the error is on actual power
        IF(cont_meth .EQ. 2) THEN
            Error1 = (cSteamD - RelPwr)
        END IF

        ! Error on Reactor Power to the Steam Demand
        Error2 = (cSteamD - Pnew/Nom_Power)

        ! WRITE(*,*)
        ! WRITE(*,*) "Error1 --> ", Error1
        ! WRITE(*,*) "Error2 --> ", Error2
        ! WRITE(*,*)

        ! Calculating a new Reactivity Demand
        reactD   = (Error1 * G1 + Error2 * G2)

        CALL React_Dem(el_time, tloop, dt, avg_void, RelPwr, RelXe, react_CR, reactD, react, &
                       rho_ex, rho_fuel, rho_void, TFD, cont_meth, G_feed)

        IF(ABS(react) .LT. 1e-12) react = 0

        !---> Calculating a New Reactor Power
        ! Reactivity terms
        betaR   = 0.0065
        lambdaR = 0.08  ! [1/s]
        l       = 1e-4    ! [s] ?
        Tou1    = ((1 - react) * l)/(react - betaR)
        Tou2    = ((1 - react) * l)/(betaR)
        bR      =  lambdaR - 1/Tou1
        cR      = -lambdaR*(1/Tou1 + 1/Tou2)
        alpha1  = (1/2.0)*(-bR + SQRT(bR**2.0 - 4.0*cR))
        alpha2  = (1/2.0)*(-bR - SQRT(bR**2.0 - 4.0*cR))

        ! New Power and Precursor Concentration
        IF(tloop .EQ. 1) THEN
            Cnew = (betaR*Pnew)/(l*lambdaR)
        END IF
        P0    = Pnew
        C0    = Cnew ! (betaR*P0)/(l*lambdaR)
        Pnew  = ((P0)/(alpha2-alpha1)) * &
                (alpha1*EXP(alpha1*dt*3600) - alpha2*EXP(alpha2*dt*3600)) + &
                (lambdaR*(P0+C0)/(alpha2-alpha1)) * &
                (EXP(alpha1*dt*3600) - EXP(alpha2*dt*3600))
        Cnew  = ((C0)/(alpha2-alpha1)) * &
                (alpha1*EXP(alpha1*dt*3600) - alpha2*EXP(alpha2*dt*3600)) - &
                ((C0/Tou1-P0/Tou2)/(alpha2-alpha1)) * &
                (EXP(alpha1*dt*3600) - EXP(alpha2*dt*3600))

        Pnew = ABS(Pnew)
        Cnew = ABS(Cnew)
        ! WRITE(*,*) "Reactivity ->", react*1e5, "pcm"
        ! WRITE(*,*) "  Tou1 & 2 ->", Tou1, Tou2
        ! WRITE(*,*) "alpha1 & 2 ->", alpha1, alpha2
        ! WRITE(*,*) "   b and c ->", bR, cR
        ! WRITE(*,*) "     Power ->", P0, Pnew, (Pnew - P0)/Nom_Power
        ! WRITE(*,*) "     Conc  ->", C0, Cnew
        ! WRITE(*,*) "        CR ->", react_CR * 1e5
        ! STOP

        ! Ensuring reactor power can't surpass 900MW
        ! IF(Pnew .GT. Nom_Power + 1e3) THEN
        !     Pnew = Nom_Power + 1e3
        !     ! WRITE(*,*) "Steam Demand is too high "
        ! END IF

        ! Making sure the reactor can't change power level faster than the ramping limit
        cRamp     = ((Pnew - P0)/(dt * 60 * Nom_Power)) ! [%/min]
        IF(ABS(cRamp) .GT. RampLimit) THEN
            Pnew = P0 + RampLimit * dt * Nom_Power * 60 * (cRamp/ABS(cRamp))
            WRITE(*,*)
            WRITE(*,*) "Ramp Limit Reached"
            WRITE(*,*)
        END IF

        ! Calculating the Values needed for heat input in each node
        Q0    = (Pnew/(shape_int(Hfuel))) * (3.41214)
        qflux = Pnew*3.41214*gf/(n_frods*pi*Hfuel*D_o)
        q0hot = Fq*qflux

    END SUBROUTINE

    SUBROUTINE React_Dem(el_time, tloop, dt, avg_void, RelPwr, RelXe, rho_CR, rho_demand, &
                         rho_net, rho_ex, rho_fuel, rho_void, TFD, cont_meth, G_feed)
    ! This subroutine makes up for a reactivity demand using a method of reactivity control
    IMPLICIT NONE

    ! Inputs
    REAL(8),INTENT(IN):: avg_void, RelPwr, dt, RelXe, el_time
    INTEGER,INTENT(IN):: tloop
    INTEGER,INTENT(INOUT):: cont_meth
    REAL(8),INTENT(INOUT):: rho_CR, rho_ex, rho_fuel, rho_void, TFD, G_feed, rho_demand
    REAL(8),INTENT(OUT):: rho_net

    ! Local
    REAL(8):: alpha_fuel, alpha_void, alpha_Xe, init_CR, drhodt_lim, drhodt, &
              rho_CR0, rho_Xe
    LOGICAL:: cb = .FALSE.

    !-------Variable Definitions-------!
    !**INTEGER**   ==>
    !**REAL**      ==>  avg_void     == Average void in the core
    !                   RelPwr       == Relative Reactor Power
    !                   rho_ex       == Excess Reactivity [dK/K]
    !                   rho_fuel     == Reactivity from the fuel temperature [dK/K]
    !                   alpha_fuel   == Reactivity Coefficient from fuel temp increase [rho/relPwr]
    !                   rho_void     == Reactivity from the avg core void [dK/K]
    !                   alpha_void   == Reactivity coefficient from the average core void [rho/voidFract]
    !                   rho_net      == final reactivity
    !                   init_CR      == Initial Control Rod Reactivity [dK/K]
    !                   drhodt_limit == Limit on the speed of the reactivity insertion [(dK/K)/hour]
    !                   TFD          == Feed Temperature
    !                   cont_meth    == the method desired to be used to follow the demand
    !                                   1 = Control Rods
    !                                   2 = Feed Temperature
    !                                   Other = Control Rods
    !**CHARACTER** ==>
    !**LOGICAL**   ==>
    !----------------------------------!

    ! Assigning Values to the Coefficients to perform a reactivity balance
    alpha_void = -9000 * 1e-5 ! -3000
    alpha_fuel = -1200 * 1e-5
    alpha_Xe   = -2650 * 1e-5
    init_CR    = -4500  * 1e-5
    drhodt_lim =  CR_INS_LIMIT  * 1e-5 ! 4000 works 2800

    ! Limit on how small of a change can be made
    IF(ABS(rho_demand) * 1e5 .LT. 0.01) THEN
        rho_demand = rho_demand * 0
    END IF

    ! Limits on the feed temperature band
    IF((TFD .GT. 500) .OR. (TFD .LT. 350) ) THEN
        cont_meth = 1
    ! ELSEIF(ABS(rho_demand) .LT. 1e-8)THEN !<--- Testing switch back to other method
    !     cont_meth = 1
    !     cb        = .TRUE.
    END IF

    IF(tloop .EQ. 1) THEN ! If this is the first time the subroutine is entered
        rho_CR   = init_CR
        rho_fuel = RelPwr * alpha_fuel
        rho_void = avg_void * alpha_void
        rho_Xe   = RelXe * alpha_Xe
        rho_ex   = -(rho_void + rho_fuel + rho_CR + rho_Xe)

    ELSE
        !--> Feedback from void and fuel power
        rho_void = avg_void * alpha_void
        rho_fuel = alpha_fuel * RelPwr
        rho_Xe   = alpha_Xe   * RelXe

        !--> Changing an input to meet the desired reactivity demand
        drhodt = ABS(rho_demand)/dt
        IF(cont_meth .NE. 2) THEN ! Using Control rods to follow steam demand
            ! Emplimenting a limit on how fast the control rods can insert reactivity
            IF(drhodt .GT. drhodt_lim ) THEN
                rho_CR = rho_CR + drhodt_lim * dt * (ABS(rho_demand)/rho_demand)
            ELSE
                rho_CR = rho_CR + rho_demand
            END IF
        ELSE ! Using feed temperature
            TFD = TFD + (-1.0)*rho_demand*G_feed
        END IF
    END IF

    !--> Final Net Reactivity
    rho_net = rho_void + rho_fuel + rho_ex + rho_CR + rho_Xe
    ! rho_net = rho_demand

    ! Writing the Reactivity
    CALL writeSingle(22,rho_net,rho_ex,rho_void,rho_fuel,rho_CR,rho_Xe,el_time)

    IF(cb) cont_meth = 2

    END SUBROUTINE

    SUBROUTINE Steam_Pwr(cSteamD, msteam, RxPwr, pPower, Q0, q0hot, dt,&
                        SteamFollowG, nomSteam )
        IMPLICIT NONE

        !Inputs
        REAL(8), INTENT(IN):: msteam, cSteamD, pPower, dt
        REAL(8), OPTIONAL,INTENT(IN):: nomSteam, SteamFollowG
        REAL(8), INTENT(INOUT):: RxPwr

        !Outpus
        REAL(8), INTENT(OUT):: Q0, q0hot

        !Calculations
        ! INTEGER ::
        REAL(8) :: Error1, G1, nSteam, qflux, cmSteamD, cRamp

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  cSteamD   == Current steam Demand as a fraction of nominal (input)
        !                   msteam    == Current Steam mass flow rate (input)   [lbm/hr]
        !                   RxPwr     == Current Reactor Power (input/output) [W]
        !                   Q0        == Heat Constant applied to the shaping function (output) [btu]
        !                   q0hot     == Hot channel heat flux constant [btu/ft^2]
        !                  [nomSteam] == Optional input for the full power steam flow [lbm/hr]
        !
        !                   Error1    == Error on the steam demand versus steam output
        !                   G1        == Gain on the Error1
        !                   nSteam    == Local nominal steam flow rate [lbm/hr]
        !                   qflux     == max average heat flux [btu/ft^2*hr]
        !                   cmSteam   == Current steam mass flow rate [lbm/hr]
        !
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        !Checking if there is an input for the nominal full power steam flow rate
        IF(PRESENT(nomSteam))THEN
            nSteam = nomSteam
        ELSE
            nSteam = NomSteamM
        END IF
        IF(PRESENT(SteamFollowG))THEN
            G1 = SteamFollowG
        ELSE
            G1 = 0.008
        END IF

        ! Calculating a steam flow rate demand based on the relative steam flow desired
        cmSteamD = nSteam*cSteamD

        ! Error on the Desired steam flow rate versus the actual steam flow rate
        Error1 = (cmSteamD - msteam) / nSteam

        ! Calculating a new Reactor Power Level
        RxPwr  = RxPwr + RxPwr * Error1 * G1


        ! Ensuring reactor power can't surpass 900MW
        IF(RxPwr .GT. Nom_Power + 1e3) THEN
            RxPwr = Nom_Power + 1e3
            ! WRITE(*,*) "Steam Demand is too high "
        END IF

        ! Making sure the reactor can't change power level faster than the ramping limit
        cRamp     = ((RxPwr - pPower)/(dt * 60 * Nom_Power)) ! [%/min]
        IF(ABS(cRamp) .GT. RampLimit) THEN
            RxPwr = pPower + RampLimit * dt * Nom_Power * 60 * (cRamp/ABS(cRamp))
            WRITE(*,*)
            WRITE(*,*) "Ramp Limit Reached"
            WRITE(*,*)
        END IF

        ! Calculating the Values needed for heat input in each node
        Q0    = (RxPwr/(shape_int(Hfuel))) * (3.41214)
        qflux = RxPwr*3.41214*gf/(n_frods*pi*Hfuel*D_o)
        q0hot = Fq*qflux

    END SUBROUTINE

    SUBROUTINE XenonTrack(N_Xe, N_I, POWER, dt )
    ! This tracks the xenon concentration
        IMPLICIT NONE

        ! Input Variables
        REAL(8),INTENT(INOUT):: N_Xe, N_I
        REAL(8),INTENT(IN):: POWER, dt

        ! Local Variables
        REAL(8):: lambda_I, lambda_Xe, gamma_I, gamma_Xe, nFlux, rho_f, sigma_f, &
                  M_235, NA, Xe0, I0, sigma_Xe, rho_235, Eta, Vol_f, N_235


        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        ! Initial Values for the Concentrations
        Xe0 = N_Xe
        I0  = N_I

        ! Constants
        NA        = 6.0221409e23
        gamma_I   = 0.064
        gamma_Xe  = 0.0661 !0.0025
        lambda_I  = LOG(2.0)/6.57 ! [hr^-1]
        lambda_Xe = LOG(2.0)/9.10 ! [hr^-1]
        sigma_Xe  = 2.6e6  * 1e-24 * 0.00107639 ! barn --> ft^2
        sigma_f   = 585.   * 1e-24 * 0.00107639 ! barn --> ft^2
        rho_f     = 685.  ![lbm/ft^3] (FROM OAK RIDGE REPORT ORNL/TM-2000/351 TABLE 3.1)
        Eta       = 200. * 1.51857e-16 ! MeV --> btu
        M_235     = 235. * 0.00220462 ! grams/mol --> lbm/mol
        Vol_f     = pi*(D_o**2.0)/4.0 * Hfuel * n_frods
        N_235     = rho_235 * NA/M_235

        ! Calculating the Neutron Flux
        nFlux = ((POWER*W2Btu)/(Eta * N_235 * sigma_f * Vol_f))
        WRITE(*,*) "Flux --> ", nFlux

        ! Calculating Concentrations
        ! N_I = (exp(-lambda_I*dt)*((I0*M_235*lambda_I - NA*gamma_I*phi*rho_f*sigma_f)/ &
        !       (M_235*(lambda_Xe - lambda_I + phi*sigma_Xe)) + (NA*gamma_I*phi*rho_f* &
        !       sigma_f*exp(lambda_I*dt))/(M_235*(lambda_Xe - lambda_I + phi*sigma_Xe)))* &
        !       (lambda_Xe - lambda_I + phi*sigma_Xe))/lambda_I
        !
        ! N_Xe = exp(-lambda_I*dt)*((I0*M_235*lambda_I - NA*gamma_I*phi*rho_f*sigma_f)/ &
        !        (M_235*(lambda_Xe - lambda_I + phi*sigma_Xe)) + (NA*gamma_I*phi*rho_f* &
        !        sigma_f*exp(lambda_I*dt))/(M_235*(lambda_Xe - lambda_I + phi*sigma_Xe))) &
        !        - exp(-dt*(lambda_Xe + phi*sigma_Xe))*((I0*M_235*lambda_I*lambda_Xe - &
        !        M_235*Xe0*phi**2*sigma_Xe**2 - M_235*Xe0*lambda_Xe**2 + M_235*Xe0*lambda_I*&
        !        lambda_Xe + I0*M_235*lambda_I*phi*sigma_Xe + M_235*Xe0*lambda_I*phi*sigma_Xe &
        !        - 2*M_235*Xe0*lambda_Xe*phi*sigma_Xe + NA*gamma_Xe*phi**2*rho_f*sigma_f*sigma_Xe &
        !        - NA*gamma_I*lambda_I*phi*rho_f*sigma_f - NA*gamma_Xe*lambda_I*phi*rho_f*sigma_f &
        !        + NA*gamma_Xe*lambda_Xe*phi*rho_f*sigma_f)/(M_235*(lambda_Xe + phi*sigma_Xe)* &
        !        (lambda_Xe - lambda_I + phi*sigma_Xe)) + (NA*phi*rho_f*sigma_f* &
        !        exp(lambda_Xe*dt + phi*sigma_Xe*dt)*(gamma_I*lambda_I + gamma_Xe*lambda_I &
        !        - gamma_Xe*lambda_Xe - gamma_Xe*phi*sigma_Xe))/(M_235*(lambda_Xe + &
        !        phi*sigma_Xe)*(lambda_Xe - lambda_I + phi*sigma_Xe)))

        WRITE(*,*) "Before Calc N_Xe --> ", Xe0

        N_I = (exp(-lambda_I*dt)*((I0*M_235*lambda_I - NA*gamma_I*nFlux*rho_f*sigma_f)/ &
              (M_235*(lambda_Xe - lambda_I + nFlux*sigma_Xe)) + (NA*gamma_I*nFlux*rho_f* &
              sigma_f*exp(lambda_I*dt))/(M_235*(lambda_Xe - lambda_I + nFlux*sigma_Xe)))* &
              (lambda_Xe - lambda_I + nFlux*sigma_Xe))/lambda_I
        ! WRITE(*,*) "After Calc N_I --> ", N_I, " Delta of ", N_I - I0

        N_Xe = exp(-lambda_I*dt)*((I0*M_235*lambda_I - NA*gamma_I*nFlux*rho_f*sigma_f)/ &
               (M_235*(lambda_Xe - lambda_I + nFlux*sigma_Xe)) + (NA*gamma_I*nFlux*rho_f* &
               sigma_f*exp(lambda_I*dt))/(M_235*(lambda_Xe - lambda_I + nFlux*sigma_Xe))) &
               - exp(-dt*(lambda_Xe + nFlux*sigma_Xe))*((I0*M_235*lambda_I* &
               lambda_Xe - M_235*Xe0*(nFlux**2)*sigma_Xe**2 - M_235*Xe0*lambda_Xe**2 + M_235* &
               Xe0*lambda_I*lambda_Xe + I0*M_235*lambda_I*nFlux*sigma_Xe + M_235*Xe0* &
               lambda_I*nFlux*sigma_Xe - 2*M_235*Xe0*lambda_Xe*nFlux*sigma_Xe +  &
               NA*gamma_Xe*nFlux**2*rho_f*sigma_f*sigma_Xe - NA*gamma_I*lambda_I*nFlux* &
               rho_f*sigma_f - NA*gamma_Xe*lambda_I*nFlux*rho_f*sigma_f + NA*gamma_Xe* &
               lambda_Xe*nFlux*rho_f*sigma_f)/(M_235*(lambda_Xe + nFlux*sigma_Xe)* &
               (lambda_Xe - lambda_I + nFlux*sigma_Xe)) + (NA*nFlux*rho_f*sigma_f* &
               exp(lambda_Xe*dt + nFlux*sigma_Xe*dt)*(gamma_I*lambda_I + gamma_Xe*lambda_I - &
               gamma_Xe*lambda_Xe - gamma_Xe*nFlux*sigma_Xe))/(M_235*(lambda_Xe + nFlux*sigma_Xe)*&
               (lambda_Xe - lambda_I + nFlux*sigma_Xe)))
    END SUBROUTINE

    SUBROUTINE XenonTrack2(N_Xe2, N_I2, POW, dt, tloop, PoisonEquil )
    ! This tracks the xenon concentration
        IMPLICIT NONE

        ! Input Variables
        REAL(8),INTENT(INOUT):: N_Xe2, N_I2
        REAL(8),INTENT(IN):: POW, dt
        INTEGER,INTENT(IN):: tloop
        LOGICAL,INTENT(IN):: PoisonEquil

        ! Local Variables
        REAL(8):: Avag, gamma_I2, gamma_Xe2, lambda_I2, lambda_Xe2, &
                  sigma_Xe2, sigma_f2, rho_U235, Eta2, M_U235, Vol_f2, N_U235, &
                  nflux2, Xe02, I02, y
        ! lambda_I, lambda_Xe, gamma_I, gamma_Xe, nFlux, rho_f, sigma_f, &
        !           M_235, NA, Xe0, I0, sigma_Xe, rho_235, Eta, Vol_f, N_235


        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        ! Initial Values for the Concentrations
        Xe02 = N_Xe2
        I02  = N_I2

        ! Constants
        Avag       = 6.0221409e23
        gamma_I2   = 0.0639
        gamma_Xe2  = 0.0025
        lambda_I2  = LOG(2.0)/6.58
        lambda_Xe2 = LOG(2.0)/9.14
        sigma_Xe2  = 2.65e6 * 1e-24 * 0.00107639
        sigma_f2   = 585.1  * 1e-24 * 0.00107639
        rho_U235   = 684.83473 ![lbm/ft^3]10.97 g/cm³
        Eta2       = 200. * 1.51857e-16 ! MeV to btu
        M_U235     = 235 * 0.00220462 !gram/mol to lbm/mol
        Vol_f2     = (pi/4.0)*(D_o**2) * Hfuel * n_frods
        N_U235     = (rho_U235 * Avag)/M_U235
        ! ! NA        = 6.0221409e23
        ! ! gamma_I   = 0.064
        ! ! gamma_Xe  = 0.0025
        ! ! lambda_I  = LOG(2.0)/6.57 ! [hr^-1]
        ! ! lambda_Xe = LOG(2.0)/9.10 ! [hr^-1]
        ! ! sigma_Xe  = 2.6e6  * 1e-24 * 0.00107639 ! barn --> ft^2
        ! ! sigma_f   = 585.   * 1e-24 * 0.00107639 ! barn --> ft^2
        ! ! rho_235   = 685.  ![lbm/ft^3] (FROM OAK RIDGE REPORT ORNL/TM-2000/351 TABLE 3.1)
        ! ! Eta       = 200. * 1.51857e-16 ! MeV --> btu
        ! ! M_235     = 235. * 0.00220462 ! grams/mol --> lbm/mol
        ! ! Vol_f     = pi*(D_o**2.0)/4.0 * Hfuel * n_frods
        ! ! N_235     = rho_235 * NA/M_235
        !
        ! y = 1e1

        ! Calculating the Neutron Flux
        nflux2 = ((POW)/(Eta2 * N_U235 * sigma_f2 * Vol_f2))

        IF((PoisonEquil) .AND. (tloop .EQ. 1)) THEN
            I02  = (Avag*gamma_I2*nflux2*rho_U235*sigma_f2)/(M_U235*lambda_I2)
            Xe02 = (Avag*gamma_I2*nflux2*rho_U235*sigma_f2 + Avag*gamma_Xe2*nflux2*&
                    rho_U235*sigma_f2)/(M_U235*(lambda_Xe2 + nflux2*sigma_Xe2))
        END IF
        ! WRITE(*,*) "Iodine ---> ", I02
        ! WRITE(*,*) "Xenon  ---> ", Xe02
        ! STOP
        !
        ! ! Calculating Concentrations
        N_I2  = (exp(-dt*lambda_I2)*((I02*M_U235*lambda_I2 - Avag*gamma_I2*nflux2* &
                rho_U235*sigma_f2)/(M_U235*(lambda_Xe2 - lambda_I2 + nflux2*sigma_Xe2)) &
                + (Avag*gamma_I2*nflux2*rho_U235*sigma_f2*exp(dt*lambda_I2))/(M_U235*(lambda_Xe2 &
                - lambda_I2 + nflux2*sigma_Xe2)))*(lambda_Xe2 - lambda_I2 + nflux2*sigma_Xe2))/lambda_I2

        N_Xe2 = exp(-dt*lambda_I2)*((I02*M_U235*lambda_I2 - Avag*gamma_I2*nflux2* &
                rho_U235*sigma_f2)/(M_U235*(lambda_Xe2 - lambda_I2 + nflux2*sigma_Xe2)) &
                + (Avag*gamma_I2*nflux2*rho_U235*sigma_f2*exp(dt*lambda_I2))/(M_U235*(lambda_Xe2 &
                - lambda_I2 + nflux2*sigma_Xe2))) - exp(-dt*(lambda_Xe2 + nflux2*sigma_Xe2))* &
                ((I02*M_U235*lambda_I2*lambda_Xe2 - M_U235*Xe02*nflux2**2*sigma_Xe2**2 - &
                M_U235*Xe02*lambda_Xe2**2 + M_U235*Xe02*lambda_I2*lambda_Xe2 + I02*M_U235*lambda_I2*&
                nflux2*sigma_Xe2 + M_U235*Xe02*lambda_I2*nflux2*sigma_Xe2 - 2*M_U235*Xe02*lambda_Xe2*nflux2*&
                sigma_Xe2 + Avag*gamma_Xe2*nflux2**2*rho_U235*sigma_f2*sigma_Xe2 - Avag*gamma_I2*lambda_I2*nflux2*&
                rho_U235*sigma_f2 - Avag*gamma_Xe2*lambda_I2*nflux2*rho_U235*sigma_f2 + Avag*gamma_Xe2*lambda_Xe2*&
                nflux2*rho_U235*sigma_f2)/(M_U235*(lambda_Xe2 + nflux2*sigma_Xe2)*(lambda_Xe2 - &
                lambda_I2 + nflux2*sigma_Xe2)) + (Avag*nflux2*rho_U235*sigma_f2*exp(dt*lambda_Xe2 &
                + dt*nflux2*sigma_Xe2)*(gamma_I2*lambda_I2 + gamma_Xe2*lambda_I2 - &
                gamma_Xe2*lambda_Xe2 - gamma_Xe2*nflux2*sigma_Xe2))/(M_U235*(lambda_Xe2 + &
                nflux2*sigma_Xe2)*(lambda_Xe2 - lambda_I2 + nflux2*sigma_Xe2)))
    END SUBROUTINE

END MODULE
