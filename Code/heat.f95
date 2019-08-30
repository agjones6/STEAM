!Heat Transfer Module
MODULE heat

    USE globals
    USE state
    IMPLICIT NONE

CONTAINS

    REAL(8) FUNCTION HeatTransfer(Q_Rx,P_fuel,dt)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: Q_Rx, P_fuel, dt
        REAL(8)           :: T_fuel, UA, fuel_dens, fuel_Cp, V_fuel, M_fuel, Tau!qflux_avg, q0_a, q0_h, q_vol

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  Q_Rx      == Current reactor power [Btu/hr]
        !                   P_fuel    == Next time-step fuel power (from neutronics model) [Btu/hr]
        !                   dt        == Current time step [hr]
        !                   T_fuel    == Average fuel tempearture [°F]
        !                   UA        == Reactor heat transfer constant [Btu/(hr*°F)]
        !                   fuel_dens == Density of the fuel [lbm/ft^3]
        !                   fuel_Cp   == Cp of the fuel [Btu/(lbm*°F)]
        !                   V_fuel    == Volume of the fuel [ft^3]
        !                   M_fuel    == Mass of the fuel [lbm]
        !                   Tau       == Reactor constant
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        ! CALL Getq0(Q_Rx,qflux_avg,q0_a,q0_h)
        !
        ! q_vol = 4*qflux_avg/D_o ![Btu/(hr*ft^3)]

        T_fuel       = 1200. ![°F]

        fuel_dens    = 685. ![lbm/ft^3] (FROM OAK RIDGE REPORT ORNL/TM-2000/351 TABLE 3.1)
        fuel_Cp      = Cp_fuel(T_fuel)
        UA           = GetUA(T_fuel)
        V_fuel       = pi*(D_o**2.0)/4.0 * Hfuel * n_frods
        M_fuel       = V_fuel * fuel_dens

        Tau          = (M_fuel * fuel_Cp) / UA

        Tau          = Tau

        HeatTransfer = Q_Rx*EXP(-dt/Tau) + P_fuel*(1 - EXP(-dt/Tau))

    END FUNCTION

    REAL(8) FUNCTION GetUA(T_Rx)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T_Rx
        REAL(8)           :: Q, T_sat

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T_Rx  == Average fuel temperature [F]
        !                   Q     == Nominal reactor power [Btu/hr]
        !                   T_sat == Saturation temperature [°F]
        !**LOGICAL**   ==>
        !----------------------------------!

        Q     = Nom_Power*W2Btu
        T_sat = Tsat(Pref)

        GetUA = Q / (T_Rx-T_sat)

    END FUNCTION

    !Calculates Cp of UO2 fuel (ORNL/TM-2000/351 Equation 4.2)
    REAL(8) FUNCTION Cp_Fuel(T)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T
        REAL(8)           :: T_K, C1, C2, C3, Theta, Ea, K2F, J2Btu, kg2lbm

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T        == Average fuel temperature [°F]
        !                   T_K      == Average fuel temperature [K]
        !                   C1,C2,C3 == Constants
        !                   Theta    == Einstein temperature [K]
        !                   Ea       == Electron activation energy divided by Boltzman constant [K]
        !                   K2F      == Conversion from CHANGE in Kelvin to CHANGE in Fahrenheit
        !                   J2Btu    == Conversion from Joule to Btu
        !                   kg2lbm   == Converstion from kg to lbm
        !**LOGICAL**   ==>
        !----------------------------------!

        !!Conversions from Chart of the Nuclides 17th Edition
        K2F     = 1.8
        J2Btu   = 1./1.05587e3
        kg2lbm  = 1./4.5359237e-1

        T_K     = (T-32) * (5/9) + 273.15 ![°F -> K]

        !!Equation and constants from ORNL/TM-2000/351 Equation and Table 4.2
        C1      = 302.27   ![J/(kg*K)]
        C2      = 8.463e-3 ![J/(kg*K^2)]
        C3      = 8.741e7  ![J/kg]
        Theta   = 548.68   ![K]
        Ea      = 18531.7  ![K]

        Cp_Fuel = C1 * (Theta/T_K)**2.0 * EXP(Theta/T_K) / ((EXP(Theta/T_K) - 1.)**2.) &
                  + 2.0*C2*T_K + C3*Ea*EXP(-Ea/T_K)/(T_K**2.0) ![J/(kg*K)]

        !Convert to [Btu/(lbm*°F)]
        Cp_fuel = Cp_Fuel * J2Btu / (K2F * kg2lbm)

    END FUNCTION

    SUBROUTINE Getq0(Power,qflux_avg,q0_avg,q0_hot)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: Power
        REAL(8),INTENT(OUT) :: qflux_avg,q0_avg, q0_hot

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  Power == Current reactor power [W]
        !                   qflux_avg == Average reactor heat flux [Btu/(hr*ft^2)]
        !                   q0_avg    == q0 term for the average channel [Btu/(hr*ft^2)]
        !                   q0_hot    == q0 term for the hot channel [Btu/(hr*ft^2)]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        qflux_avg = Power*W2Btu*gf/(n_frods*pi*Hfuel*D_o) ![Btu/(hr*ft^2)]
        q0_hot = Fq*qflux_avg                             ![Btu/(hr*ft^2)]
        q0_avg = Fz*qflux_avg                             ![Btu/(hr*ft^2)]

    END SUBROUTINE

    !Calculates outer cladding temperature based on the Thom Correlation (fully-developed nucleate boiling)
    REAL(8) FUNCTION FDNB(z,corr,pressure)
        IMPLICIT NONE

        INTEGER,INTENT(IN):: corr
        REAL(8),INTENT(IN):: z, pressure
        REAL(8):: hc, y
        INTEGER:: m


        !-------Variable Definitions-------!
        !**INTEGER**   ==>  corr     == Correlation to be used: 1 for Thom, 2 for Jens-Lottes
        !                   m        == FDNB Correlation exponent
        !**REAL**      ==>  z        == Axial location
        !                   pressure == Current system pressure [psia]
        !                   hc       == Convective heat transfer coefficient
        !                   y        == FDNB correlation terms
        !**LOGICAL**   ==>
        !----------------------------------!

        IF(corr .NE. 2) THEN !Default is Thom Correlation
            m = 2
            y = EXP(2*pressure/1260) / (72**2)
        ELSE
            m = 4
            y = EXP(4*pressure/900) / (60**4)
        END IF

        !FDNB = Tinf(z) + q(z)/hc

    END FUNCTION

    REAL(8) FUNCTION T_inf(z,h_in,q0,mdot)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: z, h_in, q0, mdot
        REAL(8)           :: h_inf

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  z     == Axial location [z]
        !                   h_in  == Core inlet enthalpy [Btu/lbm]
        !                   q0    == Either average or hot channel q0 [Btu/(hr*ft^2)]
        !                   mdot  == Core mass flow rate [lbm/hr]
        !                   h_inf == Bulk fluid enthalpy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        h_inf = h_in + (pi * D_o) / mdot * q0 * shape_int(z)
        !T_inf =

    END FUNCTION

    REAL(8) FUNCTION T_surface(z,Tco,q0)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: z, Tco, q0
        REAL(8)           :: R_o, R_i

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  z    == Axial location [z]
        !                   Tco  == Cladding temperature at z [°F]
        !                   q0   == Either average or hot channel q0 [Btu/(hr*ft^2)]
        !                   R_o  == Cladding outer radius [ft]
        !                   R_i  == Cladding inner radius [ft]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        R_o = D_o / 2.
        R_i = D_o - 2*t_clad

        T_surface = Tco + q0*shapez(z)*R_o*(1/(H_G*R_i) + 1/(k_clad)*LOG(R_o/R_i))

    END FUNCTION

    ! !Calculates claddding surface temperature
    ! REAL(8) FUNCTION T_co()
    !     IMPLICIT NONE
    !
    !
    !
    !
    ! END FUNCTION

    !Calculates fuel centerline temperature
    REAL(8) FUNCTION T_centerline(z,T_s,q0)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: z, T_s, q0
        REAL(8)           :: T0, tol, correct, error, LHS, RHS
        INTEGER           :: a, loop_s, looplim

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  a       == Counter variable
        !                   loop_s  == Counter variable
        !                   looplim == Maximum number of allowed loops
        !**REAL**      ==>  z       == Axial location [z]
        !                   T_s     == Fuel surface temperature at z [°F]
        !                   q0      == Either average or hot channel q0 [Btu/(hr*ft^2)]
        !                   T0      == Fuel centerline temperature guess
        !                   tol     == Convergence tolerance
        !                   correct == Guess variable correction factor
        !                   error   == Convergence error
        !                   LHS     == Left-hand side of the equation
        !                   RHS     == Right-hand side of the equation
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        tol = 1e-3

        correct = 0.1
        error   = 100
        a       = 0
        loop_s  = 0
        loopLim = 100

        !Fuel centerline temperature guess
        T0 = 1000

        DO WHILE((ABS(error) .GT. tol) .AND. (loop_s .LT. loopLim))
            RHS = q0*shapez(z)*(0.5*D_o) / 2.
            LHS = ( 3978.1*LOG((692.6 + T0)/692.6) + 6.02366e-12/4.*((T0+460.0)**4 - 460.0**4)) - &
                  ( 3978.1*LOG((692.6 + T_s)/692.6) + 6.02366e-12/4.*((T_s+460.0)**4 - 460.0**4))

            error = LHS - RHS
            IF(error .GT. 0) THEN
                T0 = T0 +  1 * correct;
                a  = 0;
            ELSEIF(error .LT. 0) THEN
                T0 = T0 -  1 * correct;
                IF(a .EQ. 0) THEN
                    correct = correct * 0.1;
                    a = a + 1;
                    T0 = T0 + 1 * correct;
                END IF
            ELSE
                EXIT
            END IF

            loop_s = loop_s + 1
        END DO
        IF(loop_s .EQ. loopLim) WRITE(*,*) "CPR Loop Limit Reached!!!!"

        T_centerline = T0

    END FUNCTION

    !Calculates the Nusselt number using the Weisman correlation
    REAL(8) FUNCTION Nu_Weisman(S,D_o)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: S, D_o

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  S   == Rod pitch
        !                   D_o == Cladding outer diameter
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        REAL(8):: C

        C = 0.042*(S/D_o) - 0.024

    END FUNCTION

END MODULE
