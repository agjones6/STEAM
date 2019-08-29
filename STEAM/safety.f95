MODULE safety
! Module is made to do safety Analysis for:
!   Critical Heat Flux (CHF)
!   Fuel Centerline Temperature
!   Max Cladding Temperature

    !Subroutines used
    USE state
    USE globals
    USE pressuredrop
    IMPLICIT NONE

CONTAINS

    SUBROUTINE CHF(G_CORE, PRESS, q0_HOT, U0, CPR, TOLERANCE)
    ! Subroutine is used to calculate the Critical Power Ratio.
    !   The form of CHF being analysed is Dry Out because it is a BWR

        IMPLICIT NONE

        ! Input Variables
        REAL(8), INTENT(IN):: G_CORE, PRESS, q0_HOT, U0
        REAL(8), OPTIONAL, INTENT(IN):: TOLERANCE

        ! Output Variables
        REAL(8), INTENT(OUT):: CPR

        ! Local Variables - Water Proprties
        REAL(8):: Cp, mu, conduct, T_avg, T, T0, Tsat, &
                  h0, m_core, Ho
        ! Local Variables - Geometry
        REAL(8):: zcrit, z_dist, De, Ax
        ! Local Variables - SI Unit Variables for CISE-4 Correlation
        REAL(8):: G_SI, Do_SI, De_SI, H_SI, P_SI, P_c_SI, a_CISE, b_CISE, Xc_SI, qc
        ! Local Varaibles - Conditional criteria
        REAL(8):: criteria_a
        ! Local Variables - Miscellaneous
        REAL(8):: LH, correct, tol, check
        INTEGER:: a, loop_s, loopLim


        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  G_CORE     == Mass Flux of the Core [lbm/hr*ft^2]
        !                   PRESS      == System Pressure [psia]
        !                   q0_HOT     == Max Heat Flux in the Hot Channel [btu/ft^2]
        !                   U0         == Core Inlet Internal Energy
        !                   TOLERANCE  == Optional Tolerance on the convergence. Default is 1e-3
        ! -- Water Props -- Ho         == Non-Boiling Height [ft]
        !                   Cp         == Water Specific Heat
        !                   mu         == Water Viscosity
        !                   conduct    == Water thermal conductivity
        !                   T          == Temperature [°F]
        !                   T0         == Core Inlet Temperature [°F]
        !                   Tsat       == Saturation Temperature [°F]
        !                   T_avg      == Core average Temperature [°F]
        !                   h0         == Inlet Enthalpy [btu/lbm]
        !                   G          == Core Mass Flux [lbm/hr*ft^2]
        ! -- Geometry ----- zcrit      == Critical Boiling length. It is desired to
        !                                   Observe conditions at the top of the fuel so zcrit=hfuel
        !                   De         == Equivalent Diameter [ft]
        !                   Ax         == Cross Sectional Area of a Channel [ft]
        !
        !                   CPR        == Outputted Critical Power Ratio
        ! -- Criteria ----- criteria_a == Conditional criteria for the a value of the CISE-4 correlation
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        ! Geometry Declarations
        Ax = s**2 - (pi/4)*D_o**2 ![ft]
        De = (4*Ax)/(pi*D_o)      ![ft]

        ! Initial Properties
        T0     = Temp(U0)        ![°F]
        Tsat   = Temp(uf(PRESS)) ![°F]
        T_avg  = (T0 + Tsat)/2   ![°F]
        h0     = hliq(T0,PRESS)  ![Btu/lbm]
        m_core = G_CORE * Ax     ![lbm/hr]

        ! Water Properties
        Cp      = Cpl(T_avg, PRESS)       ![Btu / (lbm * °F)]
        conduct = kl(T_avg, PRESS)        ![Btu / (hr * ft * °F)]
        mu      = Viscosity(T_avg, PRESS) ![lbm / (hr * ft)]

        ! Looking Core exit for Dry-Out conditions
        zcrit = Hfuel ![ft]
        Ho    = 0.2   ! Initial Guess [ft]

        ! --> Evaluating the terms in the CISE-4 Correlation
        ! SI Unit Variables for the CISE-4 Correlation
        G_SI   = G_CORE/737.3560421 ! kg/m^2*s
        Do_SI  = D_o/3.2808         ! m
        De_SI  = De/3.2808          ! m
        H_SI   = Hfuel/3.2808       ! m
        P_SI   = PRESS/145.0377439  ! MPa
        P_c_SI = 22.064             ! MPa "Critical Pressure"

        ! Coefficients for the CISE-4 Correlation
        criteria_a = 3375. * (1 - P_SI/P_c_SI)**3
        IF(G_SI .LE. criteria_a) THEN
            a_CISE = 1 / ( 1 + 1.481e-4 * (1-(P_SI / P_c_SI))**(-3) * G_SI)
        ELSE
             a_CISE = (1-(P_SI / P_c_SI)) / ( (G_SI / 1000)**(1/3))
        END IF

        b_CISE = 0.199*((P_c_SI/P_SI) -1)**(0.4) * G_SI * Do_SI**(1.4)

        ! Iterating on the Critical Non-Boiling Height
        IF(PRESENT(TOLERANCE)) THEN
            tol = TOLERANCE
        ELSE
            tol = 1e-3
        END IF

        correct = 0.1
        check   = 100
        a       = 0
        loop_s  = 0
        loopLim = 100

        DO WHILE((ABS(check) .GT. tol) .AND. (loop_s .LT. loopLim))
            ! Critical Quality
            Xc_SI = (Do_SI/De_SI) * ((a_CISE * (H_SI - Ho/3.2808)) &
                /( (H_SI - Ho/3.2808) + b_CISE))

            ! The Check Will Compare
            LH    = (shape_int(zcrit)/shape_int(Ho) - 1)*((hf(PRESS) - h0)/hfg(PRESS))

            check = LH - Xc_SI
            IF(check .GT. 0) THEN
                Ho = Ho +  1 * correct;
                a  = 0;
            ELSEIF(check .LT. 0) THEN
                Ho = Ho -  1 * correct;
                IF(a .EQ. 0) THEN
                    correct = correct * 0.1;
                    a = a + 1;
                    Ho = Ho + 1 * correct;
                ELSE
                END IF
            END IF

            loop_s = loop_s + 1
        END DO

        IF(loop_s .EQ. loopLim) WRITE(*,*) "CPR Loop Limit Reached!!!!"

        qc  = m_core * (hf(PRESS) - h0)/(shape_int(Ho)*pi*D_o)
        CPR = qc/q0_HOT

        IF(CPR/CPR .NE. 1) CPR = 0



    END SUBROUTINE

    REAL(8) FUNCTION T_dist(m,T0,P,q0,Cp,z)
    ! Calculation of Temperature Distribution
    IMPLICIT NONE

    ! Input Variables
    REAL(8), INTENT(IN):: m,T0,P,q0,Cp,z

    !-------Variable Definitions-------!
    !**INTEGER**   ==>
    !**REAL**      ==>  m  == Mass Flow Rate
    !                   T0 == Core Inlet Temperature
    !                   z  == Axial Height
    !                   q0 == Maximum Heat Flux
    !**CHARACTER** ==>
    !**LOGICAL**   ==>
    !----------------------------------!

    T_dist = T0 + ((pi*D_o)/(m*Cp))*q_int(z,q0)

    IF(T_dist .GT. Temp(uf(P))) THEN
        T_dist = Temp(uf(P))
    END IF

    END FUNCTION

    REAL(8) FUNCTION Pc(G,De,TEM,PRESS)
    ! Calculation of the Peclet Number
    IMPLICIT NONE

    ! Input Variables
    REAL(8), INTENT(IN):: G, De, TEM, PRESS

    !-------Variable Definitions-------!
    !**INTEGER**   ==>
    !**REAL**      ==>  G       == Mass Flux
    !                   De      == Equilivent flow diameter
    !                   TEM     == Temperature of fluid
    !                   PRESS   == Current Pressure
    !**CHARACTER** ==>
    !**LOGICAL**   ==>
    !----------------------------------!
        ! Calculation
        Pc = Re(G,De,TEM,PRESS) * Pr(TEM,PRESS)

    END FUNCTION

    !Calculation of the Prandlt Number
    REAL(8) FUNCTION Pr(T,P)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P
        REAL(8)           :: Cp, kw, mu

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T == Temperature [°F]
        !                   P == Pressure    [°F]
        !                   Cp == Specific heat of water ![Btu/(lbm°F)]
        !                   kw == Thermal conductivity of water ![Btu/(hr*ft*°F)]
        !                   mu == Dynamic viscosity of water ![lbm/(hr*ft)]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!


        Pr = Cpl(T,P) * Viscosity(T,P) / kl(T,P)

        END FUNCTION

END MODULE
