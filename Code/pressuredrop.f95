! PRESSURE DROP CALCULATION MODULE

MODULE pressuredrop

    USE globals
    USE state
    USE conditional
    IMPLICIT NONE

CONTAINS

    !Reynolds Number
    REAL(8) FUNCTION Re(G,De,T,P)
        IMPLICIT NONE

        !Inputs
        REAL(8):: G, De, T, P
        !Calculations
        REAL(8):: mu_l
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  G    == Mass Flux [lbm/hr]
        !                   De   == Equillivent Diameter
        !                   rho  == Density
        !                   mu_l == Water specific viscosity [lbm/hr*ft]
        !                   T    == Temperature
        !                   P    == Pressure
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        mu_l = Viscosity(T,P)
        Re = (G*De)/mu_l
        
    END FUNCTION

    !Two phase multiplier for friction
    REAL(8) FUNCTION Philo2(X,T,P)
        IMPLICIT NONE

        !Inputs
        REAL(8):: X, T, P
        !Calculations
        REAL(8):: Ki, mu_l, mu_g, rho_g, rho_f
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  X       == Quality
        !                   T       == Temperature
        !                   P       == Pressure
        !                   Ki      == Variable for calculation
        !                   mu_l    == Viscosity of liquid [lbm/hr*ft]
        !                   mu_g    == Viscosity of gas [lbm/hr*ft]
        !                   rho_g   == Saturated gas density
        !                   rho_f   == Saturated liquid density
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        mu_l  = Viscosity(T,P)
        mu_g  = Mug(P)
        rho_g = rhog(P)
        rho_f = rhof(P)

        IF(X .GT. 0) THEN !Two phase condition
            Ki =  (mu_l/mu_g)**(0.2) * ((1-X)/X)**1.8 * (rho_g/rho_f);

            Philo2 = (1 + 20/((Ki)**0.5) + 1/Ki) * (1-X)**1.8;
        ELSE !Single phase condition
            Philo2 = 1
        END IF


    END FUNCTION

    !Two phase multiplier for local losses
    REAL(8) FUNCTION Psi_loss(X,P)
        IMPLICIT NONE

        !Inputs
        REAL(8):: X, P
        !Calculations
        REAL(8):: rho_g, rho_f

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  X       == Quality
        !                   T       == Temperature
        !                   P       == Pressure
        !                   rho_g   == Saturated gas density
        !                   rho_f   == Saturated liquid density
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        rho_g = rhog(P)
        rho_f = rhof(P)

        IF(X .GT. 0) THEN
            Psi_loss = 1 + (rho_f/rho_g - 1) * X
        ELSE
            Psi_loss = 1
        END IF
        ! IF(Psi_loss .GT. 5) Psi_loss = 5

    END FUNCTION

    !Friction pressure loss
    !Note: Only multiplied by ONE veloctiy in the function. This needs to be
    !      multiplied by a velocity outside the function
    REAL(8) FUNCTION DP_friction(vel2,rho1,rho2,X,T,P,De,L) !psi
        IMPLICIT NONE

        !Inputs
        REAL(8):: vel2, De, rho1(1,2), rho2, X, T(1,2), P, L
        !Calculations
        REAL(8):: f,G,Rey
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  vel2  == Velocity of node [ft/hr]
        !                   rho1 == Density of node (rho_n)[lbm/ft^3]
        !                   rho2 == Density of subsequent junction (rho_jn) [lbm/ft^3]
        !                   X    == Quality
        !                   T    == Temperature [°F]
        !                   P    == Pressure [psia]
        !                   De   == Equillivent Diameter [ft]
        !                   L    == Length [ft]
        !                   f    == Friction factor
        !                   G    == Mass Flux
        !
        !                   gc   == Global variable for conversion
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        G    = vel2 * rho1(1,2)
        Rey  = Re(ABS(G),De,T(1,2),P)

        !Determining Frition Factor
        IF(Rey .LE. 0) THEN !To Prevent Dividing by 0
            f = 0.
        ELSE !Smooth Friction Factor
            f = 0.184*Rey**(-0.2)
        END IF
        !Ensuring that the friction factor isn't incredibly large
        ! IF(f .GT. 0.1) THEN
        !     f = 0.1
        ! END IF

        !Friction Pressure Drop
        DP_friction = (f * L / De) * (rho2**2/rho1(1,1)) * ABS(vel2) * Psi_loss(X,P)


    END FUNCTION

    !Local pressure loss
    !Note: Only multiplied by ONE veloctiy in the function. This needs to be
    !      multiplied by a velocity outside the function
    REAL(8) FUNCTION DP_loss(vel,rho1,rho2,X,T,P,K) !psi
        IMPLICIT NONE

        !Inputs
        REAL(8):: vel, rho1, rho2, X, T, P, K
        !Calculations
        !REAL(8)::
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  vel  == Velocity [ft/hr]
        !                   rho1 == Density of node (rho_n)[lbm/ft^3]
        !                   rho2 == Density of subsequent junction (rho_jn) [lbm/ft^3]
        !                   X    == Quality
        !                   T    == Temperature [°F]
        !                   P    == Pressure [psia]
        !                   K    == Loss Coefficient
        !
        !                   gc   == Global variable for conversion
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        DP_loss = K * (rho2**2/rho1) * ABS(vel) * Psi_loss(X,P)

        !DP_loss = (K * G**2 * Psi_loss(X,P)) / (2 * rho * gc * 144)\

    END FUNCTION

    !Gravitational pressure Loss
    !Note: does not get multiplied by a velocity term
    REAL(8) FUNCTION DP_elev(rho2,height) !psi
        IMPLICIT NONE

        !Inputs
        REAL(8):: rho2, height
        !Calculations
        !REAL(8)::z
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  rho2   == Density entering region (n)[lbm/ft^3]
        !                   height == Height of region
        !
        !                   gc   == Global variable for conversion
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        DP_elev = (rho2 * height) * (grav)

    END FUNCTION

    !Phase Velocity pressure loss
    !Note: Includes both velocity terms
    REAL(8) FUNCTION DP_vel(alpha1,alpha2,rho1,rho2,vg1,vg2,vl1,vl2,T1,T2,Psys) !psi
        !Inputs
        REAL(8):: alpha1, alpha2, rho1, rho2, vg1, vg2, vl1, vl2, T1, T2, Psys
        !Calculations
        REAL(8):: rho_l1, rho_l2, rho_g1, rho_g2
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  alpha1 == Vapor void fraction entering the region (jn-1)
        !                   alpha2 == Vapor void fraction leaving the region (jn)
        !                   rho1   == Density entering region [lbm/ft^3]
        !                   rho2   == Density leaving region [lbm/ft^3]
        !                   vg1    == Velocity of gas entering node (jn-1) [ft/hr]
        !                   vg2    == Velocity of gas leaving node (jn) [ft/hr]
        !                   vl1    == Velocity of liquid entering node (jn-1) [ft/hr]
        !                   vl2    == Velocity of liquid leaving node (jn) [ft/hr]
        !                   T1     == Temperature entering node (jn-1) [°F]
        !                   T2     == Temperature leaving node (jn) [°F]
        !                   P1     == Pressure entering node (jn-1) [°F]
        !                   P2     == Pressure leaving node (jn) [°F]
        !                   rho_l1 == Density of liquid entering node (jn-1) [lbm/ft^3]
        !                   rho_l2 == Density of liquid leaving node (jn) [lbm/ft^3]
        !
        !                   gc    == Global variable for conversion
        !**CHARACTER** ==>

        !**LOGICAL**   ==>
        !----------------------------------!

        !Properties of fluid entering
        rho_l1 = rhol(uliq(T1),Psys)
        rho_g1 = rhog(Psys)
        !Properties of fluid exiting
        rho_l2 = rhol(uliq(T2),Psys)
        rho_g2 = rhog(Psys)

        !Calculation of Pressure change
        DP_vel = ( ((1-alpha2)*rho_l2*alpha2*rho_g2)/rho2 ) * (vg2 - vl2)**2 &
               - ( ((1-alpha1)*rho_l1*alpha1*rho_g1)/rho1 ) * (vg1 - vl1)**2

    END FUNCTION

    !Acceleration pressure loss
    !Note: Includes both velocity terms
    REAL(8) FUNCTION DP_acc(vel2,vel1,rho2) !psi
        IMPLICIT NONE

        !Inputs
        REAL(8):: vel2, vel1, rho2
        !Calculations
        !REAL(8)::
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  vel1 == Velocity entering region (jn-1) [ft/hr]
        !                   vel2 == Velocity leaving region (jn) [ft/hr]
        !                   rho1  == Density entering region (jn-1) [lbm/ft^3]
        !                   gc    == Global variable for conversion
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        DP_acc = rho2 * vel2 * (vel2 - vel1)
        !DP_acc = (G2**2 / rho2 - G1**2 / rho1) / (gc * 144)

    END FUNCTION

    !Pressure drop with time
    !Note: This needs to multiplied by a velocity outside the function!!
    REAL(8) FUNCTION DP_time(rho2,L,dt) !psi
        IMPLICIT NONE

        !Inputs
        REAL(8):: rho2, L, dt
        !Calculations
        !REAL(8)::
        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  rho2 == Density of subsequent junction (rho_jn) [lbm/ft^3]
        !                   L    == Length of node (dz_n) [ft]
        !                   dt   == Time Step (hr?)
        !
        !                   gc   == Global variable for conversion
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        DP_time = L * rho2 / (dt)

    END FUNCTION

END MODULE
