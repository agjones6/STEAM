! Conditional Equations Module
MODULE conditional

    USE globals
    USE state
    USE steamdome
    USE files
    IMPLICIT NONE

CONTAINS

    !Finds the new density based on number of phases present
    REAL(8) FUNCTION rho(test,u,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: test, u, P, ag
        ! REAL(8),INTENT(OUT)::
        REAL(8):: crit

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  test == Phase criteria: density
        !                   T    == Current temperature
        !                   P    == Current pressure
        !                   ag   == Void fraction
        !                   u    == Specific Internal Energy
        !                   crit == Criteria for checking single phase or 2 phase
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        !Checking Phase
        crit = rhof(P)

        IF(test .GE. crit) THEN
            !Single phase
            rho = rhol(u,P)                     ![lbm/ft^3]
        ELSE
            !Two-phase
            rho = (1 - ag)*rhof(P) + ag*rhog(P) ![lbm/ft^3]
        END IF

    END FUNCTION


    !Finds the density*specifc internal energy based on number of phases present
    REAL(8) FUNCTION rho_u(test,u,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN)::  test, u, P, ag
        ! REAL(8),INTENT(OUT):: y
        REAL(8):: crit

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  test == Phase criteria: density
        !                   T    == Current temperature
        !                   P    == Current pressure
        !                   ag   == Void fraction
        !                   crit == Single- / two-phase criteria
        !                   u    == Specific internal energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        crit = rhof(P)

        IF(test .GE. crit) THEN
            rho_u  = rhol(u,P) * u                               ![Btu/ft^3]
        ELSE
            rho_u  = (1 - ag) * rhof(P)*uf(P) + ag*rhog(P)*ug(P) ![Btu/ft^3]
        END IF

    END FUNCTION


    !Partial derivative of density with respect to either specifc internal
    !energy or void fraction, depending on number of phases present
    SUBROUTINE drhodpsi(test,y,T,P)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: test, T, P
        REAL(8),INTENT(OUT) :: y
        REAL(8)             :: crit, u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  test == Current value of single- /
        !                           two-phase criteria
        !                   y    == drhodpsi output
        !                   T    == Current temperature
        !                   P    == Current pressure
        !                   crit == Single- / two-phase criteria
        !                   u    == Specific internal energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        crit = rhof(P)

        IF(test .GE. crit) THEN
            u = uliq(T)
            y = drholdu(u,P)         ![lbm^2/(Btu*ft^3)]
        ELSE
            y = -(rhof(P) - rhog(P)) ![lbm/ft^3]
        END IF

    END SUBROUTINE


    !Partial derivative of density*specific internal energy with respect to
    !either specifc internal energy or void fraction, depending on number of
    !phases present
    SUBROUTINE drho_udpsi(test,y,T,P)
        IMPLICIT NONE

        REAL(8),INTENT(IN)::  test, T, P
        REAL(8),INTENT(OUT):: y
        REAL(8):: crit, u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  test == Current value of single- /
        !                           two-phase criteria
        !                   y    == drhodpsi output
        !                   T    == Current temperature
        !                   P    == Current pressure
        !                   crit == Single- / two-phase criteria
        !                   u    == Specific internal energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        crit = rhof(P)

        IF(test .GE. crit) THEN
            u = uliq(T)
            y = rhol(u,P) + u * drholdu(u,P)     ![lbm/ft^3]
        ELSE
            y = -(rhof(P)*uf(P) - rhog(P)*ug(P)) ![Btu/ft^3]
        END IF

    END SUBROUTINE

    !Partial derivative of density with respect to Pressure evaluated at either
    !internal energy at the given node or void fraction at the node
    SUBROUTINE drhodP(test,y,T,P,alpha)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  ::  test, T, P, alpha
        REAL(8),INTENT(OUT) :: y
        REAL(8)             :: crit, u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  test  == Current value of single- /
        !                            two-phase criteria
        !                   y     == drhodpsi output
        !                   T     == Current temperature
        !                   P     == Current pressure
        !                   crit  == Single- / two-phase criteria
        !                   alpha == vapor volume void fraction
        !                   u    == Specific internal energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        crit = rhof(P)

        IF(test .GE. crit) THEN
            u = uliq(T)
            y = drholdP(u,P)                            ![lbm/(ft^3*psia)]
        ELSE
            y = (1-alpha)*drhofdP(P) + alpha*drhogdP(P) ![lbm/ft^3*psia)]
        END IF

    END SUBROUTINE

    !Partial derivative of density*specific internal energy with respect to
    !Pressure evaluated at either internal energy at the given node or void
    !fraction at the node
    SUBROUTINE drho_udP(test,y,T,P,alpha)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  ::  test, T, P, alpha
        REAL(8),INTENT(OUT) :: y
        REAL(8)             :: crit, u_f, u_g, u_l

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  test  == Current value of single- /
        !                            two-phase criteria
        !                   y     == drhodpsi output
        !                   T     == Current temperature
        !                   P     == Current pressure
        !                   crit  == Single- / two-phase criteria
        !                   alpha == vapor volume void fraction
        !                   u_l   == Specific internal energy of liquid
        !                   u_f   == Specific internal energy of saturated liquid
        !                   u_g   == Specific internal energy of saturated gas
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        crit = rhof(P)
        u_l = uliq(T)
        u_f = uf(P)
        u_g = ug(P)

        IF(test .LE. crit) THEN
            y = drholdP(u_l,P) * u_l                      ![dimensionless]
        ELSE
            y = (1-alpha)*(u_f*drhofdP(P) + rhof(P)*dufdP(P)) + &
                alpha*(u_g*drhogdP(P) + rhog(P)*dugdP(P)) ![dimensionless]


        END IF


    END SUBROUTINE


    REAL(8) FUNCTION S_n(alpha,rho,G,X,T,P,Ax,j1,j2) !psi

        !Inputs
        INTEGER,INTENT(IN) :: j1, j2
        REAL(8),INTENT(IN) :: alpha(:,:), rho(:,:), G(:,:), X(:,:), T(:,:), P ,Ax(:)
        !Calculations
        REAL(8)            :: rho_l1, rho_l2, rho_g1, rho_g2, h_fg1, h_fg2, vg1, vg2, vl1, vl2
        !-------Variable Definitions-------!
        !**INTEGER**   ==>  j1     == Number of previous junction
        !                   j2     == Number of current junction
        !**REAL**      ==>  alpha  == Vapor void fraction matrix
        !                   rho    == Density  [lbm/ft^3]
        !                   G      == Mass Flux [lbm/ft^2*hr]
        !                   X      == Quality
        !                   T      == Temperature  [°F]
        !                   P      == Current system pressure [psia]
        !                   Ax     == Cross sectional area  [ft^2]
        !                   rho_l1 == Density of liquid entering node (jn-1) [lbm/ft^3]
        !                   rho_l2 == Density of liquid leaving node (jn) [lbm/ft^3]
        !                   h_fg1  == Saturated gas-liquid (jn-1) [btu/lbm]
        !                   h_fg2  == Saturated gas-liquid (jn-1) [btu/lbm]
        !                   vg1    == Velocity of gas entering node (jn-1) [ft/hr]
        !                   vg2    == Velocity of gas leaving node (jn) [ft/hr]
        !                   vl1    == Velocity of liquid entering node (jn-1) [ft/hr]
        !                   vl2    == Velocity of liquid leaving node (jn) [ft/hr]
        !                   n      == Size of the matrices (dim 1)
        !                   gc     == Global variable for conversion
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        !Properties of fluid entering
        rho_l1 = rhol(uliq(T(j1,1)),P)
        rho_g1 = rhog(P)
        h_fg1  = hfg(P)
        CALL phase_vel(vg1, vl1, G(j1,1), X(j1,1), alpha(j1,1), P, j1)

        !Properties of fluid exiting
        rho_l2 = rhol(uliq(T(j2,1)),P)
        rho_g2 = rhog(P)
        h_fg2  = hfg(P)
        CALL phase_vel(vg2, vl2, G(j2,1), X(j2,1), alpha(j2,1), P, j2)

        !Calculation of Pressure change
        S_n = ( ((1-alpha(j2,1))*rho_l2*alpha(j2,1)*rho_g2)/rho(j2,1) ) * h_fg2 * &
                 (vg2 - vl2) * Ax(j2) - ( ((1-alpha(j1,1))*rho_l1*alpha(j1,1)*rho_g1)/rho(j1,1) ) &
                 * h_fg1 * (vg1 - vl1) * Ax(j1) ![Btu/hr]

        IF(j2 .EQ. 11) THEN
            ! WRITE(*,*) alpha(j1,1), rho_l1, rho_g1, rho(j1,1), h_fg1
            S_n = ( ((1-alpha(j1,1))*rho_l1*alpha(j1,1)*rho_g1)/rho(j1,1) ) &
                     * h_fg1 * (vg2 - vl2) * Ax(j1) ![Btu/hr]
        END IF

    END FUNCTION

    SUBROUTINE phase_vel(vg,vl,G,X,alpha,P,junction)
              !phase_vel(vg1,vl1,Gj(n,1),Xj(n-1,1),agj(n-1,1),Psys,n-1)

        !Inputs
        REAL(8),INTENT(IN):: G, X, alpha, P
        INTEGER,INTENT(IN):: junction
        REAL(8),INTENT(INOUT):: vg, vl
        !Calculation
        !REAL(8)::

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  junction == Current junction index
        !**REAL**      ==>  G        == Mass Flux matrix [lbm/hr*ft^2]
        !                   X        == Quality matrix
        !                   alpha    == Vapor void fraction matrix
        !                   P        == Pressure [psia]
        !                   vg       == Gas phase velocity [ft/hr]
        !                   vl       == liquid phase velociy [ft/hr]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------
        IF( (junction .EQ. 9) .OR. (alpha .LE. 0) ) THEN
            vl = 0
            vg = 0
        ! ELSE IF(alpha .LE. 0) THEN
        !     vl = G / rhof(P)
        !     vg = 0
        ELSE
            vl = ( G*(1-X)) / ( rhof(P)*(1-alpha) ) ![ft/hr]
            vg = ( G*X ) / ( rhog(P)*alpha )        ![ft/hr]
        END IF

    END SUBROUTINE

    !Determines all of the Temperature, Enthalpy, Void Fraction,
    ! Liquid Internal Energy, and Densities values for each junction and node
    SUBROUTINE InitialProps(v_j,T_in,Psys,al_rhol11,Tn,Tj,agn,agj,rhon,rhoj,uln,ulj, &
                            rho_un,rho_uj,Xn,Xj,Gj,q0,dt,el_time,tcvPos,mfeed,UNIT)

        !Inputs
        REAL(8),INTENT(INOUT) :: T_in, Psys(:,:), q0, dt, el_time, tcvPos, mfeed
        REAL(8),INTENT(INOUT) :: Tn(:,:),     Tj(:,:), &
                                 agn(:,:),    agj(:,:), &
                                 rhon(:,:),   rhoj(:,:), &
                                 uln(:,:),    ulj(:,:), &
                                 rho_un(:,:), rho_uj(:,:), &
                                 Xn(:,:),     Xj(:,:), &
                                 Gj(:,:),     v_j(:,:)
        REAL(8),INTENT(INOUT) :: al_rhol11(:,:)
        INTEGER,OPTIONAL,INTENT(IN) :: UNIT
        !Calculation
        INTEGER:: n, my_unit
        ! REAL(8)::

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  Tn(),Tj()          == Junction and Node Temperatures [°F]
        !                   agn(),agj()        == Junction and Node Void Fractions
        !                   rhon(),rhoj()      == Junction and Node Densities [lbm/ft^3]
        !                   uln(),ulj()        == Junction and Node Internal Energies
        !                   rho_un(),rho_uj()  == Junction and Node Density * Internal Energy
        !                   T_in               == Inlet Temperature [°F]
        !                   Psys()             == System Pressure [psia]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        IF(PRESENT(UNIT) .EQV. .FALSE.) THEN

            !---> Nodal Properties
            ! Normal Conditions
            Tn             = T_in
            uln            = uliq(T_in)
            rhon(1:10,:)   = rhol(uliq(T_in),Psys(1,1))
            rho_un(1:10,:) = rhol(uliq(T_in),Psys(1,1)) * uliq(T_in)
            Xn(1:10,:)     = 0
            agn(1:10,:)    = 0

            ! SteamDome
            agn(11,:)    = Init_Void_SD
            rhon(11,:)   = SD_rho(T_in,Psys(1,1),agn(11,1))
            Xn(11,:)     = (rhof(Psys(1,1)) - rhon(11,1)) / &
                           (rhof(Psys(1,1)) - rhog(Psys(1,1)))
            rho_un(11,:) = SD_rho_u(T_in,Psys(1,1),agn(11,1))
            ! al_rhol11
            al_rhol11(1,1) = (1-agn(11,1)) * rhol(uln(11,1),Psys(1,1))
            al_rhol11(1,2) = (1-agn(11,2)) * rhol(uln(11,2),Psys(1,2))

            ! Junction Properties
            Tj             = T_in
            ulj            = uliq(T_in)
            rhoj(1:11,:)   = rhol(uliq(T_in),Psys(1,1))
            rho_uj(1:11,:) = rhol(uliq(T_in),Psys(1,1)) * uliq(T_in)
            Xj(1:11,:)     = 0
            agj(1:11,:)    = 0
            Gj             = 0
            v_j            = 0

            !Setting Junction 12
            agj(12,:)    = 1
            Xj(12,:)     = 1
            rhoj(12,:)   = rhog(Psys(1,1))
            rho_uj(12,:) = agj(12,1)*rhog(Psys(1,1))*ug(Psys(1,1))
            ulj(12,:)    = 0




        ELSE
            CALL ReadSnap(VEL = v_j, DENS_n = rhon, DENS_j = rhoj, DENSu_n = rho_un, &
                          DENSu_j = rho_uj, VOID_n = agn, VOID_j = agj, intENERGY_n = uln, &
                          intENERGY_j = ulj, QUAL_n = Xn, QUAL_j = Xj, TEM_n = Tn, TEM_j = Tj, &
                          PRESS = Psys, VOID_DENS11 = al_rhol11, HEATq = q0, &
                          TIME = el_time, VALVE = tcvPos, DELTA_t = dt, FEEDFLOW = mfeed, &
                          FILE_NUM = UNIT)

            v_j(:,2) = v_j(:,1)
            rhon(:,2) = rhon(:,1)
            rhoj(:,2) = rhoj(:,1)
            rho_un(:,2) = rho_un(:,1)
            rho_uj(:,2) = rho_uj(:,1)
            agn(:,2) = agn(:,1)
            agj(:,2) = agj(:,1)
            uln(:,2) = uln(:,1)
            ulj(:,2) = ulj(:,1)
            Xn(:,2) = Xn(:,1)
            Xj(:,2) = Xj(:,1)
            Tn(:,2) = Tn(:,1)
            Tj(:,2) = Tj(:,1)

            !Mass Flux
            DO n = 1,SIZE(v_j,1)
                Gj(n,2) = v_j(n,2) * rhoj(n,2)
            END DO
            Gj(:,1) = Gj(:,2)

        END IF

    END SUBROUTINE

    !Determines all of the Temperature, Enthalpy, Void Fraction,
    ! Liquid Internal Energy, and Densities values for each junction and node
    SUBROUTINE UpdateProps_k(dP,Psys,dPsi_n,d_ul11,al_rhol11,Tn,Tj,agn,agj,rhon,rhoj,uln,ulj, &
                            rho_un,rho_uj,Xn,Xj,Gj,v_j,qn,Vn,dt,dt_change,tcvPos)
        !Inputs
        REAL(8),INTENT(IN):: dP, qn(:,:), Vn(:), dt, d_ul11, tcvPos
        REAL(8),INTENT(INOUT):: Psys(:,:), dPsi_n(:), al_rhol11(:,:)
        REAL(8),INTENT(INOUT)::     Tn(:,:),     Tj(:,:), &
                                   agn(:,:),    agj(:,:), &
                                  rhon(:,:),   rhoj(:,:), &
                                   uln(:,:),    ulj(:,:), &
                                rho_un(:,:), rho_uj(:,:), &
                                    Xn(:,:),     Xj(:,:), &
                                                  Gj(:,:), &
                                                 v_j(:,:)
        LOGICAL,INTENT(IN):: dt_change
        !Calculation
        INTEGER:: n, s1, s2
        REAL(8):: y, w, rho_check, rhoj_dum, uln_dum

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n                  == Counter variable
        !                   s1                 == number of nodes
        !                   s2                 == number of junctions
        !**REAL**      ==>  dt                 == Current time step
        !                   d_ul11             == Change in specific interanl energy in node 11
        !                   dP                 == Change in pressure [psia]
        !                   Tn(),Tj()          == Junction and Node Temperatures [°F]
        !                   agn(),agj()        == Junction and Node Void Fractions
        !                   rhon(),rhoj()      == Junction and Node Densities [lbm/ft^3]
        !                   uln(),ulj()        == Junction and Node Internal Energies
        !                   rho_un(),rho_uj()  == Junction and Node Density * Internal Energy
        !                   T_in               == Inlet Temperature [°F]
        !                   Psys               == System Pressure [psia]
        !                   dPsi_n             == Change in internal energy or void fraction
        !                   qn()               == Values of heat generation rate [Btu/lbm]
        !                   dt_change          == if dt changed due to convergence
        !                   tcvPos             == Turbine control valve position
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        s1 = SIZE(Tn,1)
        s2 = SIZE(Tj,1)

        !If dt changed, the k value is set to the t value
        IF(dt_change)THEN
            !Setting Nodes
            DO n = 1,s1
                Tn(n,2)     = Tn(n,1)
                agn(n,2)    = agn(n,1)
                rhon(n,2)   = rhon(n,1)
                uln(n,2)    = uln(n,1)
                rho_un(n,2) = rho_un(n,1)
                Xn(n,2)     = Xn(n,1)
            END DO

            !Setting Junctions
            DO n = 1,s2
                Tj(n,2)     = Tj(n,1)
                agj(n,2)    = agj(n,1)
                rhoj(n,2)   = rhoj(n,1)
                ulj(n,2)    = ulj(n,1)
                rho_uj(n,2) = rho_uj(n,1)
                Xj(n,2)     = Xj(n,1)
                Gj(n,2)     = Gj(n,1)
                v_j(n,2)    = v_j(n,1)
            END DO

            !Setting pressure and alpha*rho11 back to the previous t value
            al_rhol11(1,2)  = al_rhol11(1,1)
            Psys(1,2)       = Psys(1,1)
        ELSE !When dt doesnt change and the time step is kept constant

            !--> Updating System Pressure
            ! Psys(1,2) = Psys(1,2) + dP/144 ![lbf / in^2]

            !---> Updating Each Property at the node and junction
            DO n = 1,s1

                rho_check = rhon(n,2)
                rhoj_dum  = rhoj(n,2)
                uln_dum   = uln(n,2)

                IF((n .NE. 11)) THEN ! .OR. (tcvPos .LE. 0)) THEN
                    !rho * u
                    CALL drho_udpsi(rho_check,y,Tn(n,2),Psys(1,2))
                    CALL   drho_udP(rho_check,w,Tn(n,2),Psys(1,2),agn(n,2))
                    rho_un(n,2) = rho_u(rho_check, uln(n,2), Psys(1,2), agn(n,2)) &
                                + dPsi_n(n) * y + (dP/144) * w

                    !Defining a new density for the phase criteria
                    CALL drhodpsi(rho_check,y,Tn(n,2),Psys(1,2))
                    CALL   drhodP(rho_check,w,Tn(n,2),Psys(1,2),agn(n,2))
                    rhon(n,2)   = rho(rho_check, uln(n,2), Psys(1,2), agn(n,2)) &
                                  + dPsi_n(n) * y + (dP/144) * w
                    ! IF(n .EQ. 11) al_rhol11(1,2) = SD_al_rhol(Tn(11,1),Psys(1,2),agn(11,1))

                    rhoj(n,2)   = rhon(n,2)

                ELSE !Steam Dome (n = 11)
                    IF(agn(n,2) .LE. 0.) dPsi_n(n) = 0
                    rhon(n,2)      = SD_rho(Tn(n,2),Psys(1,2),agn(n,2)) &
                                   + dPsi_n(n) * SD_drhodag(Tn(n,2),Psys(1,2)) &
                                   + d_ul11 * SD_drhodul(Tn(n,2),Psys(1,2),agn(n,2)) &
                                   + (dP/144) * SD_drhodP(Tn(n,2),Psys(1,2),agn(n,2))
                    !rho * u
                    rho_un(n,2)    = SD_rho_u(Tn(n,2),Psys(1,2),agn(n,2)) &
                                   + dPsi_n(n) * SD_drho_udag(Tn(n,2),Psys(1,2)) &
                                   + d_ul11 * SD_drho_udul(Tn(n,2),Psys(1,2),agn(n,2)) &
                                   + (dP/144) * SD_drho_udP(Tn(n,2),Psys(1,2),agn(n,2))
                    !alphal * rhol
                    al_rhol11(1,2) = SD_al_rhol(Tn(n,2),Psys(1,2),agn(n,2)) &
                                   + dPsi_n(n) * SD_dalrholdag(Tn(n,2),Psys(1,2)) &
                                   + d_ul11 * SD_dalrholdul(Tn(n,2),Psys(1,2),agn(n,2)) &
                                   + (dP/144) * SD_dalrholdP(Tn(n,2),Psys(1,2),agn(n,2))

                    rhoj(11,2)      = rhol(uln(11,2),Psys(1,2))
                END IF

                !Updating Internal Energy or Void Fraction
                IF((n .NE. 11)) THEN ! .OR. (tcvPos .LE. 0)) THEN
                    IF(rho_check .GE. rhof(Psys(1,2))) THEN !single Phase
                        ! IF(n .EQ. 8) WRITE(*,*) "Single"
                        uln(n,2) = uln(n,2) + dPsi_n(n) !
                        ! agn(n,2) = 0
                    ELSE !two phase condition
                        ! IF(n .EQ. 8) WRITE(*,*) "Two <<<-----"
                        agn(n,2) = agn(n,2) + dPsi_n(n)
                        uln(n,2) = uf(Psys(1,2))
                    END IF
                ELSE !Case for the 11th node (steam dome)
                    uln(n,2) = uln(n,2) + d_ul11
                    agn(n,2) = agn(n,2) + dPsi_n(n)
                END IF

                !--> Error Checks
                !Ensuring there is no negative void fraction
                IF(agn(n,2) .LE. 0.) agn(n,2) = 0
                IF(agn(n,2) .GT. 1.) agn(n,2) = 1
                ! Ensuring the Internal energy can't be greater than saturation
                IF(uln(n,2) .GE. uf(Psys(1,2))) uln(n,2) = uf(Psys(1,2))
                !Ensuring Densities dont get out of control
                IF(rhon(n,2) .LT. rhog(Psys(1,2))) THEN
                    rhon(n,2)   = rhog(Psys(1,2))
                    rho_un(n,2) = agn(n,2)*rhog(Psys(1,2))*ug(Psys(1,2))
                END IF
                ! IF(rhon(n,2) .GT. rhog(Psys(1,2))) THEN
                !     rhon(n,2)   = rhog(Psys(1,2))
                !     rho_un(n,2) = agn(n,2)*rhog(Psys(1,2))*ug(Psys(1,2))
                ! END IF
                IF((rho_check .LT. rhof(Psys(1,2))) .AND. (rhon(n,2) .GT. rhof(Psys(1,2)))) THEN
                    rhon(n,2) = rhof(Psys(1,2))
                END IF

                !Bulk Fluid Temperature
                Tn(n,2)     = bulktemp(Psys(1,2),uln(n,2))

                !Quality
                IF(n .EQ. 9) THEN !Seperator
                    Xn(n,2) = agn(n,2)/(agn(n,2) + (1-agn(n,2))*(rhof(Psys(1,2)) &
                              /rhog(Psys(1,2))))
                ELSE IF(n .EQ. 11) THEN !Steam Dome
                    Xn(n,2) = (agn(n,2)/(agn(n,2) + (1-agn(n,2))*(rhof(Psys(1,2)) &
                              /rhog(Psys(1,2)))))
                ELSE
                    Xn(n,2) = X(Psys(1,2),agn(n,2),Gj(n,2))
                END IF

                !Mass Flux
                Gj(n,2) = v_j(n,2) * rhoj(n,2)

                ! Updating junctions
                IF(n .NE. 11) THEN
                    agj(n,2)    = agn(n,2)
                    rhoj(n,2)   = rhon(n,2)
                    ulj(n,2)    = uln(n,2)
                    rho_uj(n,2) = rho_un(n,2)
                    Xj(n,2)     = Xn(n,2)
                    Tj(n,2)     = bulktemp(Psys(1,2),ulj(n,2)) !Tn(n,2)
                ELSE
                    Tj(11,2)     = Tn(n,2)
                    agj(11,2)    = 0
                    rhoj(11,2)   = rhol(uln(11,2),Psys(1,2))
                    ulj(11,2)    = uln(n,2)
                    rho_uj(11,2) = rhoj(11,2)*ulj(11,2)
                    Xj(11,2)     = 0
                END IF

                !Checking for NANs
                IF((Tn(n,2) .NE. Tn(n,2)) .OR. &
                        (agn(n,2)    .NE. agn(n,2)) .OR. &
                        (rhon(n,2)   .NE. rhon(n,2)) .OR. &
                        (uln(n,2)    .NE. uln(n,2)) .OR. &
                        (rho_un(n,2) .NE. rho_un(n,2)) .OR. &
                        (Xn(n,2)     .NE. Xn(n,2)) .OR. &
                        (Gj(n,2)     .NE. Gj(n,2)) ) THEN
                   WRITE(*,*)
                   WRITE(*,*) "There Are NaNs. Calculations have been stopped."
                   WRITE(*,*)
                   STOP
                END IF

            END DO


            !Junction 12: Steam Outlet
            agj(12,2)    = 1
            Tj(12,2)     = Tj(9,2)
            rhoj(12,2)   = rhog(Psys(1,2))
            ulj(12,2)    = uliq(Tj(12,2))
            rho_uj(12,2) = rho_u(rhoj(12,2), Tj(12,2), Psys(1,2), agj(12,2))
            Gj(12,2)     = v_j(12,2) * rhoj(12,2)
            Xj(12,2)     = X(Psys(1,2),agj(12,2),Gj(12,2))

            !--> Updating System Pressure
            Psys(1,2) = Psys(1,2) + dP/144 ![lbf / in^2]

        END IF
    END SUBROUTINE


    !Sets time t Variales equal to the past time k values
    SUBROUTINE UpdateProps_t(Psys,al_rhol11,Tn,Tj,agn,agj,rhon,rhoj,uln,ulj, &
                            rho_un,rho_uj,Xn,Xj,Gj,v_j,Pprev)
        !Inputs
        REAL(8),INTENT(INOUT)::     Tn(:,:),     Tj(:,:), &
                                   agn(:,:),    agj(:,:), &
                                  rhon(:,:),   rhoj(:,:), &
                                   uln(:,:),    ulj(:,:), &
                                rho_un(:,:), rho_uj(:,:), &
                                    Xn(:,:),     Xj(:,:), &
                                                 Gj(:,:), &
                                                v_j(:,:)
        REAL(8),INTENT(INOUT):: Psys(:,:), al_rhol11(:,:), Pprev
        !Calculation
        INTEGER:: n, s1, s2
        ! REAL(8)::

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n                  == Counter variable
        !                   s1                 == number of nodes
        !                   s2                 == number of junctions
        !**REAL**      ==>  Tn(),Tj()          == Junction and Node Temperatures [°F]
        !                   agn(),agj()        == Junction and Node Void Fractions
        !                   rhon(),rhoj()      == Junction and Node Densities [lbm/ft^3]
        !                   uln(),ulj()        == Junction and Node Internal Energies
        !                   rho_un(),rho_uj()  == Junction and Node Density * Internal Energy
        !                   Gj                 == Junction Mass Flux values
        !                   v_j                == junction velocities
        !                   Psys()             == System Pressure (t,k)
        !                   al_rhol11()        == alpha*rho for the 11th node (t,k)
        !                   Pprev              == Previous pressure for determining Pressure change in time
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        s1 = SIZE(Tn,1)
        s2 = SIZE(Tj,1)

        !Setting Nodes
        DO n = 1,s1
            Tn(n,1)     = Tn(n,2)
            agn(n,1)    = agn(n,2)
            rhon(n,1)   = rhon(n,2)
            uln(n,1)    = uln(n,2)
            rho_un(n,1) = rho_un(n,2)
            Xn(n,1)     = Xn(n,2)
        END DO

        !Setting Junctions
        DO n = 1,s2
            Tj(n,1)     = Tj(n,2)
            agj(n,1)    = agj(n,2)
            rhoj(n,1)   = rhoj(n,2)
            ulj(n,1)    = ulj(n,2)
            rho_uj(n,1) = rho_uj(n,2)
            Xj(n,1)     = Xj(n,2)
            Gj(n,1)     = Gj(n,2)
            v_j(n,1)    = v_j(n,2)
        END DO

        al_rhol11(1,1) = al_rhol11(1,2)
        Pprev          = Psys(1,1)
        Psys(1,1)      = Psys(1,2)
    END SUBROUTINE
END MODULE
