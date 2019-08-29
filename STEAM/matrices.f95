! Matrix-filling Module
MODULE matrices

    USE globals
    USE state
    USE conditional
    USE pressuredrop
    IMPLICIT NONE

CONTAINS

    REAL(8) FUNCTION a_value(rhoj,Aj,Vn,dt,j,n)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: j, n
        REAL(8),INTENT(IN) :: rhoj(:,:), Aj(:), Vn(:), dt

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  j       == Evaluation junction
        !                   n       == Evaluation node
        !**REAL**      ==>  rhoj()  == Junction densities
        !                   Aj()    == Junction cross-sectional areas
        !                   Vn()    == Nodal volumes
        !                   dt      == Current time step size
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        a_value = rhoj(j,1) * Aj(j) * dt / Vn(n)

    END FUNCTION


    REAL(8) FUNCTION b_value(rho_uj,Aj,Vn,dt,P,j,n)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: j, n
        REAL(8),INTENT(IN) :: rho_uj(:,:), Aj(:), Vn(:), dt, P

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  j       == Evaluation junction
        !                   n       == Evaluation node
        !**REAL**      ==>  rhoj()  == Junction densities
        !                   Aj()    == Junction cross-sectional areas
        !                   Vn()    == Nodal volumes
        !                   dt      == Current time step size
        !                   P       == Current system pressure
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        b_value = (rho_uj(j,1) + P*ff) * Aj(j) * dt / Vn(n)

    END FUNCTION


    SUBROUTINE Assign_a(a,aprev,rhoj,Aj,Vn,dt)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: rhoj(:,:), Aj(:), Vn(:), dt
        REAL(8),INTENT(OUT) :: a(:), aprev(:)

        INTEGER             :: i

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  i       == Counter variable
        !**REAL**      ==>  rhoj()  == Junction densities
        !                   Aj()    == Junction cross-sectional areas
        !                   Vn()    == Nodal volumes
        !                   dt      == Current time step size
        !                   a()     == a_n matrix
        !                   aprev() == a_n-1 matrix
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        DO i = 1, SIZE(a)
            a(i) = a_value(rhoj,Aj,Vn,dt,i,i)
            IF(i .EQ. 1) THEN !Node 1 (inlet plenum)
                aprev(i) = a_value(rhoj,Aj,Vn,dt,10,i)
            ELSE IF(i .EQ. 10) THEN !Node 10 (downcomer)
                aprev(i) = a_value(rhoj,Aj,Vn,dt,11,i)
            ELSE IF( i .EQ. 11) THEN !Node 11 (steam dome)
                aprev(i) = a_value(rhoj,Aj,Vn,dt,9,i)
            ELSE !Everywhere else. These follow the uniform pattern
                aprev(i) = a_value(rhoj,Aj,Vn,dt,i-1,i)
            END IF

            !! EX: aprev(2) is a(n,n-1) = a(2,1) in the notes
        END DO

    END SUBROUTINE


    SUBROUTINE Assign_b(b,bprev,rho_uj,Aj,Vn,dt,P)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: rho_uj(:,:), Aj(:), Vn(:), dt, P
        REAL(8),INTENT(OUT) :: b(:), bprev(:)

        INTEGER             :: i

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  i       == Counter variable
        !**REAL**      ==>  rhoj()  == Junction densities
        !                   Aj()    == Junction cross-sectional areas
        !                   Vn()    == Nodal volumes
        !                   dt      == Current time step size
        !                   P       == Current system pressure
        !                   b()     == b_n matrix
        !                   bprev() == b_n-1 matrix
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        DO i = 1, SIZE(b)
            b(i) = b_value(rho_uj,Aj,Vn,dt,P,i,i)
            IF(i .EQ. 1) THEN !Node 1 (inlet plenum)
                bprev(i) = b_value(rho_uj,Aj,Vn,dt,P,10,i)
            ELSE IF(i .EQ. 10) THEN !Node 10 (downcomer)
                bprev(i) = b_value(rho_uj,Aj,Vn,dt,P,11,i)
            ELSE IF( i .EQ. 11) THEN !Node 11 (steam dome)
                bprev(i) = b_value(rho_uj,Aj,Vn,dt,P,9,i)
            ELSE !Everywhere else. These follow the uniform pattern
                bprev(i) = b_value(rho_uj,Aj,Vn,dt,P,i-1,i)
            END IF
            !! EX: bprev(2) is b(n,n-1) = b(2,1) in the notes
        END DO

    END SUBROUTINE


    SUBROUTINE Assign_alpha(alpha,rhon,Tn,P,a,aprev,b,bprev)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: rhon(:,:), Tn(:,:), P, a(:), aprev(:), &
                               b(:), bprev(:)
        REAL(8),INTENT(OUT) :: alpha(:,:)

        INTEGER             :: n, sn
        REAL(8)             :: y, w

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n       == Counter variable
        !                   sn      == Number of nodes
        !**REAL**      ==>  alpha() == Alpha matrix
        !                   rhon()  == Nodal densities
        !                   Tn()    == Nodal temperatures
        !                   P       == Current system pressure
        !                   a()     == a_n matrix
        !                   aprev() == a_n-1 matrix
        !                   aspec() == a values that don't fit the general pattern
        !                   b()     == b_n matrix
        !                   bprev() == b_n-1 matrix
        !                   bspec() == b values that don't fit the general pattern
        !                   y,w     == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        sn = SIZE(rhon,1)

        DO n = 1,sn

            CALL drhodpsi(rhon(n,2),y,Tn(n,2),P)   !y is drhodpsi
            CALL drho_udpsi(rhon(n,2),w,Tn(n,2),P) !w is drho_udpsi

            !Main diagonal in 'A' of alpha(n,n)
            IF(n .LT. sn) alpha(n,n) = a(n)/y - b(n)/w

            !Second diagonal in 'A' of alpha(n,n-1)
            IF((n .GT. 1) .AND. (n .LT. sn-1)) &
                alpha(n,n-1) = -(aprev(n)/y - bprev(n)/w)

        END DO

        !alpha(1,10): node 1 junction 10
        CALL drhodpsi(rhon(1,2),y,Tn(1,2),P)   !y is drhodpsi
        CALL drho_udpsi(rhon(1,2),w,Tn(1,2),P) !w is drho_udpsi
        alpha(1,10) = -(aprev(1)/y - bprev(1)/w)

        !alpha(10,11): node 10 junction 11
        CALL drhodpsi(rhon(10,2),y,Tn(10,2),P)   !y is drhodpsi
        CALL drho_udpsi(rhon(10,2),w,Tn(10,2),P) !w is drho_udpsi
        alpha(10,11) = -(aprev(10)/y - bprev(10)/w)

    END SUBROUTINE


    SUBROUTINE Assign_beta(beta,rhon,Tn,P,agn)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: rhon(:,:), Tn(:,:), P, agn(:,:)
        REAL(8),INTENT(OUT) :: beta(:)

        INTEGER             :: n, sn
        REAL(8)             :: y, y2, w, w2

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n           == Counter variable
        !                   sn          == Number of nodes
        !**REAL**      ==>  beta()      == Beta matrix
        !                   rhon()      == Nodal densities
        !                   Tn()        == Nodal temperatures
        !                   P           == Current system pressure
        !                   agn()       == Nodal void fractions
        !                   y,y2,w,w2   == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        sn = SIZE(rhon,1)

        DO n = 1,sn
            ! IF(n .NE. sn) THEN !nodes 1-10 are the same equation
                CALL   drhodpsi(rhon(n,2), y,  Tn(n,2), P)
                CALL drho_udpsi(rhon(n,2), w,  Tn(n,2), P)
                CALL     drhodP(rhon(n,2), y2, Tn(n,2), P, agn(n,2)) ![lbm / ft^3 * psi]
                CALL   drho_udP(rhon(n,2), w2, Tn(n,2), P, agn(n,2)) ![lbm / ft^3 * psi]


                beta(n) = y2/y - w2/w

            ! END IF

        END DO

    END SUBROUTINE


    SUBROUTINE Assign_SMSE(SM,SE,rhon,rho_un,qn,dt,Vn,agj,rhoj,Gj,Xj,Tj,Psys,Aj,mfeed,hFD)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: rhon(:,:), rho_un(:,:), qn(:,:), dt, Vn(:), &
                                agj(:,:), rhoj(:,:), Gj(:,:), Xj(:,:), &
                                Tj(:,:), Psys(:,:), Aj(:), mfeed, hFD
        REAL(8),INTENT(OUT) :: SM(:), SE(:)

        INTEGER             :: n

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n           == Counter variable
        !**REAL**      ==>  rhon()      == Nodal densities
        !                   rho_un()    == Nodal densities * internal energy
        !                   qn()        == Nodal heat generation rate
        !                   dt          == Current time step
        !                   Vn()        == Nodal volumes
        !                   agj()       == Junction void fractions
        !                   rhoj()      == Junction densities
        !                   Gj()        == Junction mass fluxes
        !                   Xj()        == Junction qualities
        !                   Tj()        == Junction temperatures
        !                   Psys()      == Current system pressure
        !                   Aj()        == Junction cross-secitonal areas
        !                   mfeed       == Feed mass flow rate
        !                   hFD         == Feed enthalpy
        !                   SM()        == SM matrix
        !                   SE()        == SE matrix
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        DO n = 1,SIZE(rhon,1)
                SM(n) = rhon(n,1) ![lbm/ft^3]
            IF(n .EQ. 1) THEN
                SE(n) = rho_un(n,1) + (qn(n,1) - S_n(agj,rhoj,Gj,Xj,Tj,Psys(1,1),Aj,&
                                                     10,n)) * (dt / Vn(n)) ![Btu/ft^3]
            ELSE IF(n .EQ. 10) THEN
                SE(n) = rho_un(n,1) + (qn(n,1) - S_n(agj,rhoj,Gj,Xj,Tj,Psys(1,1),Aj,&
                                                     11,n)) * (dt / Vn(n)) ![Btu/ft^3]
            ELSE IF(n .EQ. 11) THEN
                SM(n) = rhon(n,1) + mfeed*dt/Vn(11) ![lbm/ft^3]
                SE(n) = rho_un(n,1) + (S_n(agj,rhoj,Gj,Xj,Tj,Psys(1,1),Aj,9,11) &
                                    + mfeed*hFD) * (dt / Vn(n)) ![Btu/ft^3]
                !!Due to the homogenous flow assumption in the separator, S_n is
                !! zero for this case

            ELSE
                SE(n) = rho_un(n,1) + (qn(n,1) - S_n(agj,rhoj,Gj,Xj,Tj,Psys(1,1),Aj,&
                                                     n-1,n)) * (dt / Vn(n)) ![Btu/ft^3]
            END IF

        END DO

    END SUBROUTINE


    SUBROUTINE Assign_M(M,Tn,Psys,agn,a9,a11,a12,SM)
        IMPLICIT NONE

        REAL(8),INTENT(IN)    :: Tn(:,:), Psys(:,:), agn(:,:), a9, a11, a12, SM(:)
        REAL(8),INTENT(OUT)   :: M(:)

        REAL(8)               :: y, w

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  Tn()    == Nodal temperatures
        !                   Psys()  == System pressure
        !                   agn()   == Nodal void fractions
        !                   a9()    == Special 'a' value, steam dome
        !                   a11()   == Special 'a' value, steam dome
        !                   a12()   == Special 'a' value, steam dome
        !                   SM()    == SM matrix (steam dome equations)
        !                   M()     == M matrix (steam dome equations)
        !                   y,w     == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        y    = SD_drhodag(Tn(11,2),Psys(1,2)) !Drho/Dalpha @ j11, k [lbm / ft^3]

        M(1) = SD_drhodul(Tn(11,2),Psys(1,2),agn(11,2)) / y ![ft^3 / Btu]

        M(2) = (SD_drhodP(Tn(11,2),Psys(1,2),agn(11,2))/144) / y  ![ft^2 / lbf]

        M(3) = a12 / y ![hr / ft]
        M(4) = a11 / y ![hr / ft]
        M(5) = a9  / y ![hr / ft]

        w    = SD_rho(Tn(11,2),Psys(1,2),agn(11,2)) ![lbm / ft^3]
        M(6) = ((SM(11) - w) / y) ![dimensionless]

    END SUBROUTINE


    SUBROUTINE Assign_E(E,Tn,Psys,agn,b9,b11,b12,SE)
        IMPLICIT NONE

        REAL(8),INTENT(IN)    :: Tn(:,:), Psys(:,:), agn(:,:), b12, b11, b9, SE(:)
        REAL(8),INTENT(OUT)   :: E(:)

        REAL(8)               :: y, w

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  Tn()    == Nodal temperatures
        !                   Psys()  == System pressure
        !                   agn()   == Nodal void fractions
        !                   b9()    == Special 'b' value, steam dome
        !                   b11()   == Special 'b' value, steam dome
        !                   b12()   == Special 'b' value, steam dome
        !                   SE()    == SE matrix (steam dome equations)
        !                   E()     == E matrix (steam dome equations)
        !                   y,w     == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        y    = SD_drho_udag(Tn(11,2),Psys(1,2)) !D(rho*u)/Dalphag [Btu / ft^3]
        E(1) = (SD_drho_udul(Tn(11,2),Psys(1,2),agn(11,2)) / y) ![ft^3 / Btu]

        E(2) = (SD_drho_udP(Tn(11,2),Psys(1,2),agn(11,2))/144) / y  ![ft^2 / lbf]

        E(3) = b12 / y ![hr / ft]
        E(4) = b11 / y ![hr / ft]
        E(5) = b9  / y ![hr / ft]

        w = SD_rho_u(Tn(11,2),Psys(1,2),agn(11,2))
        E(6) = ((SE(11) - w) / y) ![dimensionless]

    END SUBROUTINE


    SUBROUTINE Assign_LM(LM,Tn,Tj,Psys,agn,agj,Aj,al9,a11,S_Ml11)
        IMPLICIT NONE

        REAL(8),INTENT(IN)    :: Tn(:,:), Tj(:,:), Psys(:,:), agn(:,:), &
                                 agj(:,:), Aj(:), al9, a11, S_Ml11
        REAL(8),INTENT(OUT)   :: LM(:)

        REAL(8)               :: y, w

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  Tn()    == Nodal temperatures
        !                   Tj()    == Juntion temperatures
        !                   Psys()  == System pressure
        !                   agn()   == Nodal void fractions
        !                   agj()   == Junction void fractions
        !                   aj()    == Junction cross-sectional areas
        !                   al9()   == Special 'a' value, steam dome
        !                   a11()   == Special 'a' value, steam dome
        !                   S_Ml11  == Steam dome value
        !                   LM()    == LM matrix (steam dome equations)
        !                   y,w     == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        y     = SD_dalrholdag(Tn(11,2),Psys(1,2)) !D(alphal*rhol)/Dalphag [lbm / ft^3]
        LM(1) = SD_dalrholdul(Tn(11,2),Psys(1,2),agn(11,2)) / y ![ft^3 / Btu]

        LM(2) = (SD_dalrholdP(Tn(11,2),Psys(1,2),agn(11,2))/144) / y  ![ft^2 / lbf]

        LM(3) = 0.

        LM(4) = a11 / y ![hr / ft]

        LM(5) = al9 / y ![hr / ft]

        w = SD_al_rhol(Tn(11,2),Psys(1,2),agn(11,2)) ![lbm / ft^3]
        LM(6)  = ((S_Ml11 - w) / y) ![dimensionless]

    END SUBROUTINE


    SUBROUTINE Assign_Psi(psi,M,E,LM)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: M(:), E(:), LM(:)
        REAL(8),INTENT(OUT) :: psi(:)

        REAL(8)             :: y, w

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  M()     == M matrix   (steam dome equations)
        !                   E()     == E matrix   (steam dome equations)
        !                   LM()    == LM matrix  (steam dome equations)
        !                   psi()   == Psi matrix (steam dome equations)
        !                   y,w     == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        y        = (M(1) - LM(1)) ![ft^3 / Btu]
        w        = (E(1) - LM(1)) ![ft^3 / Btu]
        psi(1)   = ((M(3) - LM(3))/y - (E(3) - LM(3))/w) ![Btu*hr / ft^4]
        psi(2)   = ((M(4) - LM(4))/y - (E(4) - LM(4))/w) ![Btu*hr / ft^4]
        psi(3)   = ((M(5) - LM(5))/y - (E(5) - LM(5))/w) ![Btu*hr / ft^4]

    END SUBROUTINE


    SUBROUTINE Assign_Xi(Xi,M,E,LM,tcvpos,rhon,Tn,Psys,agn,SM,SE,uln)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: M(:), E(:), LM(:), tcvpos, rhon(:,:), Tn(:,:), &
                               Psys(:,:), agn(:,:), SM(:), SE(:), uln(:,:)
        REAL(8),INTENT(OUT) :: Xi(:)

        INTEGER             :: i, n
        REAL(8)             :: y, y2, w

        !-------Variable Definitions-------!
        !**INTEGER**   ==> i, n     == Counter variables
        !**REAL**      ==>  M()     == M matrix   (steam dome equations)
        !                   E()     == E matrix   (steam dome equations)
        !                   LM()    == LM matrix  (steam dome equations)
        !                   tcvpos  == Current position of the TCV
        !                   rhon()  == Nodal densities
        !                   Tn()    == Nodal temperatures
        !                   Psys()  == System pressure
        !                   agn()   == Nodal void fractions
        !                   SM()    == SM matrix
        !                   SE()    == SE matrix
        !                   uln()   == Nodal liquid specific internal energies
        !                   Xi()    == Xi matrix
        !                   y,y2,w  == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        Xi = 0.

        i = SIZE(Tn,1)
        DO n = 1,i
            IF(n .EQ. i) THEN ! .AND. (tcvPos .GT. 0)) THEN !node 11 being the SteamDome
                Xi(n) = (M(6) - LM(6))/(M(1)-LM(1)) - (E(6) - LM(6))/(E(1)-LM(1)) ![Btu / ft^3]
            ELSE
                CALL drhodpsi(rhon(n,2),y,Tn(n,2),Psys(1,2))
                CALL drho_udpsi(rhon(n,2),w,Tn(n,2),Psys(1,2))
                Xi(n) = (SM(n) - rho(rhon(n,2), uln(n,2), Psys(1,2), agn(n,2))) / y &
                        - (SE(n) - rho_u(rhon(n,2), uln(n,2), Psys(1,2), agn(n,2))) / w
            END IF
        END DO
        ! WRITE(*,*) "11 -> ", Xi(11)
    END SUBROUTINE


    SUBROUTINE Assign_SME(SME,Gj,Xj,Xn,agj,Psys,rhoj,rhon,Ln,v_j,Tj,Tn,De_n,Kn,dt)
        IMPLICIT NONE

        REAL(8),INTENT(IN)  :: Gj(:,:), Xj(:,:), Xn(:,:), agj(:,:), Psys(:,:), &
                               rhoj(:,:), rhon(:,:), Ln(:), v_j(:,:), Tj(:,:), &
                               Tn(:,:), dt, De_n(:), Kn(:)
        REAL(8),INTENT(OUT) :: SME(:,:)

        REAL(8)             :: vg1, vg2, vl1, vl2
        INTEGER             :: n

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  i, n     == Counter variables
        !**REAL**      ==>  M()     == M matrix   (steam dome equations)
        !                   E()     == E matrix   (steam dome equations)
        !                   LM()    == LM matrix  (steam dome equations)
        !                   tcvpos  == Current position of the TCV
        !                   rhon()  == Nodal densities
        !                   Tn()    == Nodal temperatures
        !                   Psys()  == System pressure
        !                   agn()   == Nodal void fractions
        !                   SM()    == SM matrix
        !                   SE()    == SE matrix
        !                   uln()   == Nodal liquid specific internal energies
        !                   Xi()    == Xi matrix
        !                   y,y2,w  == Placeholder variables
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        SME = 0
        DO n = 1,SIZE(Xn,1)
            CALL phase_vel(vg1,vl1,Gj(n-1,1),Xj(n-1,1),agj(n-1,1),Psys(1,1), n-1)
            CALL phase_vel(vg2,vl2,Gj(n,1),Xj(n,1),agj(n,1),Psys(1,1),n)


            IF(n .EQ. 10) THEN !Downcomer Negative Length
                CALL phase_vel(vg1,vl1,Gj(11,1),Xj(11,1),agj(11,1),Psys(1,1), 11)
                CALL phase_vel(vg2,vl2,Gj(n,1),Xj(n,1),agj(n,1),Psys(1,1),n)
                SME(1,1) = SME(1,1) + (Ln(n))*rhoj(n,1)*v_j(n,1)/dt - DP_acc(v_j(n,1),v_j(11,1),rhoj(n,1)) &
                           - DP_elev(rhoj(n,1), -Ln(n)) &
                           - DP_vel(agj(11,1),agj(n,1),rhoj(11,1) &
                           ,rhoj(n,1),vg1,vg2,vl1,vl2,Tj(11,1),Tj(n,1),Psys(1,1))
            ELSE IF(n .EQ. 1) THEN !Inlet Plenum has a Delta H of 0
                CALL phase_vel(vg1,vl1,Gj(10,1),Xj(10,1),agj(10,1),Psys(1,1), 10)
                CALL phase_vel(vg2,vl2,Gj(n,1),Xj(n,1),agj(n,1),Psys(1,1),n)
                SME(1,1) = SME(1,1) + Ln(n)*rhoj(n,1)*v_j(n,1)/dt - DP_acc(v_j(n,1),v_j(10,1),rhoj(n,1)) &
                           - DP_vel(agj(10,1),agj(n,1),rhoj(10,1) &
                           ,rhoj(n,1),vg1,vg2,vl1,vl2,Tj(10,1),Tj(n,1),Psys(1,1))
            ELSE IF(n .EQ. 9) THEN !Seperator has a Delta H of 0
                SME(1,1) = SME(1,1) + Ln(n)*rhoj(n,1)*v_j(n,1)/dt - DP_acc(v_j(n,1),v_j(n-1,1),rhoj(n,1)) &
                           - DP_vel(agj(n-1,1),agj(n,1),rhoj(n-1,1) &
                           ,rhoj(n,1),vg1,vg2,vl1,vl2,Tj(n-1,1),Tj(n,1),Psys(1,1))
            ELSE IF(n .EQ. 11) THEN !Steam Dome has a Delta H of 0
                CALL phase_vel(vg1,vl1,Gj(9,1),Xj(9,1),agj(9,1),Psys(1,1), 9)
                CALL phase_vel(vg2,vl2,Gj(n,1),Xj(n,1),agj(n,1),Psys(1,1),n)
                SME(1,1) = SME(1,1) + Ln(n)*rhoj(n,1)*v_j(n,1)/dt - DP_acc(v_j(n,1),v_j(9,1),rhoj(n,1)) &
                           - DP_vel(agj(9,1),agj(n,1),rhoj(9,1) &
                           ,rhoj(n,1),vg1,vg2,vl1,vl2,Tj(9,1),Tj(n,1),Psys(1,1))
            ELSE ! This includes the core and chimney
                SME(1,1) = SME(1,1) + Ln(n)*rhoj(n,1)*v_j(n,1)/dt - DP_acc(v_j(n,1),v_j(n-1,1),rhoj(n,1)) &
                           - DP_elev(rhoj(n,1),Ln(n)) &
                           - DP_vel(agj(n-1,1),agj(n,1),rhoj(n-1,1) &
                           ,rhoj(n,1),vg1,vg2,vl1,vl2,Tj(n-1,1),Tj(n,1),Psys(1,1))
            END IF

            SME(2,2) = SME(2,2) + DP_friction(v_j(n,2),rhon(n,:),rhoj(n,1),Xn(n,1),Tn(n,:), &
                       Psys(1,2),De_n(n),Ln(n)) * v_j(n,2) + &
                       DP_loss(v_j(n,2),rhon(n,1),rhoj(n,1),Xn(n,1),Tn(n,1), &
                       Psys(1,2),Kn(n)) * (v_j(n,2)/2)

        END DO


    END SUBROUTINE


END MODULE
