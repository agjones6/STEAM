! Steam Dome Functions Module
MODULE steamdome

    USE state
    IMPLICIT NONE

CONTAINS

    REAL(8) FUNCTION SD_rho(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   ==Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_rho = (1 - ag)*rhol(u,P) + ag*rhog(P) ![lbm/ft^3]

    END FUNCTION


    REAL(8) FUNCTION SD_rho_u(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_rho_u = (1 - ag)*rhol(u,P)*u + ag*rhog(P)*ug(P) ![Btu/ft^3]

    END FUNCTION


    REAL(8) FUNCTION SD_al_rhol(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_al_rhol = (1 - ag)*rhol(u,P) ![lbm/ft^3]

    END FUNCTION


    REAL(8) FUNCTION SD_drhodag(T,P)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_drhodag = -(rhol(u,P) - rhog(P)) ![lbm/ft^3]

    END FUNCTION


    REAL(8) FUNCTION SD_drhodul(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_drhodul = (1 - ag)*drholdu(u,P) ![lbm^2/(Btu*ft^3)]

    END FUNCTION


    REAL(8) FUNCTION SD_drhodP(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_drhodP = (1 - ag)*drholdP(u,P) + ag*drhogdP(P) ![lbm/(ft^3*psia)]

    END FUNCTION


    REAL(8) FUNCTION SD_drho_udag(T,P)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_drho_udag = -(rhol(u,P)*u - rhog(P)*ug(P)) ![Btu/ft^3]

    END FUNCTION


    REAL(8) FUNCTION SD_drho_udul(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_drho_udul = (1 - ag)*(rhol(u,P) + u*drholdP(u,P)) ![]

    END FUNCTION


    REAL(8) FUNCTION SD_drho_udP(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_drho_udP = (1 - ag)*(u*drholdP(u,P)) + ag*(ug(P)*drhogdP(P) + &
                                                     rhog(P)*dugdP(P))

    END FUNCTION


    REAL(8) FUNCTION SD_dalrholdag(T,P)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_dalrholdag = -rhol(u,P) ![lbm/ft^3]

    END FUNCTION


    REAL(8) FUNCTION SD_dalrholdul(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_dalrholdul = (1 - ag)*drholdu(u,P) ![Btu/ft^3]

    END FUNCTION


    REAL(8) FUNCTION SD_dalrholdP(T,P,ag)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: T, P, ag
        REAL(8):: u

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  T   == Current temperature
        !                   P   == Current pressure
        !                   ag  == Void fraction
        !                   u   == Internal Energy
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        u = uliq(T)
        SD_dalrholdP = (1 - ag)*drholdP(u,P) ![lbm/(ft^3*psia)]

    END FUNCTION

END MODULE
