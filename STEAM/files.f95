! File Management Module
MODULE files

    IMPLICIT NONE

CONTAINS

    SUBROUTINE openfile(u,filename)
        IMPLICIT NONE

        INTEGER,INTENT(IN):: u
        INTEGER:: stat
        CHARACTER(LEN=*),INTENT(IN):: filename

        !-------Variable Definitions-------!
        !**INTEGER**   ==> u        = File ID
        !                  length   = Length of file name string
        !                  stat     = Info (error) variable
        !**REAL**      ==>
        !**CHARACTER** ==> filename = Name of the file
        !**LOGICAL**   ==>
        !----------------------------------!

        OPEN(UNIT=u,FILE=filename,STATUS="REPLACE",ACTION="WRITE",IOSTAT=stat)
        IF(stat .NE. 0) THEN
            WRITE(*,*)
            WRITE(*,'(3A)') "!!THE FILE '", filename, &
        					"' DID NOT OPEN PROPERLY!!"
            WRITE(*,'(A,I4)') "The error code is: ", stat
            STOP
        END IF

    END SUBROUTINE

    SUBROUTINE writeVel(Bmat,n1,el_time,loopde,q0,Psys,n2,n3,MOI)
        IMPLICIT NONE

        INTEGER,INTENT(IN):: n1, loopde, n2
        INTEGER, OPTIONAL, INTENT(IN):: n3
        INTEGER:: i, n4
        REAL(8),INTENT(IN):: Bmat(n1,1), el_time, q0, Psys
        REAL(8),OPTIONAL,INTENT(IN):: MOI(n2,2)

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n1      == Number of junctions
        !                   loopde  == Current loop
        !                   i       == Counter variable
        !                   n2      == size of MOI
        !                   n3      == old time (t) 1, new time (k) 2
        !**REAL**      ==>  Bmat    == B matrix (velocities, dP)
        !                   el_time == Elapsed time
        !                   q0      == Heat input
        !                   Psys    == Current system pressure
        !                   MOI     == Matrix of interest
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!
        n4 = 2

        WRITE(13,'(A,F10.6,A,I6,A,F6.3,A,F11.6)') "Elapsed Time = ", el_time,&
                                      "  ||  k loop = "   , loopde &
                                    , "  ||  q0 = "       , q0 , &
                                      "  ||  Psys = ", Psys

        DO i = 1, n1
            WRITE(13,'(F15.6)',ADVANCE='NO') Bmat(i,1)
            IF(i .EQ. n1-1) WRITE(13,'(A)',ADVANCE='NO') " | "
        END DO
        WRITE(13,*)
        IF(PRESENT(MOI))THEN
            DO i = 1,n2
                IF(PRESENT(n3)) THEN
                    WRITE(13,'(F15.6)',ADVANCE='NO') MOI(i,n3)
                ELSE
                    WRITE(13,'(F10.6)',ADVANCE='NO') MOI(i,n4)
                END IF
            END DO
            WRITE(13,*)
        END IF

        WRITE(13,*)
        ! IF(loopde .EQ. 20) THEN
        !     WRITE(13,'(2A)') "-------------------------------------------------------------------------------------------", &
        !                      "-------------------------------------------------------------------------------------------"
        !     WRITE(13,*)
        ! END IF

    END SUBROUTINE

    SUBROUTINE writeVel2(my_mat,c,el_time,Psys,dt, my_file)
        IMPLICIT NONE

        INTEGER,INTENT(IN):: my_file, c
        INTEGER:: i, n1
        REAL(8),INTENT(IN):: el_time, Psys, dt
        REAL(8),INTENT(INOUT):: my_mat(:,:)

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n1      == Number of junctions
        !                   my_file == File Number to write to
        !**REAL**      ==>  my_mat     == array of interest
        !                   el_time == elapsed time
        !                   Psys    == System Pressure
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        n1 = SIZE(my_mat,1)
        WRITE(my_file,'(F25.10,A)',ADVANCE='NO') el_time, ","
        WRITE(my_file,'(F25.10,A)',ADVANCE='NO') Psys, ","
        WRITE(my_file,'(F25.10,A)',ADVANCE='NO') dt, ","

        DO i = 1, n1
            IF(i .LT. n1) THEN
                WRITE(my_file,'(F25.10,A)',ADVANCE='NO') my_mat(i,c), ","
            ELSE
                WRITE(my_file,'(F25.10,A)',ADVANCE='NO') my_mat(i,c)
            END IF
        END DO
        WRITE(my_file,*)
        ! WRITE(my_file,*)


    END SUBROUTINE

    SUBROUTINE writeSingle(my_file,val1,val2,val3,val4,val5,val6,val7,val8,val9,val10)
        IMPLICIT NONE

        ! Filename
        INTEGER,INTENT(IN):: my_file
        ! Up to 10 Inputs to be printed
        REAL(8),OPTIONAL,INTENT(IN):: val1,val2,val3,val4,val5,val6,val7,val8,val9,val10


        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n1      == Number of junctions
        !                   my_file == File Number to write to
        !**REAL**      ==>  my_mat     == array of interest
        !                   el_time == elapsed time
        !                   Psys    == System Pressure
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!


        IF(PRESENT(val1)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val1, ","
        IF(PRESENT(val2)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val2, ","
        IF(PRESENT(val3)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val3, ","
        IF(PRESENT(val4)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val4, ","
        IF(PRESENT(val5)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val5, ","
        IF(PRESENT(val6)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val6, ","
        IF(PRESENT(val7)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val7, ","
        IF(PRESENT(val8)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val8, ","
        IF(PRESENT(val9)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val9, ","
        IF(PRESENT(val10)) WRITE(my_file,'(F25.10,A)',ADVANCE='NO') val10, ","
        WRITE(my_file,*)

    END SUBROUTINE

    SUBROUTINE writeTDVars(Tj,Tn,agj,agn,rhoj,rhon,ulj,uln,rho_uj,rho_un)
        IMPLICIT NONE

        REAL(8),INTENT(IN):: Tj(:,:), Tn(:,:), &
                             agj(:,:), agn(:,:), rhoj(:,:), rhon(:,:), &
                             ulj(:,:), uln(:,:), rho_uj(:,:), rho_un(:,:)

        !-------Variable Definitions-------!
        !**INTEGER**   ==>
        !**REAL**      ==>  Tn()    ==  Temperatures:      node     [°F]
        !                   Tj()    ==  Temperatures:      junction [°F]
        !                   agn()   ==  Void Fractions:    node     []
        !                   agj()   ==  Void Fractions:    junction []
        !                   rhon()  ==  Densities:         node     [lbm/ft^3]
        !                   rhoj()  ==  Densities:         junction [lbm/ft^3]
        !                   uln()   ==  Internal Energies: node     [Btu/lbm]
        !                   ulj()   ==  Internal Energies: junction [Btu/lbm]
        !                   rho_un()==  rho*u:             node     [Btu/ft^3]
        !                   rho_uj()==  rho*u:             junction [Btu/ft^3]
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        CALL openfile(12,'TDVars.txt')

        WRITE(12,'(A)') "________________________TEMPERATURE________________________"
        CALL loopwrite(Tj,Tn)
        WRITE(12,'(A)') "_______________________VOID FRACTION_______________________"
        CALL loopwrite(agj,agn)
        WRITE(12,'(A)') "__________________________DENSITY__________________________"
        CALL loopwrite(rhoj,rhon)
        WRITE(12,'(A)') "___________________LIQUID INTERNAL ENERGY__________________"
        CALL loopwrite(ulj,uln)
        WRITE(12,'(A)') "___________________________RHO*U___________________________"
        CALL loopwrite(rho_uj,rho_un)

        CLOSE(12)

    END SUBROUTINE

    SUBROUTINE SnapShot(VEL,DENS_n,DENS_j,DENSu_n,DENSu_j,VOID_n,VOID_j, &
                        intENERGY_n,intENERGY_j,QUAL_n,QUAL_j,TEM_n,TEM_j,PRESS, &
                        VOID_DENS11,HEATq, TIME, DELTA_t, VALVE, FEEDFLOW, &
                        MY_COLUMN, FILE_NUM)
        !This subroutine will pull all of the variables needed to start the reactor at the given
        ! point and put them in a file of choice. Each variable is a COLUMN of values
        IMPLICIT NONE

        !Inputs
        INTEGER,INTENT(IN):: FILE_NUM, MY_COLUMN

        REAL(8),INTENT(IN):: VEL(:,:), DENS_n(:,:), DENS_j(:,:), DENSu_n(:,:), &
                             DENSu_j(:,:), VOID_n(:,:), VOID_j(:,:), intENERGY_n(:,:), &
                             intENERGY_j(:,:), QUAL_n(:,:), QUAL_j(:,:), &
                             TEM_n(:,:), TEM_j(:,:), PRESS(:,:), VOID_DENS11(:,:), &
                             HEATq, TIME, DELTA_t, VALVE, FEEDFLOW
        !Calculations
        INTEGER:: i, n1, c
        REAL(8):: y

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n1          == Number of junctions
        !                   my_file     == File Number to write to
        !                   MY_COLUMN   == Column desired to write
        !                   c           == Dummy variable to make the MY_COLUMN name smaller
        !**REAL**      ==>  VEL         == Junction Velocity array
        !                   DENS        == Nodal/Junction Density array
        !                   DENSu       == Nodal/Junction Density*internal energy array
        !                   VOID        == Nodal/Junction Void array
        !                   intEnergy   == Nodal/Junction Internal Energy Array
        !                   QUAL        == Nodal/Junction Quality Array
        !                   TEM         == Nodal/Junction Temperature Array
        !                   PRESS       == System Pressure Array
        !                   VOID_DENS11 == al_rho11 array
        !                   HEATq       == q0 value
        !                   TIME        == Current Elapsed time of the system
        !                   DELTA_t     == Current dt of the system
        !                   y           == dummy variable
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        !----------------------------------- OUTPUT ------------------------------------
        ! v_j   |rhon  rhoj|  rho_u    ag      ul       X       T     Psysk  al_rhol11     q0     el_time    dt     tcvPos
        !  .    |  .     . |    .       .       .       .       .     Psyst      .          .         .       .       .
        !  .    |  .     . |    .       .       .       .       .       .        .          .         .       .       .
        !  .    |  .     . |    .       .       .       .       .       .        .          .         .       .       .
        !  .    | n11    . |   n/j     n/j     n/j     n/j     n/j      .        .          .         .       .       .
        ! j12   |  0    j12|    j       j       j       j       j       0        0          0         0       0       0
        !-------------------------------------------------------------------------------
        ! Junction and Node Properties are both captured shown with |rhon  rhoj| column.
        !   All of the other junction/nodes follow the same format


        n1 = SIZE(VEL,1)
        y  = 0 !Dummy variable to place 0s in the 12th colum
        c  = MY_COLUMN !dummy variable just to make the variable name smaller

        DO i = 1, n1
            IF(i .LT. n1) THEN !(1-11)
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VEL(i,c)         !Column 1 *ONLY JUNCTION*
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') DENS_n(i,c)      !Column 2
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') DENS_j(i,c)      !Column 3
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') DENSu_n(i,c)     !Column 4
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') DENSu_j(i,c)     !Column 5
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VOID_n(i,c)      !Column 6
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VOID_j(i,c)      !Column 7
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') intENERGY_n(i,c) !Column 8
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') intENERGY_j(i,c) !Column 9
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') QUAL_n(i,c)      !Column 10
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') QUAL_j(i,c)      !Column 11
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') TEM_n(i,c)       !Column 12
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') TEM_j(i,c)       !Column 13
            ELSE ! (12)
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VEL(i,c)         !Column 1 *ONLY JUNCTION*
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 2
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') DENS_j(i,c)      !Column 3
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 4
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') DENSu_j(i,c)     !Column 5
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 6
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VOID_j(i,c)      !Column 7
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 8
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') intENERGY_j(i,c) !Column 9
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 10
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') QUAL_j(i,c)      !Column 11
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 12
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') TEM_j(i,c)       !Column 13
            END IF
            IF(i .LE. 1) THEN !(1 Value)
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') PRESS(1,1)       !Column 14
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VOID_DENS11(1,1) !Column 15
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') HEATq            !Column 16
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') TIME             !Column 17
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') DELTA_t          !Column 18
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VALVE            !Column 19
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') FEEDFLOW         !Column 20
            ELSEIF(i .EQ. 2) THEN
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') Press(1,2)       !Column 14
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') VOID_DENS11(1,2) !Column 15
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 16
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 17
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 18
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 19
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 20
            ELSE
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 14
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 15
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 16
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 17
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 18
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 19
                WRITE(FILE_NUM,'(F50.35)',ADVANCE='NO') y                !Column 20
            END IF
            WRITE(FILE_NUM,*)
        END DO

    END SUBROUTINE

    SUBROUTINE ReadSnap(VEL,DENS_n,DENS_j,DENSu_n,DENSu_j,VOID_n,VOID_j, &
                        intENERGY_n,intENERGY_j,QUAL_n,QUAL_j,TEM_n,TEM_j,PRESS, &
                        VOID_DENS11,HEATq, TIME, DELTA_t, VALVE, FEEDFLOW, &
                        FILE_NUM)
        ! This subroutine reads in snapshots that are taken

        IMPLICIT NONE

        !Inputs
        INTEGER,INTENT(IN):: FILE_NUM

        REAL(8),INTENT(INOUT):: VEL(:,:), DENS_n(:,:),DENS_j(:,:), DENSu_n(:,:), &
                                DENSu_j(:,:), VOID_n(:,:),VOID_j(:,:), &
                                intENERGY_n(:,:),intENERGY_j(:,:), QUAL_n(:,:),QUAL_j(:,:), &
                                TEM_n(:,:),TEM_j(:,:), PRESS(:,:), VOID_DENS11(:,:), HEATq, &
                                TIME, DELTA_t, VALVE, FEEDFLOW

        !Calculations
        INTEGER:: i, n1, c
        REAL(8):: y, dum

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  n1          == Number of junctions
        !                   my_file     == File Number to write to
        !                   c           == Dummy variable to make the MY_COLUMN name smaller
        !**REAL**      ==>  VEL         == Junction Velocity array
        !                   DENS        == Nodal/Junction Density array
        !                   DENSu       == Nodal/Junction Density*internal energy array
        !                   VOID        == Nodal/Junction Void array
        !                   intEnergy   == Nodal/Junction Internal Energy Array
        !                   QUAL        == Nodal/Junction Quality Array
        !                   TEM         == Nodal/Junction Temperature Array
        !                   PRESS       == System Pressure Array
        !                   VOID_DENS11 == al_rho11 array
        !                   HEATq       == q0 value
        !                   TIME        == Current Elapsed time of the system
        !                   DELTA_t     == Current dt of the system
        !                   VALVE       == Current tcv valve position
        !                   y           == dummy variable
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        n1 = SIZE(VEL,1)
        y  = 0 !Dummy variable to place 0s in the 12th colum
        c  = 1 !dummy variable just to make the variable name smaller


        DO i = 1, n1
            IF(i .LE. 1) THEN !(1)
                READ(FILE_NUM,*) VEL(i,c), DENS_n(i,c),DENS_j(i,c), DENSu_n(i,c), &
                                 DENSu_j(i,c), VOID_n(i,c),VOID_j(i,c), &
                                 intENERGY_n(i,c),intENERGY_j(i,c), &
                                 QUAL_n(i,c),QUAL_j(i,c), TEM_n(i,c),TEM_j(i,c), &
                                 PRESS(1,1), VOID_DENS11(1,1), HEATq, &
                                 TIME, DELTA_t, VALVE, FEEDFLOW
            ELSEiF(i .EQ. 2) THEN !(3)
                READ(FILE_NUM,*) VEL(i,c), DENS_n(i,c),DENS_j(i,c), DENSu_n(i,c), &
                                 DENSu_j(i,c), VOID_n(i,c),VOID_j(i,c), &
                                 intENERGY_n(i,c),intENERGY_j(i,c), &
                                 QUAL_n(i,c),QUAL_j(i,c), TEM_n(i,c),TEM_j(i,c), &
                                 PRESS(1,2), VOID_DENS11(1,2)
            ELSEIF(i .LT. n1) THEN  !(3-11)
                READ(FILE_NUM,*) VEL(i,c), DENS_n(i,c),DENS_j(i,c), DENSu_n(i,c), &
                                 DENSu_j(i,c), VOID_n(i,c),VOID_j(i,c), &
                                 intENERGY_n(i,c),intENERGY_j(i,c), &
                                 QUAL_n(i,c),QUAL_j(i,c), TEM_n(i,c),TEM_j(i,c)
            ELSE !(12)
                READ(FILE_NUM,*) VEL(i,c),dum,DENS_j(i,c),dum, DENSu_j(i,c),dum,VOID_j(i,c), &
                                 dum,intENERGY_j(i,c), dum, QUAL_j(i,c), dum, TEM_j(i,c)
            END IF

        END DO


    END SUBROUTINE

    SUBROUTINE loopwrite(Aj,An)
        IMPLICIT NONE

        INTEGER:: i, sj
        REAL(8),INTENT(IN):: Aj(:,:), An(:,:)

        !-------Variable Definitions-------!
        !**INTEGER**   ==>  sj      == Number of junctions
        !                   i       == Counter variable
        !**REAL**      ==>  Aj()    == Junction array
        !                   An()    == Node array
        !**CHARACTER** ==>
        !**LOGICAL**   ==>
        !----------------------------------!

        sj = SIZE(Aj,1)

        WRITE(12,*)
        WRITE(12,'(A)') "|------------ t ------------| |------------ k ------------|"
        WRITE(12,'(A)') "|  n  | Junction |   Node   | |  n  | Junction |   Node   |"
        WRITE(12,'(A)') "|---------------------------| |---------------------------|"

        DO i = 1,sj
            IF(i .EQ. sj) THEN
                WRITE(12,'(A,2(I2,A,D8.3,A))') &
                    "| ",i,"  | ",Aj(i,1)," |          | | ",i,&
                    "  | ",Aj(i,2)," |          |"
            ELSE
                WRITE(12,'(A,2(I2,A,2(D8.3,A)))') &
                    "| ",i,"  | ",Aj(i,1)," | ",An(i,1)," | | ",i,"  | ",&
                    Aj(i,2)," | ",An(i,2)," |"
            END IF
        END DO

        WRITE(12,'(A)') "|---------------------------------------------------------|"
        WRITE(12,*)
        WRITE(12,*)

    END SUBROUTINE

    SUBROUTINE ReadProf(MATRIX, TIME, FILE, TIME_MULT  )
    ! This subroutine will read in a profile in the form of a tsv file
    IMPLICIT NONE

    ! Inputs
    REAL(8),ALLOCATABLE,INTENT(INOUT):: MATRIX(:), TIME(:)
    REAL(4), OPTIONAL, INTENT(IN):: TIME_MULT
    INTEGER, INTENT(IN):: FILE

    ! Local
    INTEGER:: s, i, io
    REAL(8):: dum1, dum2, tm

    !-------Variable Definitions-------!
    !**INTEGER**   ==>
    !**REAL**      ==>  MATRIX(:)  == The matrix values wish to be read in
    !                   TIME(:)    == The time at which the matrix values occur
    !                   FILE       == The file number desired to be read
    !                  [TIME_MULT] == This is an optional multiplier on all of the time values
    !                   tm         == Permanent Time Multiplier that is set to 1 without TIME_MULT input
    !                   s          == The number of points in the input file
    !**CHARACTER** ==>
    !**LOGICAL**   ==>
    !----------------------------------!

    ! Determining the time multiplier was actually inputted
    IF(PRESENT(TIME_MULT)) THEN
        tm = TIME_MULT
    ELSE
        tm = 1
    END IF

    ! Allocating sizes for the arrays. This should be made dynamic
    s = 0
    DO
        READ(FILE,*,iostat=io)
        IF (io/=0) EXIT
        s = s + 1
    END DO
    REWIND(UNIT = FILE)
    ALLOCATE(TIME(s))
    ALLOCATE(MATRIX(SIZE(TIME)))

    ! Reading in values from the FILE given
    DO i = 1,s
        READ(FILE,*) dum1, dum2
        TIME(i)   = dum1 * tm
        MATRIX(i) = dum2
    END DO

    END SUBROUTINE

    

!     !NOT CURRENTLY IN USE
!     SUBROUTINE writeoutput(zj,n1,zn,n2,Aj,Tj,hj,agj,n3,An,Vn,Tn,hn,agn,&
!                            a,aprev,b,bprev,alpha,n4,aspec,bspec,n5)
!         IMPLICIT NONE
!
!         INTEGER,INTENT(IN):: n1, n2, n3, n4, n5
!         INTEGER:: i, j
!         REAL(8),INTENT(IN):: zj(n1), zn(n2), Aj(n3), Tj(n3), hj(n3), agj(n3), &
!                              An(n4), Vn(n4), Tn(n4), hn(n4), agn(n4), a(n4), &
!                              aprev(n4), b(n4), bprev(n4), alpha(n4,n4), &
!                              aspec(n5), bspec(n5)
!
!         !-------Variable Definitions-------!
!         !**INTEGER**   ==>  n1      == Number of core junctions
!         !                   n2      == Number of core nodes
!         !                   n3      == Number of junctions
!         !                   n4      == Number of nodes
!         !                   n5      == Number of 'special' alphas (2)
!         !                   i,j     == Counter variable
!         !**REAL**      ==>  zj      == Core junction axial positions
!         !                   zn      == Core node axial positions
!         !                   Aj      == Junction cross-sectional areas
!         !                   Tj      == Junction bulk fluid temperatures
!         !                   hj      == Junction enthalpies
!         !                   agj     == Junction void fractions
!         !                   An      == Node cross-sectional areas
!         !                   Vn      == Node volumes
!         !                   Tn      == Node bulk fluid temperatures
!         !                   hn      == Node enthalpies
!         !                   agn     == Node void fractions
!         !                   a       == 'a_n' values, used in the alpha matrix
!         !                   aprev   == 'a_n-1' values, used in the alpha matrix
!         !                   b       == 'b_n' values, used in the alpha matrix
!         !                   bprev   == 'b_n-1' values,used in the alpha matrix
!         !                   alpha   == 'alpha' matrix, used in the 'A' matrix
!         !                   aspec   == Special a values: for alpha(1,10),(10,11)
!         !                   bspec   == Special b values: for alpha(1,10),(10,11)
!         !**CHARACTER** ==>
!         !**LOGICAL**   ==>
!         !----------------------------------!
!
!         WRITE(12,'(A)') "|---------------------------------------------------|"
!         WRITE(12,'(A)') "|---------------CORE: AXIAL POSITIONS---------------|"
!         WRITE(12,'(A)') "|---------------------------------------------------|"
!         WRITE(12,'(A)') "|  n  |       JUNCTION       |         NODE         |"
!         WRITE(12,'(A)') "|---------------------------------------------------|"
!         DO i = 1,SIZE(zj)
!             IF(i .EQ. 1) THEN
!                 WRITE(12,'(A,I1,A,F6.3,A)') &
!                 "|  ",i, "  |        ",zj(i),"        |                      |"
!             ELSE
!                 WRITE(12,'(A,I1,A,2(F6.3,A))') &
!                 "|  ",i, "  |        ",zj(i),"        |        ",zn(i-1),&
!                 "        |"
!             END IF
!         END DO
!         WRITE(12,'(A)') "|---------------------------------------------------|"
!
!         WRITE(12,*)
!
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|"
!         WRITE(12,'(A)') "|----------------------------------GEOMETRY-----------------------------------|"
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|" !51
!         WRITE(12,'(A)') "|  n  | JUNCTION AREA [ft^2] |   NODE AREA [ft^2]   |   NODE VOLUME [ft^2]    |"
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|"
!         DO i = 1,SIZE(Aj)
!             IF(i .EQ. SIZE(Aj)) THEN
!                 WRITE(12,'(A,I2,A,D9.4,A)') &
!                 "| ",i, "  |      ",Aj(i),"       |                      |                         |"
!             ELSE
!                 WRITE(12,'(A,I2,A,3(D9.4,A))') &
!                 "| ",i, "  |      ",Aj(i),"       |      ",An(i),"       |        ", &
!                 Vn(i),"        |"
!             END IF
!         END DO
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|"
!
!         WRITE(12,*)
!
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|"
!         WRITE(12,'(A)') "|---------------------------THERMODYNAMIC VARIABLES---------------------------|"
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|"
!         WRITE(12,'(A)') "|     |    TEMPERATURE [°F]   |   ENTHALPY [Btu/lbm]  |     VOID FRACTION     |" !17each
!         WRITE(12,'(A)') "|  n  | JUNCTION  |    NODE   | JUNCTION  |    NODE   | JUNCTION  |    NODE   |"
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|"
!         DO i = 1,SIZE(Aj)
!             IF(i .EQ. SIZE(Aj)) THEN
!                 WRITE(12,'(A,I2,A,3(D9.4,A))') &
!                 "| ",i, "  | ",Tj(i)," |           | ",hj(i)," |           | ",agj(i)," |           |"
!             ELSE
!                 WRITE(12,'(A,I2,A,6(D9.4,A))') &
!                 "| ",i, "  | ",Tj(i)," | ",Tn(i)," | ",hj(i)," | ",hn(i)," | ",agj(i)," | ",agn(i)," |"
!             END IF
!         END DO
!         WRITE(12,'(A)') "|-----------------------------------------------------------------------------|"
!
!         WRITE(12,*)
!
!         DO i = 1,SIZE(alpha,1) !Row
!             DO j = 1,SIZE(alpha,1) !Column
!                 WRITE(12,'(F11.3,A)',ADVANCE='NO') alpha(i,j), ", "
!             END DO
!             WRITE(12,*)
!         END DO
!     END SUBROUTINE
!
END MODULE
