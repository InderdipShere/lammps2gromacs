!!  Calculates probability distribution of all bonds and all angles
!!  Takes the bond and angle list from lammps data file and then calculates bond distance and angles 
!! accordingly.
! INPUT(21) input_Validate_flex.dat
! OUTPUT(31) bond_validate_flex.dat
! OUTPUT(32) ang_validate_flex.dat
!! 27April22: Code improvement, integer*8 was used to define bond and ang (lots of memory)
!! 13June22:  Code updated for the dihedral and inproper angles as well
!! 20March24: Update the angles from cos(deg) to deg 

PROGRAM VALIDATE_FLEX
  IMPLICIT NONE
  REAL*8, ALLOCATABLE, DIMENSION(:)::RX, RY, RZ
  REAL*8::  X, Y, Z, RIJ(3), RIJ2, RJK(3), BOX(3), BOX_INV(3), RIJ_MAX,&
            DBOND, DANG, DDHI, DIMP, DBOND_INV, DANG_INV, DDHI_INV, DIMP_INV,&
            DBOND_2, DANG_2, BOX_VEC(3,3), XY, XZ, YZ, &
            RKL(3), RIJK(3), RJKL(3), XJ,YJ,ZJ, XK,YK,ZK, BOX_VEC_INV(3,3), &
            TEMP_A(9), TEMP_B(9), G(3), RAD2DEG
  INTEGER, ALLOCATABLE, DIMENSION(:):: COUNT_B, COUNT_A,COUNT_D, COUNT_I,&
            B_LIST_TYPE, A_LIST_TYPE, D_LIST_TYPE, I_LIST_TYPE ,&
            B_LIST_AI, A_LIST_AI, D_LIST_AI, I_LIST_AI ,&
            B_LIST_AJ, A_LIST_AJ, D_LIST_AJ, I_LIST_AJ , &
            A_LIST_AK, D_LIST_AK, I_LIST_AK,&
            D_LIST_AL, I_LIST_AL 
  INTEGER*8, ALLOCATABLE::BOND_DIS(:,:), ANG_DIS(:,:), DHI_DIS(:,:), IMP_DIS(:,:)
  INTEGER:: I, J, K, L, HEAD_TRAJ, TAIL_TRAJ, HEAD_LIST, NFRAME, NSKIP, NATOM,&
            NBOND, NANG, NDIHED, NIMPROP, NBIN_B, NBIN_A, NBIN_D, NBIN_I,&
            NBONDT, NANGT, NDIHEDT, NIMPROPT, TRAJ_COUNT, TOTAL_SKIP,&
            IY, AI, AJ, AK, AL, IG, ABOX(3), SKIP_FRAME,&
            debug

  CHARACTER*50:: TRAJFILE, LISTFILE , GROFORMAT, BOND_FORMAT, ANG_FORMAT
  CHARACTER*4::  RAN_C1, RAN_C2

  GROFORMAT  = '(i5,2a5,i5,3f8.3)'
  OPEN(21,FILE="input_Validate_flex.dat")
  READ(21,*) TRAJFILE
  READ(21,*) HEAD_TRAJ, TAIL_TRAJ
  READ(21,*) NFRAME, NSKIP, SKIP_FRAME
  READ(21,*) NATOM 
  READ(21,*) (BOX(I), I = 1, 3)
  READ(21,*) XY, XZ, YZ
  READ(21,*) NBONDT, NANGT, NDIHEDT, NIMPROPT
  READ(21,*) NBOND,  NANG,  NDIHED,  NIMPROP
  READ(21,*) RIJ_MAX
  READ(21,*) NBIN_B, NBIN_A, NBIN_D, NBIN_I
  READ(21,*) LISTFILE
  READ(21,*) HEAD_LIST
  
  DBOND = RIJ_MAX/DBLE(NBIN_B)
  DANG  = 180.0/DBLE(NBIN_A)
  DDHI  = 180.0/DBLE(NBIN_D)
  DIMP  = 180.0/DBLE(NBIN_I)
  RAD2DEG = 180.0/3.141592653589793238462
  
  DBOND_INV = 1.0D0 /DBOND
  DANG_INV  = 1.0D0/DANG
  DDHI_INV  = 1.0D0/DDHI 
  DIMP_INV  = 1.0D0/DIMP
  DBOND_2   = DBOND/2.0D0
  BOX_INV   = 1.0D0/BOX ! BOX VECTOR
  BOX_VEC   = 0.0
  DO I = 1, 3
    BOX_VEC(I,I)= BOX(I)
  END DO
  BOX_VEC(1,2) = XY ; BOX_VEC(1,3)= XZ
  BOX_VEC(2,3) = YZ
  PRINT*, "BOX_VEC" 
  PRINT*, BOX_VEC(1,:)
  PRINT*, BOX_VEC(2,:)
  PRINT*, BOX_VEC(3,:)
  TEMP_A =  RESHAPE(BOX_VEC,(/9/))
  PRINT*, "BOX VECTOR 1D", TEMP_A
  CALL INVERSE_3_3 (TEMP_A, TEMP_B)
  PRINT*, TEMP_B
  BOX_VEC_INV = RESHAPE(TEMP_B,(/3,3/))
  PRINT*, "INVERSE BOX_VECTOR"
  PRINT*, BOX_VEC_INV(1,:)
  PRINT*, BOX_VEC_INV(2,:)
  PRINT*, BOX_VEC_INV(3,:)
  WRITE(BOND_FORMAT,'(A11,I0,A10)') '"(F10.8,2X,', NBONDT,'F10.8,1X)"'
  WRITE(ANG_FORMAT,*) "'(F10.8, 2X,", NANGT, "(F10.8, 1X))'"
  DEBUG = 100
  !!!! NEED TO ADD HERE
  
  
  PRINT*, "TRAJECTORY FILE ", TRIM(TRAJFILE)
  PRINT*, "TOTAL NO OF HEAD FRAME SKIP", SKIP_FRAME
  PRINT*, "NO OF BONDS, ANGLE, DIHEDRAL, IMPROPER"
  PRINT*, NBOND, NANG,  NDIHED, NIMPROP
  PRINT*, "TYPE OF BONDS, ANGLE, DIHEDRAL, IMPROPER"
  PRINT*, NBONDT, NANGT, NDIHEDT, NIMPROPT
  PRINT*, "NO OF BIN: ",     NBIN_B, NBIN_A, NBIN_D, NBIN_I
  PRINT*, "BIN THICKNESS",   DBOND, DANG, DDHI, DIMP
  PRINT*, "MAX BOND STRECH", RIJ_MAX
  PRINT*, "LIST FILE ",      TRIM(LISTFILE)
  PRINT*, "BOND_FROMAT",     TRIM(BOND_FORMAT)
  PRINT*, "ANG_FROMAT",      TRIM(ANG_FORMAT)

  ALLOCATE( COUNT_B(NBONDT),  BOND_DIS(NBONDT,0:NBIN_B), B_LIST_TYPE(NBOND),  B_LIST_AI(NBOND),  B_LIST_AJ(NBOND), &
            COUNT_A(NANGT),   ANG_DIS(NANGT,0:NBIN_A),        A_LIST_TYPE(NANG),   A_LIST_AI(NANG),   A_LIST_AJ(NANG),   A_LIST_AK(NANG), &
            COUNT_D(NDIHEDT), DHI_DIS(NDIHEDT, 0:NBIN_D),     D_LIST_TYPE(NDIHED), D_LIST_AI(NDIHED), D_LIST_AJ(NDIHED), D_LIST_AK(NDIHED), D_LIST_AL(NDIHED),  &
            COUNT_I(NIMPROPT),IMP_DIS(NIMPROPT,0:NBIN_I),     I_LIST_TYPE(NIMPROP),I_LIST_AI(NIMPROP),I_LIST_AJ(NIMPROP),I_LIST_AK(NIMPROP),I_LIST_AL(NIMPROP), &
            RX(NATOM), RY(NATOM), RZ(NATOM) )
  OPEN (22,FILE=TRIM(LISTFILE), STATUS='OLD')
  DO I = 1, HEAD_LIST
    READ(22,*) ! SKIP
  END DO
  READ(22,*) ! BONDLIST
  READ(22,*) ! BLANK
  !write(debug,*) "## Bond reading type"
  DO I =1 , NBOND
    READ(22,*) J, B_LIST_TYPE(I), B_LIST_AI(I), B_LIST_AJ(I)
    !write(debug,*) i, B_LIST_AI(I), B_LIST_AJ(I),  B_LIST_TYPE(I)
  END DO 
  READ(22,*) !BLANK
  READ(22,*) !ANGLE LIST
  READ(22,*) !BLANK
  DO I = 1, NANG 
    READ(22,*) J, A_LIST_TYPE(I), A_LIST_AI(I), A_LIST_AJ(I), A_LIST_AK(I)
  END DO
  READ(22,*) !BLANK
  READ(22,*) !DIHEDRAL HEADING
  READ(22,*) !BLANK
  DO I = 1, NDIHED
    READ(22,*) J, D_LIST_TYPE(I), D_LIST_AI(I), D_LIST_AJ(I), D_LIST_AK(I), D_LIST_AL(I)
  END DO
  READ(22,*) !BLANK
  READ(22,*) !IMPROPER HEADING
  READ(22,*) !BLANK
  DO I = 1, NIMPROP
    READ(22,*) J, I_LIST_TYPE(I), I_LIST_AI(I), I_LIST_AJ(I), I_LIST_AK(I), I_LIST_AL(I)
  END DO

  CLOSE(22)
  PRINT*, "PARAMETER READ PROPERLY: BONDS, ANGLE, DIHEDRAL, IMPROPER"  
 
  COUNT_B = 0 ;  COUNT_A = 0 ;  COUNT_D = 0 ;  COUNT_I = 0
  BOND_DIS= 0 ;  ANG_DIS = 0 ;  DHI_DIS = 0 ;  IMP_DIS = 0  
  TOTAL_SKIP = (NSKIP-1)*(HEAD_TRAJ+TAIL_TRAJ+NATOM)
  SKIP_FRAME =  SKIP_FRAME*(HEAD_TRAJ+TAIL_TRAJ+NATOM)

  OPEN(23,FILE=TRIM(TRAJFILE),STATUS='OLD')
  IF (NSKIP == 0) THEN
    PRINT*, "ERROR: NSKIP IS 0, MIN 1 IS REQUIRED"
    STOP
  END IF
  NFRAME =  INT(NFRAME/NSKIP)
  
  OPEN(31,FILE="bond_raw.dat")
  OPEN(32,FILE="angle_raw.dat")
  OPEN(33,FILE="dihedral_raw.dat")
  OPEN(34,FILE="improper_raw.dat")
  DO TRAJ_COUNT = 1, SKIP_FRAME
    READ(23,*) ! SKIP
  END DO 
  DO TRAJ_COUNT = 1, NFRAME
    PRINT*, "TRAJ_COUNT ", TRAJ_COUNT, "/",  NFRAME 
    DO I = 1, HEAD_TRAJ
      READ(23,*) ! SKIP
    END DO
    !write(debug,*) "# traj" , traj_count
    DO I = 1, NATOM  
      READ(23,GROFORMAT) J,RAN_C1,RAN_C2,K, RX(I), RY(I), RZ(I)
      !write(debug,groformat)  J,RAN_C1,RAN_C2,K, RX(I), RY(I), RZ(I)
    END DO    
    DO I =1, TAIL_TRAJ
      READ(23,*) !! SKIP
    END DO

    !! Bond distribution
    !write(debug,*) "# Bond_calculations"
    WRITE(31,*) "#", TRAJ_COUNT !!! BONDS
    DO I = 1, NBOND
      AI = B_LIST_AI(I)
      AJ = B_LIST_AJ(I)
      IY = B_LIST_TYPE(I)      
      RIJ(1) =  RX(AI) - RX(AJ)
      RIJ(2) =  RY(AI) - RY(AJ)
      RIJ(3) =  RZ(AI) - RZ(AJ)
      
      ABOX(3) = ANINT(RIJ(3)*BOX_INV(3))
      RIJ     =  RIJ - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RIJ(2)*BOX_INV(2))
      RIJ     = RIJ - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RIJ(1)*BOX_INV(1))
      RIJ     = RIJ - ABOX(1)*BOX_VEC(:,1)      
      RIJ2   = DSQRT(DOT_PRODUCT(RIJ,RIJ))
      IG     = INT(RIJ2*DBOND_INV)
      COUNT_B(IY)      =  COUNT_B(IY) + 1
      BOND_DIS(IY, IG) =  BOND_DIS(IY,IG) + 1
      WRITE(31,*) IY, AI, AJ, RIJ2
    END DO
    WRITE(31,*)

    !! Angle distribution
    WRITE(32,*) "#", TRAJ_COUNT
    DO I = 1, NANG
      AI = A_LIST_AI(I)
      AJ = A_LIST_AJ(I)
      AK = A_LIST_AK(I)
      IY = A_LIST_TYPE(I)
      X = RX(AJ)
      Y = RY(AJ)
      Z = RZ(AJ)
      RIJ(1) =  RX(AI) - X
      RIJ(2) =  RY(AI) - Y
      RIJ(3) =  RZ(AI) - Z
      RJK(1) =  RX(AK) - X
      RJK(2) =  RY(AK) - Y
      RJK(3) =  RZ(AK) - Z
      ABOX(3) = ANINT(RIJ(3)*BOX_INV(3))
      RIJ     = RIJ - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RIJ(2)*BOX_INV(2))
      RIJ     = RIJ - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RIJ(1)*BOX_INV(1))
      RIJ     = RIJ - ABOX(1)*BOX_VEC(:,1)

      
      ABOX(3) = ANINT(RJK(3)*BOX_INV(3))
      RJK     = RJK - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RJK(2)*BOX_INV(2))
      RJK     = RJK - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RJK(1)*BOX_INV(1))
      RJK     = RJK - ABOX(1)*BOX_VEC(:,1)
      
      RIJ2   =  DOT_PRODUCT(RIJ,RJK)/DSQRT(DOT_PRODUCT(RIJ, RIJ)*DOT_PRODUCT(RJK,RJK)) != COS(THETA)
      RIJ2   =  ACOS(REAL(RIJ2,4))*RAD2DEG
      IG     =  INT(RIJ2*DANG_INV) 
      COUNT_A(IY)      =  COUNT_A(IY)    + 1
      ANG_DIS(IY, IG)  =  ANG_DIS(IY,IG) + 1
      WRITE(32,*) IY, AI, AJ, AK, RIJ2
    END DO
    WRITE(32,*)

    ! WRITE(32,*) "#", TRAJ_COUNT  !!! ANGLES
    ! DO I = 1, NANG
    !   AI = A_LIST_AI(I)
    !   AJ = A_LIST_AJ(I)
    !   AK = A_LIST_AK(I)
    !   IY = A_LIST_TYPE(I)
    !   X = RX(AJ)
    !   Y = RY(AJ)
    !   Z = RZ(AJ)
    !   RIJ(1) =  RX(AI) - X
    !   RIJ(2) =  RY(AI) - Y
    !   RIJ(3) =  RZ(AI) - Z
    !   RJK(1) =  RX(AK) - X
    !   RJK(2) =  RY(AK) - Y
    !   RJK(3) =  RZ(AK) - Z
    !   RIJ    =  RIJ -ANINT(RIJ*BOX_INV)*BOX 
    !   RJK    =  RJK -ANINT(RJK*BOX_INV)*BOX 
    !   RIJ2   =  DOT_PRODUCT(RIJ,RJK)/DSQRT(DOT_PRODUCT(RIJ, RIJ)*DOT_PRODUCT(RJK,RJK)) != COS(THETA)
    !   IG     =  INT(RIJ2*DANG_INV) 
    !   COUNT_A(IY)      =  COUNT_A(IY)    + 1
    !   ANG_DIS(IY, IG)  =  ANG_DIS(IY,IG) + 1
    !   WRITE(32,*) IY, AI, AJ, AK, RIJ2
    ! END DO
    ! WRITE(32,*)

    WRITE(33,*) "#", TRAJ_COUNT  !!! DIHEDRAL
    DO I = 1, NDIHED
      AI = D_LIST_AI(I)
      AJ = D_LIST_AJ(I)
      AK = D_LIST_AK(I)
      AL = D_LIST_AL(I)
      IY = D_LIST_TYPE(I)
      XJ = RX(AJ)
      YJ = RY(AJ)
      ZJ = RZ(AJ)
      RIJ(1) =  XJ- RX(AI)
      RIJ(2) =  YJ- RY(AI) 
      RIJ(3) =  ZJ- RZ(AI) 
      !RIJ    =  RIJ -ANINT(RIJ*BOX_INV)*BOX
      ABOX(3) = ANINT(RIJ(3)*BOX_INV(3))
      RIJ     =  RIJ - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RIJ(2)*BOX_INV(2))
      RIJ     = RIJ - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RIJ(1)*BOX_INV(1))
      RIJ     = RIJ - ABOX(1)*BOX_VEC(:,1) 
      
      XK = RX(AK) 
      YK = RY(AK) 
      ZK = RZ(AK) 
      RJK(1) =  XK - XJ
      RJK(2) =  YK - YJ
      RJK(3) =  ZK - ZJ
      !RJK    =  RJK -ANINT(RJK*BOX_INV)*BOX 
      ABOX(3) = ANINT(RJK(3)*BOX_INV(3))
      RJK     =  RJK - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RJK(2)*BOX_INV(2))
      RJK     = RJK - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RJK(1)*BOX_INV(1))
      RJK     = RJK - ABOX(1)*BOX_VEC(:,1)
      
      RIJK(1) = RIJ(2)*RJK(3)-RIJ(3)*RJK(2)
      RIJK(2) = RIJ(3)*RJK(1)-RIJ(1)*RJK(3)
      RIJK(3) = RIJ(1)*RJK(2)-RIJ(2)*RJK(1)
      !RIJK    = RIJK -ANINT(RIJK*BOX_INV)*BOX     
      ABOX(3) = ANINT(RIJK(3)*BOX_INV(3))
      RIJK     =  RIJK - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RIJK(2)*BOX_INV(2))
      RIJK     = RIJK - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RIJK(1)*BOX_INV(1))
      RIJK     = RIJK - ABOX(1)*BOX_VEC(:,1) 

      
      RKL(1) = RX(AL)-XK
      RKL(2) = RY(AL)-YK
      RKL(3) = RZ(AL)-ZK
      !RKL    =  RKL -ANINT(RKL*BOX_INV)*BOX
      ABOX(3) = ANINT(RKL(3)*BOX_INV(3))
      RKL     =  RKL - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RKL(2)*BOX_INV(2))
      RKL     = RKL - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RKL(1)*BOX_INV(1))
      RKL     = RKL - ABOX(1)*BOX_VEC(:,1) 
      
      RJKL(1) = RJK(2)*RKL(3)-RJK(3)*RKL(2)
      RJKL(2) = RJK(3)*RKL(1)-RJK(1)*RKL(3)
      RJKL(3) = RJK(1)*RKL(2)-RJK(2)*RKL(1)
      !RJKL    = RJKL -ANINT(RJKL*BOX_INV)*BOX
      ABOX(3) = ANINT(RJKL(3)*BOX_INV(3))
      RJKL     =  RJKL - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RJKL(2)*BOX_INV(2))
      RJKL     = RJKL - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RJKL(1)*BOX_INV(1))
      RJKL     = RJKL - ABOX(1)*BOX_VEC(:,1) 
      

      RIJ2   =  DOT_PRODUCT(RIJK,RJKL)/DSQRT(DOT_PRODUCT(RIJK, RIJK)*DOT_PRODUCT(RJKL,RJKL)) != COS(THETA)
      RIJ2   =  ACOS(REAL(RIJ2,4))*RAD2DEG
      IG     =  INT(RIJ2*DDHI_INV)       
      COUNT_D(IY)      =  COUNT_D(IY)    + 1
      DHI_DIS(IY, IG)  =  DHI_DIS(IY,IG) + 1
      WRITE(33,*) IY, AI, AJ, AK, AL, RIJ2 !, ig, count_D(iy), dhi_dis(iy,ig)
    END DO
    WRITE(33,*)

    WRITE(34,*) "#", TRAJ_COUNT  !!! improper
    DO I = 1, NIMPROP
      AI = I_LIST_AI(I)
      AJ = I_LIST_AJ(I)
      AK = I_LIST_AK(I)
      AL = I_LIST_AL(I)
      IY = I_LIST_TYPE(I)
      XJ = RX(AJ)
      YJ = RY(AJ)
      ZJ = RZ(AJ)
      RIJ(1) =  XJ- RX(AI)
      RIJ(2) =  YJ- RY(AI) 
      RIJ(3) =  ZJ- RZ(AI) 
      !RIJ    =  RIJ -ANINT(RIJ*BOX_INV)*BOX 
      ABOX(3) = ANINT(RIJ(3)*BOX_INV(3))
      RIJ     =  RIJ - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RIJ(2)*BOX_INV(2))
      RIJ     = RIJ - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RIJ(1)*BOX_INV(1))
      RIJ     = RIJ - ABOX(1)*BOX_VEC(:,1) 
      XK = RX(AK) 
      YK = RY(AK) 
      ZK = RZ(AK) 
      RJK(1) =  XK - XJ
      RJK(2) =  YK - YJ
      RJK(3) =  ZK - ZJ
      !RJK    =  RJK -ANINT(RJK*BOX_INV)*BOX 
      ABOX(3) = ANINT(RJK(3)*BOX_INV(3))
      RJK     =  RJK - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RJK(2)*BOX_INV(2))
      RJK     = RJK - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RJK(1)*BOX_INV(1))
      RJK     = RJK - ABOX(1)*BOX_VEC(:,1)
      RIJK(1) = RIJ(2)*RJK(3)-RIJ(3)*RJK(2)
      RIJK(2) = RIJ(3)*RJK(1)-RIJ(1)*RJK(3)
      RIJK(3) = RIJ(1)*RJK(2)-RIJ(2)*RJK(1)
      !RIJK    = RIJK -ANINT(RIJK*BOX_INV)*BOX    
      ABOX(3) = ANINT(RIJK(3)*BOX_INV(3))
      RIJK     =  RIJK - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RIJK(2)*BOX_INV(2))
      RIJK     = RIJK - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RIJK(1)*BOX_INV(1))
      RIJK     = RIJK - ABOX(1)*BOX_VEC(:,1) 
      RKL(1) = RX(AL)-XK
      RKL(2) = RY(AL)-YK
      RKL(3) = RZ(AL)-ZK
      !RKL    =  RKL -ANINT(RKL*BOX_INV)*BOX 
      ABOX(3) = ANINT(RKL(3)*BOX_INV(3))
      RKL     =  RKL - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RKL(2)*BOX_INV(2))
      RKL     = RKL - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RKL(1)*BOX_INV(1))
      RKL     = RKL - ABOX(1)*BOX_VEC(:,1) 
      RJKL(1) = RJK(2)*RKL(3)-RJK(3)*RKL(2)
      RJKL(2) = RJK(3)*RKL(1)-RJK(1)*RKL(3)
      RJKL(3) = RJK(1)*RKL(2)-RJK(2)*RKL(1)
      !RJKL    = RJKL -ANINT(RJKL*BOX_INV)*BOX 
      ABOX(3) = ANINT(RJKL(3)*BOX_INV(3))
      RJKL     =  RJKL - ABOX(3)*BOX_VEC(:,3)
      ABOX(2) = ANINT(RJKL(2)*BOX_INV(2))
      RJKL     = RJKL - ABOX(2)*BOX_VEC(:,2) 
      ABOX(1) = ANINT(RJKL(1)*BOX_INV(1))
      RJKL     = RJKL - ABOX(1)*BOX_VEC(:,1) 

      RIJ2   =  DOT_PRODUCT(RIJK,RJKL)/DSQRT(DOT_PRODUCT(RIJK, RIJK)*DOT_PRODUCT(RJKL,RJKL)) != COS(THETA)
      RIJ2   =  ACOS(REAL(RIJ2,4))*RAD2DEG
      IG     =  INT(RIJ2*DIMP_INV) 
      COUNT_I(IY)      =  COUNT_I(IY)    + 1
      IMP_DIS(IY, IG)  =  IMP_DIS(IY,IG) + 1
      WRITE(34,*) IY, AI, AJ, AK, AL, RIJ2
    END DO
    WRITE(34,*)


    DO I = 1, TOTAL_SKIP
      READ(23,*) ! SKIP
    END DO
  END DO
  CLOSE(31)
  CLOSE(32)
  CLOSE(33)
  CLOSE(34)
  PRINT*, "TRAJECTORY READ "
  


  OPEN(31,FILE="bond_validate_flex.dat")
  DO I = 0, NBIN_B
    WRITE(31,'(F10.5,2x,*(E15.10, 1X))') DBLE(I)*DBOND+DBOND_2, (DBLE(BOND_DIS(J,I))/DBLE(COUNT_B(J)), J = 1, NBONDT)
    !print*, DBLE(I)*DBOND+DBOND_2, (DBLE(BOND_DIS(J,I))/DBLE(COUNT_B(J)), J = 1, NBONDT)
  END DO
  CLOSE(31)
  PRINT*, "BONDS DONE"

  OPEN(32,FILE="ang_validate_flex.dat")
  DO I =0, NBIN_A-1
    WRITE(32,'(F10.5,2x,*(E15.10, 1X))') DBLE(I)*DANG, (DBLE(ANG_DIS(J,I))/DBLE(COUNT_A(J)), J = 1, NANGT)
  END DO
  CLOSE(32)
  PRINT*,"ANGLES DONE"

  OPEN(33,FILE="dih_validate_flex.dat")
  DO I = 0, NBIN_D-1
    WRITE(33,'(F10.5,2x,*(E15.10, 1X))') DBLE(I)*DDHI , (DBLE(DHI_DIS(J,I))/DBLE(COUNT_D(J)), J = 1, NDIHEDT)
  END DO
  CLOSE(33)
  PRINT*,"DIHEDRAL DONE"

  OPEN(34,FILE="imp_validate_flex.dat")
  DO I = 0, NBIN_I-1  
    WRITE(34,'(F10.5,2x,*(E15.10, 1X))') DBLE(I)*DIMP, (DBLE(IMP_DIS(J,I))/DBLE(COUNT_I(J)), J = 1, NIMPROPT)
  END DO
  CLOSE(34)
  PRINT*,"IMPROPER DONE"
END PROGRAM VALIDATE_FLEX
SUBROUTINE INVERSE_3_3(A, B)
  REAL*8, DIMENSION(9), INTENT(IN):: A
  REAL*8, DIMENSION(9), INTENT(OUT):: B
  REAL*8, DIMENSION(9):: C, D
  
  REAL*8:: X1, Y1, Z1 , X2Y3,  & 
           X2, Y2, Z2 , Y2Z3, &
           X3, Y3, Z3, DET_A
  X1 = A(1)
  Y1 = A(2)
  Z1 = A(3)
  X2 = A(4)
  Y2 = A(5)
  Z2 = A(6)
  X3 = A(7)
  Y3 = A(8)
  Z3 = A(9)
  Y2Z3 = Y2*Z3 - Y3*Z2
  X2Y3 = X2*Y3 - Y2*X3
  C(1) = A(5)*A(9)-A(8)*A(6)
  C(2) = A(7)*A(6)-A(4)*A(9)
  C(3) = A(4)*A(8)-A(7)*A(5)
  C(4) = A(8)*A(3)-A(9)*A(2)
  C(5) = A(1)*A(9)-A(7)*A(3)
  C(6) = A(7)*A(2)-A(8)*A(1)
  C(7) = A(6)*A(2)-A(5)*A(3)
  C(8) = A(4)*A(3)-A(1)*A(6)
  C(9) = A(5)*A(1)-A(4)*A(2)
  D(1) = C(1) ; D(5) = C(5) ; D(9)=C(9) ! TRANSPOSE
  D(2) = C(4) ; D(3) = C(7) ; D(6)=C(8)
  D(4) = C(2) ; D(7) = C(3) ; D(8)=C(6)
  DET_A = A(1)*C(1)+ A(2)*C(2)+ A(3)*C(3)
  B = 0.0 
  !B(1) = 1.0D0/X1
  !B(4) = (-X2*Y2Z3-Z2*X2Y3)/ (X1*Z2*Y2Z3)
  !B(5) = Y2Z3 /(Z2*Y2Z3)
  !B(6) = -Y2/Y2Z3
  !B(7) = X2Y3/(X1*Y2Z3)
  !B(8) = -Y3/Y2Z3
  !B(9) = Y2/Y2Z3
  ! USING COFACTOR
  B =  D/DET_A
  RETURN
END SUBROUTINE INVERSE_3_3

 


