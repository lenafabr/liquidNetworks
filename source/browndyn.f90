MODULE BROWNDYN
  ! Utilities for running brownian dynamics simulations
  IMPLICIT NONE
  
CONTAINS


  SUBROUTINE RUNBROWNDYNSIM(NETP,NSTEP,DELT,DIFF,MOB,OUTFILE,PRINTEVERY,SNAPSHOTFILE,SNAPSHOTEVERY,APPENDSNAPSHOTS,DOBROWN)
    ! run a brownian dynamics simulation of network
    ! at each time step, call LangevinStepRK4 subroutine

    ! NSTEP: number of steps to simulate
    ! DELT: timestep
    ! OUTFILE: output file for printing out info (energy, end position)
    ! PRINTEVERY: how often to print info into OUTFILE
    ! SNAPSHOTFILE: file for dumping chain configuration snapshots
    ! SNAPSHOTEVERY: how often to dump snapshots
    ! APPENDSNAPSHOTS: append snapshots to existing files (otherwise, start from scratch and append as we go)
    ! DOBROWN: include brownian forces
    
    USE MINNETUTIL, ONLY : MINNET, GETMINNETFORCE, OUTPUTSNAPSHOT, OUTPUTMINNET, & 
                         & SETEDGELENS, OUTPUTTRACK, OUTPUTPIN, OUTPUTBOUNDS, &
                         & OUTPUTEBOUNDS
    USE KEYS, ONLY : MINNETOUTFILE, VERBOSE, TRACKNODES, TRACKFILE, PINFILE, &
                   & BNDFILE, EBNDFILE, BNDDEBUG, STARTSNAPSLATE
    
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: NETP
    INTEGER, INTENT(IN) :: NSTEP
    DOUBLE PRECISION, INTENT(IN) :: DELT, DIFF, MOB
    CHARACTER(LEN=*), INTENT(IN) :: OUTFILE, SNAPSHOTFILE
    INTEGER, INTENT(IN) :: PRINTEVERY, SNAPSHOTEVERY
    LOGICAL, INTENT(IN) :: APPENDSNAPSHOTS, DOBROWN
    
    INTEGER :: STEP
    DOUBLE PRECISION :: LTOT
    INTEGER :: SHORTEDGES(NETP%NEDGE+1)
    
    OPEN(FILE=OUTFILE,UNIT=88,STATUS='UNKNOWN')
    
    CALL SETEDGELENS(NETP,SHORTEDGES,LTOT)
    PRINT*, 'STEP, NETWORK LENGTH:', 0, LTOT
    WRITE(88,*) 0, LTOT
    
    ! Initial snapshot
    CALL OUTPUTSNAPSHOT(NETP,SNAPSHOTFILE,APPENDSNAPSHOTS)
    ! if tracking nodes, output track file
    IF (TRACKNODES) THEN
      CALL OUTPUTTRACK(NETP,TRACKFILE,APPENDSNAPSHOTS)
    ENDIF

    ! save pinned nodes, to make visualizing them easier for testing
    CALL OUTPUTPIN(NETP,PINFILE,APPENDSNAPSHOTS)

    IF (BNDDEBUG) THEN
      ! for debugging purposes
      CALL OUTPUTBOUNDS(NETP,BNDFILE,APPENDSNAPSHOTS)
      CALL OUTPUTEBOUNDS(NETP,EBNDFILE,APPENDSNAPSHOTS)
    ENDIF
    
    DO STEP = 1,NSTEP

      ! dynamic step
      CALL LANGEVINSTEPRK4(NETP,DELT,DIFF,LTOT,DOBROWN)

      ! printing
      IF (MOD(STEP,PRINTEVERY).EQ.0) THEN
        ! Output information on network status
        PRINT*, 'STEP, NETWORK LENGTH, MAX LENGTH:', STEP, LTOT, NETP%EDGELEN(NETP%EIDMAX)
        WRITE(88,*) STEP, LTOT
      ENDIF

      ! saving snapshots
      IF ((MOD(STEP,SNAPSHOTEVERY).EQ.0).AND.(STEP.GT.STARTSNAPSLATE)) THEN          
        CALL OUTPUTSNAPSHOT(NETP,SNAPSHOTFILE,.TRUE.) ! output snapshot
        IF (TRACKNODES) THEN
          CALL OUTPUTTRACK(NETP,TRACKFILE,.TRUE.) ! output track nodes
        ENDIF
        CALL OUTPUTPIN(NETP,PINFILE,.TRUE.)
        IF (BNDDEBUG) THEN
          CALL OUTPUTBOUNDS(NETP,BNDFILE,.TRUE.)
          CALL OUTPUTEBOUNDS(NETP,EBNDFILE,.TRUE.)
        ENDIF
      ENDIF

    ENDDO

    ! output final network
    CALL OUTPUTMINNET(NETP,MINNETOUTFILE,.FALSE.)

    CLOSE(88)


  END SUBROUTINE RUNBROWNDYNSIM
  


  SUBROUTINE LANGEVINSTEPRK4(NETP,DELT,DIFF,LTOT,DOBROWN)
    ! propagate net forward in time, using a fourth-order Runge-Kutta method
    ! DELT: timestep
    ! LTOT: returns the final energy of the chain after the step
    ! DOBROWN: toggle whether to include brownian forces
    USE MT19937, ONLY : RNORM, GRND
    USE MINNETUTIL, ONLY : MINNET, GETMINNETFORCE, SETEDGELENSSIMPLE, SETEDGELENS, &
          & REARRANGE, GROWANDPIN, FINDMAXEDGE, UPDATEPREV
    USE KEYS, ONLY : NFIX, FIXNODES, GVEL, GROWTH, RNODES
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: NETP
    DOUBLE PRECISION, INTENT(IN) :: DELT, DIFF
    DOUBLE PRECISION, INTENT(OUT) :: LTOT
    LOGICAL, INTENT(IN) :: DOBROWN
    DOUBLE PRECISION :: S2DT, RPOS
    DOUBLE PRECISION, DIMENSION(NETP%NODEHIGH,NETP%DIM) :: POS0, K1POS, K2POS, K3POS, &
          & K4POS, FORCES, FBROWN, SARRAY
    INTEGER :: NC, I, DEG, EC, EID
    INTEGER :: SHORTEDGES(NETP%NEDGE+1)
    INTEGER :: NHIGH


    ! update previous positions and edgenodes
    CALL UPDATEPREV(NETP)

    NHIGH=NETP%NODEHIGH

    ! ZEROS = 0D0
    S2DT = SQRT(2*DIFF/DELT)

    POS0 = NETP%NODEPOS(1:NHIGH,:)

    ! get the brownian forces
    IF (DOBROWN) THEN     
      DO NC = 1,NHIGH
        ! if node active
        IF (NETP%NODEACT(NC)) THEN
          ! find brownian force
          DO I = 1,NETP%DIM
            FBROWN(NC,I) = RNORM()*S2DT
          ENDDO

          ! if node is on the boundary, dot product with edgedir to only allow sliding
          ! along boundary edges
          DEG=NETP%NODEDEG(NC)
          DO EC = 1,DEG
            EID = NETP%NODEEDGE(NC,EC)

            IF (NETP%EDGEBND(EID)) THEN
              FBROWN(NC,:) = DOT_PRODUCT(FBROWN(NC,:),NETP%EDGEDIR(EID,:))*NETP%EDGEDIR(EID,:)
              ! adjusted force, can exit this loop
              EXIT

            ENDIF

          ENDDO
        ENDIF
      ENDDO
    ELSE
      FBROWN = 0D0
    ENDIF

    ! --------- 1ST RK STEP---------------
    CALL GETMINNETFORCE(NETP,FORCES)
    K1POS = (FBROWN + FORCES)
    DO NC = 1,NHIGH
      IF (NETP%NODEGROW(NC).EQ.1) THEN
        K1POS(NC,:) = NETP%NODEDIR(NC,:)*GVEL
      ELSEIF (NETP%NODEFIX(NC)) THEN
        K1POS(NC,:) = 0D0
      ! if pinning nodes
      ELSEIF (NETP%NODEPIN(NC)) THEN
        K1POS(NC,:) = 0D0
      ENDIF
    ENDDO
    ! K1POS = MERGE(ZEROS,K1POS,NETP%NODEFIX)

    NETP%NODEPOS(1:NHIGH,:) = POS0 + DELT/2*K1POS
    CALL SETEDGELENSSIMPLE(NETP)

    ! --------- 2ND RK STEP---------------
    CALL GETMINNETFORCE(NETP,FORCES)
    K2POS = (FBROWN + FORCES)
    DO NC = 1,NHIGH
      IF (NETP%NODEGROW(NC).EQ.1) THEN
        K2POS(NC,:) = NETP%NODEDIR(NC,:)*GVEL
      ELSEIF (NETP%NODEFIX(NC)) THEN
        K2POS(NC,:) = 0D0
      ! if pinning nodes
      ELSEIF (NETP%NODEPIN(NC)) THEN
        K2POS(NC,:) = 0D0
      ENDIF
    ENDDO
    
    NETP%NODEPOS(1:NHIGH,:) = POS0 + DELT/2*K2POS
    CALL SETEDGELENSSIMPLE(NETP)

    ! --------- 3RD RK STEP---------------
    CALL GETMINNETFORCE(NETP,FORCES)
    K3POS = (FBROWN + FORCES)
    DO NC = 1,NHIGH
      IF (NETP%NODEGROW(NC).EQ.1) THEN
        K3POS(NC,:) = NETP%NODEDIR(NC,:)*GVEL
      ELSEIF (NETP%NODEFIX(NC)) THEN
        K3POS(NC,:) = 0D0
      ! if pinning nodes
      ELSEIF (NETP%NODEPIN(NC)) THEN
        K3POS(NC,:) = 0D0
      ENDIF
    ENDDO
    
    NETP%NODEPOS(1:NHIGH,:) = POS0 + DELT*K3POS
    CALL SETEDGELENSSIMPLE(NETP)

    ! --------- 4TH RK STEP---------------
    CALL GETMINNETFORCE(NETP,FORCES)
    K4POS = (FBROWN + FORCES)
    DO NC = 1,NHIGH
      IF (NETP%NODEGROW(NC).EQ.1) THEN
        K4POS(NC,:) = NETP%NODEDIR(NC,:)*GVEL
      ELSEIF (NETP%NODEFIX(NC)) THEN
        K4POS(NC,:) = 0D0
      ! if pinning nodes
      ELSEIF (NETP%NODEPIN(NC)) THEN
        K4POS(NC,:) = 0D0
      ENDIF
    ENDDO
    
    SARRAY = DELT/6*(K1POS + 2*K2POS + 2*K3POS + K4POS)

!    DO NC = 1,NHIGH
!      IF (NETP%NODEBND(NC)) THEN
!        ! PRINT*, SARRAY(NC,:)
!        ! PRINT*, NETP%NODEDIR(NC,:)
!        SARRAY(NC,:) = DOT_PRODUCT(SARRAY(NC,:),NETP%NODEDIR(NC,:))*NETP%NODEDIR(NC,:)
!        ! PRINT*, SARRAY(NC,:)
!      ENDIF
!    ENDDO

    NETP%NODEPOS(1:NHIGH,:) = POS0 + SARRAY

    ! instead of updatind sarray prior to updating nodepos
    ! now just use the calculated sarray and then project nodepos onto circle with the same angular position
    DO NC = 1,NHIGH
      IF (NETP%NODEBND(NC)) THEN
        IF (NETP%NODEACT(NC).AND.(.NOT.NETP%NODEFIX(NC))) THEN
          ! calculate radial coordinate of node
          RPOS = SQRT(SUM(NETP%NODEPOS(NC,:)**2))
          IF (RPOS.NE.RNODES) THEN
            ! shift node to be along circular boundary
            NETP%NODEPOS(NC,:) = RNODES*NETP%NODEPOS(NC,:)/RPOS
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    ! PRINT*, 'max displacement', MAXVAL(SARRAY)
    ! sets new edgelengths, keeps track of shortedges
    ! and returns the total edge length in the network
    CALL SETEDGELENS(NETP,SHORTEDGES,LTOT)

    ! if any shortedges from setedgelens call, then check for necessary rearrangements
    IF (SHORTEDGES(1).GT.0) THEN
      CALL REARRANGE(NETP,SHORTEDGES)
    ENDIF

    ! check for growth/pin/unpin events
    IF (GROWTH.GT.0D0) THEN
      CALL GROWANDPIN(NETP)
    ENDIF

    ! if we rearranged (or possibility of rearranging) should make sure EIDMAX is up to date
    CALL FINDMAXEDGE(NETP)

  END SUBROUTINE LANGEVINSTEPRK4


END MODULE BROWNDYN
