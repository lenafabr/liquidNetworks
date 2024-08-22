PROGRAM MAIN
  ! subroutines for testing parts of the code
  USE KEYS, ONLY : OUTFILE, SNAPSHOTFILE, MINNETFILE, MINNETOUTFILE, &
      & APPENDSNAPSHOTS, SNAPSHOTEVERY, & 
      & DELT, DIFF, MOB, BDSTEPS, BDPRINTEVERY, DOBROWN, CENTERENCLOSE, &
      & ACTION
  USE MINNETUTIL, ONLY : MINNET, MINNETFROMFILE, SETUPMINNET, GETMINNETFORCE, &
       & CENTENC, OUTPUTMINNET, OUTPUTSNAPSHOT, CLEANUPMINNET
  USE BROWNDYN, ONLY : RUNBROWNDYNSIM
  IMPLICIT NONE
  TYPE(MINNET), TARGET :: MINNETOBJ
  TYPE(MINNET), POINTER :: MINNETP

  MINNETP=>MINNETOBJ

  CALL READKEY

  CALL MINNETFROMFILE(MINNETP,MINNETFILE)

  IF (CENTERENCLOSE) THEN
    CALL CENTENC(MINNETP)
  ENDIF

  IF (ACTION.EQ.'BROWNDYN') THEN
    ! network dynamics, no particles on the network
    CALL RUNBROWNDYNSIM(MINNETP,BDSTEPS,DELT,DIFF,MOB,OUTFILE,BDPRINTEVERY,SNAPSHOTFILE,&
                        & SNAPSHOTEVERY,APPENDSNAPSHOTS,DOBROWN)
  ELSE
    PRINT*, 'ERROR IN MAIN: unidentified action', ACTION
    STOP 1
  ENDIF

  CALL CLEANUPMINNET(MINNETP)

END PROGRAM MAIN
