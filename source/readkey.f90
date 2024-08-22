SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM, BOOLSTARTCENT, USERSETDX
  BOOLSTARTCENT=.FALSE.
  USERSETDX=.FALSE.

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.
  BNDDEBUG = .FALSE.

  ! input/output  
  OUTFILE = '*.out'
  DUMPSNAPSHOTS = .FALSE. ! periodically dump chain snapshots
  SNAPSHOTEVERY = 1 ! how often to dump snapshots
  STARTSNAPSLATE = 0
  SNAPSHOTFILE = '*.snap.out' ! snapshot file
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  TRACKFILE = '*.track.snap.out' ! snapshot file for tracking nodes
  PINFILE = '*.pin.snap.out' ! snapshot file for pinned nodes
  BNDFILE = '*.bnd.snap.out' ! snapshot file for boundary nodes
  EBNDFILE = '*.ebnd.snap.out' ! snapshot file for boundary edges
  TRACKNODES = .FALSE. ! default do not output tracking nodes file
  CONTINUERUN = .FALSE.

  ! network parameters
  MINNETFILE = '*.net' ! File containing network structure
  MINNETOUTFILE='*.netout.out' ! File for outputting intermediate network structure
  MAXBRANCH = 3 ! max number of branches per node
  ! if positive, predefines a network dimension
  ! otherwise, dimension set to # items in NODE row of network file - 2
  MINNETDIM = 2
  NFIX = 0 ! number of fixed nodes
  FIXNODEFROMNETFILE = .FALSE. ! determine fixed nodes based on network file
  RANDFIXNODES = .FALSE.
  NPIN = 0 ! number of pinned nodes
  PINNODEFROMNETFILE = .FALSE. ! determine pinned nodes based on network file
  RANDPINNODES = .FALSE.
  EXTRANODEFACT = 100 ! by default length(nodeact) = 100 * nnode
  EXTRAEDGEFACT = 100 ! same but for edgeact and nnedge
  CENTERENCLOSE = .FALSE. ! if true center and enclose network with circular boundary
  RMULT = 2D0 ! sets radial boundary at 2*(max radial distance of initial node)
  RCIRC = 1D8 ! this is set during center and enclose routine
  RNODES = 1D8
  NBNDEDGES = 0 ! default zero, have number set by edgelen of bndary edges
  
  
  ! brownian dynamics
  DELT = 1D-3 ! time step
  BDSTEPS = 1000 ! number of brownian steps to run
  BDPRINTEVERY = 1 ! how often to print output
  DOBROWN = .TRUE. ! include brownian forces
  DONETDYN = .TRUE. ! include brownian forces
  DIFF = 1.0D-3 ! diffusivity
  MOB = 1.0 ! drift coefficient
  GVEL = 1.0 ! growth velocity
  FUSEPROB = 1.0 ! probability of fusion upon growth intersecting an edge, between 0 and 1
  GROWTH = 1.0 ! rate of growth events. Units of 1/s
  CATA = 1D-8 ! rate of catastrophe post growth. Units of 1/s
  PIN = 0.0 ! rate of pinning events. Units of 1/s
  UNPIN = 0.0 ! rate of unpinning events. Units of 1/node/s
  DX = SQRT(DIFF*DELT)! max step size/minimum edgelength (set by either sqrt(D*dt) or B*dt, whichever is larger)
  GASTD = 1.0D-2 ! width of normal distribution for angles of growth, default to basically perpendicular growth
  
  EXTRATRACKFACT = 1000 ! how much extra space for storing track information

  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
    NPARAMFILES = 1
    PARAMFILES(1) = 'param'
    ARG = ''
  ELSE
    DO I = 1,NUMARG
      CALL GETARG(I, ARG)
      NPARAMFILES = NPARAMFILES + 1
      WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
      PARAMFILES(NPARAMFILES) = DUMSTR
    ENDDO
    ! reset arg to its original value
    IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
    NPARAMREAD = NPARAMREAD + 1

    PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
    INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
    IF (.NOT.LDUM) THEN
      PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
      STOP 1
    ENDIF
    OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

    ! read in the keywords one line at a time
    DO 
      CALL READLINE(PF,FILEEND,NITEMS)
      IF (FILEEND.and.nitems.eq.0) EXIT

      ! skip empty lines
      IF (NITEMS.EQ.0) CYCLE

      ! Read in the keyword for this line
      CALL READA(WORD,CASESET=1)

      ! Skip any empty lines or any comment lines
      IF (WORD(1:1).EQ.'#') CYCLE

      SELECT CASE(WORD) ! pick which keyword
      CASE('ACTION')
        CALL READA(ACTION, CASESET=1)
      CASE('BDPRINTEVERY')
        CALL READI(BDPRINTEVERY)   
      CASE('BDSTEPS')
        CALL READI(BDSTEPS)                   
      CASE('BNDDEBUG')
        CALL READO(BNDDEBUG)
      CASE('BNDFILE')
        CALL READA(BNDFILE)
      CASE('CATA')
        CALL READF(CATA)
      CASE('CENTERENCLOSE')
        CENTERENCLOSE = .TRUE.
      CASE('DELT')
        CALL READF(DELT)        
      CASE('DIFF')
        CALL READF(DIFF)
      CASE('DX')
        USERSETDX=.TRUE.
        CALL READF(DX)
      CASE('EBNDFILE')
        CALL READA(EBNDFILE)
      CASE('EXTRAEDGEFACT')
        CALL READI(EXTRAEDGEFACT)
      CASE('EXTRANODEFACT')
        CALL READI(EXTRANODEFACT)
      CASE('EXTRATRACKFACT')
        CALL READI(EXTRATRACKFACT)
      CASE('FIXNODE')
        NFIX = NFIX + 1
      !   IF (NFIX(FC).GT.MAXNABSORBER) THEN
      !      PRINT*, 'ERROR IN READKEY: too many fixed nodes', NFIX, MAXNABSORBER
      !      STOP 1
      !   ENDIF
        CALL READI(FIXNODES(NFIX))
      CASE('FIXNODEFROMNETFILE')
        IF (NITEMS.GT.1) THEN
          CALL READO(FIXNODEFROMNETFILE)
        ELSE
          FIXNODEFROMNETFILE = .TRUE.
        ENDIF
      CASE('FUSEPROB')
        CALL READF(FUSEPROB)
      CASE('GASTD')
        CALL READF(GASTD)
      CASE('GROWTH')
        CALL READF(GROWTH)
      CASE('GVEL')
        CALL READF(GVEL)
      CASE('MAXBRANCH')
        CALL READI(MAXBRANCH)
      CASE('MINNETDIM')
        CALL READI(MINNETDIM)
      CASE('MINNETFILE')
        CALL READA(MINNETFILE)
        IF (NITEMS.GT.2) CALL READA(MINNETOUTFILE)
      CASE('MOB')
        CALL READF(MOB)
      CASE('NBNDEDGES')
        CALL READI(NBNDEDGES)
      CASE('NOBROWN')
        DOBROWN = .FALSE.
      CASE('NONETDYN')
        DONETDYN = .FALSE.
      CASE('NFIX')
        CALL READI(NFIX)
      CASE('NPIN')
        CALL READI(NPIN)
      CASE('OUTFILE')
        CALL READA(OUTFILE)
      CASE('PIN')
        CALL READF(PIN)
      CASE('PINFILE')
        CALL READA(PINFILE)
      CASE('PINNODE')
        NPIN = NPIN + 1
      !   IF (NPIN(FC).GT.MAXNABSORBER) THEN
      !      PRINT*, 'ERROR IN READKEY: too many fixed nodes', NPIN, MAXNABSORBER
      !      STOP 1
      !   ENDIF
        CALL READI(PINNODES(NPIN))
      CASE('PINNODEFROMNETFILE')
        IF (NITEMS.GT.1) THEN
          CALL READO(PINNODEFROMNETFILE)
        ELSE
          PINNODEFROMNETFILE = .TRUE.
        ENDIF
      CASE('RANDFIXNODES')
        IF (NITEMS.GT.1) THEN
          CALL READO(RANDFIXNODES)
        ELSE
          RANDFIXNODES = .TRUE.
        ENDIF
      CASE('RANDPINNODES')
        IF (NITEMS.GT.1) THEN
          CALL READO(RANDPINNODES)
        ELSE
          RANDPINNODES = .TRUE.
        ENDIF
      CASE('RMULT')
        CALL READF(RMULT)
      CASE('RNGSEED')
        CALL READI(RNGSEED)
      CASE('SNAPSHOTFILE')
        CALL READA(SNAPSHOTFILE)
      CASE('SNAPSHOTS')
        DUMPSNAPSHOTS = .TRUE.
        IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
        IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
        IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
      CASE('STARTSNAPSLATE')
        CALL READI(STARTSNAPSLATE)
      CASE('TRACKNODES')
        CALL READO(TRACKNODES)
      CASE('TRACKFILE')
        CALL READA(TRACKFILE)
      CASE('UNPIN')
        CALL READF(UNPIN)
      CASE('VERBOSE')
        CALL READO(VERBOSE)
      CASE DEFAULT
        print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
      END SELECT
    ENDDO
    CLOSE(PF)
  ENDDO

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  

  IF ((FUSEPROB.LT.0D0).OR.(FUSEPROB.GT.1D0)) THEN
    PRINT*, 'ERROR IN FUSEPROB, must be real between 0 and 1'
    STOP 1
  ENDIF

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(TRACKFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(PINFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(BNDFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(EBNDFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(MINNETFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(MINNETOUTFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
    ! use the current time of day in milliseconds
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
    ! use the last 5 characters in the command-line argument
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
    ! use the last 4 characters in the command-line argument 
    ! and additionally the millisecond time 
    CALL DATE_AND_TIME(VALUES=TIMEVAL)
    SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
    ! use this seed directly
    SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Network file: ', TRIM(MINNETFILE)
  print*, 'Output file: ', TRIM(OUTFILE)  
  IF (DUMPSNAPSHOTS) THEN
    PRINT*, 'Dumping net snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
    IF (TRACKNODES) THEN
      PRINT*, 'Dumping node tracking every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(TRACKFILE))
    ENDIF
    PRINT*, 'Dumping node pins every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(PINFILE))
    IF (BNDDEBUG) THEN
      PRINT*, 'Dumping boundary nodes every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(BNDFILE))
      PRINT*, 'Dumping boundary edges every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(EBNDFILE))
    ENDIF
  ENDIF

  ! set DX based on input values for B, D, unless user specified a DX
  IF (.NOT.USERSETDX) THEN
    IF (DIFF.GE.((MOB**2)*DELT)) THEN
      DX = SQRT(DIFF*DELT)
    ELSE
      DX = MOB*DELT
    ENDIF
  ENDIF

  PRINT*, '----------------------------------------------------'
  

END SUBROUTINE READKEY
