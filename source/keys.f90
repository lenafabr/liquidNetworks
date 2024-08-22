MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE


   ! -------- General program control ---------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE
  LOGICAL :: BNDDEBUG
  

  ! ----------------------
  ! Output / input
  ! -----------------------
  CHARACTER*100 :: OUTFILE, SNAPSHOTFILE, MINNETFILE, CONTFILE, MINNETOUTFILE
  LOGICAL :: DUMPSNAPSHOTS, RESTART, APPENDSNAPSHOTS
  INTEGER :: SNAPSHOTEVERY, STARTSNAPSLATE
  CHARACTER*100 :: TRACKFILE
  CHARACTER*100 :: PINFILE
  CHARACTER*100 :: BNDFILE
  CHARACTER*100 :: EBNDFILE
  LOGICAL :: TRACKNODES
  LOGICAL :: CONTINUERUN
  

  ! ------------
  ! network geometry and setup
  ! ------------
  INTEGER :: MINNETDIM ! predefined network dimension
  INTEGER :: MAXBRANCH ! max number of branches (per node) that can be input in the network
  ! fixed nodes
  INTEGER, PARAMETER :: MAXNFIX = 10000
  INTEGER :: NFIX
  INTEGER :: FIXNODES(MAXNFIX)
  LOGICAL :: RANDFIXNODES
  LOGICAL :: FIXNODEFROMNETFILE
  ! pinned nodes
  INTEGER :: NPIN ! just sets the initial population of pinned nodes (does not set the steady state density)
  INTEGER :: PINNODES(MAXNFIX)
  LOGICAL :: RANDPINNODES
  LOGICAL :: PINNODEFROMNETFILE
  ! active/inactive params
  INTEGER :: EXTRANODEFACT ! EXTRANODEFACT*NNODE gives size of NODEACT array, built in extras
  INTEGER :: EXTRAEDGEFACT ! same idea but for edges
  INTEGER :: EXTRATRACKFACT = 100
  LOGICAL :: CENTERENCLOSE, PERIODIC ! center and enclose the network in edges, or use PBC
  DOUBLE PRECISION :: RMULT, RCIRC, LBOX ! multiplier for boundary edges, radius of circle, or size of box
  DOUBLE PRECISION :: RNODES ! used by code to save RCIRC
  INTEGER :: NBNDEDGES ! optional user input number of boundary edges

  ! -----------------
  ! brownian dynamics
  ! -----------------
  DOUBLE PRECISION :: DELT, DIFF, MOB, DX ! timestep, diffusivity and mobility of nodes, length cutoff for edges
  DOUBLE PRECISION :: GROWTH, CATA ! growth, catastrophe rates
  DOUBLE PRECISION :: PIN, UNPIN ! pin and unpin rates, their ratio sets the steady state pin number
  DOUBLE PRECISION :: GVEL, FUSEPROB ! growth velocity and growth->fuse probability
  DOUBLE PRECISION :: GASTD ! grow angle standard deviation, assuming normal distribution of angles

  INTEGER :: BDSTEPS
  INTEGER :: BDPRINTEVERY
  LOGICAL :: DOBROWN
  LOGICAL :: DONETDYN

END MODULE KEYS
