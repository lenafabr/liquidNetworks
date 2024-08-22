MODULE MINNETUTIL
  ! Utilities for dealing with connectivity and geometry of a network
  ! including defining, manipulating, input/output

  USE STACKUTIL, ONLY : STACK
  USE KEYS, ONLY : VERBOSE

  IMPLICIT NONE
  
  TYPE MINNET
     ! define a network structure (connectivity and geometry)

     ! dimension of the space the network is embedded in
     INTEGER :: DIM 

     ! ----------------------
     ! information on network nodes
     ! ----------------------
     INTEGER :: NNODE ! number of nodes
     ! list of other node indices each node connects to
     INTEGER, POINTER :: NODENODE(:,:)
     ! degree (number of branches) of each node
     INTEGER, POINTER :: NODEDEG(:)
     ! list of branch indices each node connects to
     INTEGER, POINTER :: NODEEDGE(:,:)          
     ! spatial location of node
     DOUBLE PRECISION, POINTER :: NODEPOS(:,:)
     ! which nodes are fixed in space
     LOGICAL, POINTER :: NODEFIX(:)
     ! which nodes are active
     LOGICAL, POINTER :: NODEACT(:)
     ! how many nodes are active
     INTEGER :: NACTIVE
     ! stack of inactive nodes
     TYPE(STACK),POINTER :: NODESTACK
     ! index of highest active node
     INTEGER :: NODEHIGH
     ! which nodes are growing out of an edge
     ! 1 for growing, 2 for receding, 0 for not growing
     INTEGER, POINTER :: NODEGROW(:)
     ! growth direction
     DOUBLE PRECISION, POINTER :: NODEDIR(:,:)
     ! past node positions
     DOUBLE PRECISION, POINTER :: NODEPOSP(:,:)
     ! which nodes are on boundary
     LOGICAL, POINTER :: NODEBND(:)
     ! node tracking
     INTEGER, POINTER :: NODETRACK(:)
     ! how many nodes are being tracked
     INTEGER :: NODEON
     ! which nodes are pinned (must also be fixed)
     LOGICAL, POINTER :: NODEPIN(:)
     
     ! ------------------
     ! information on network branches
     ! ------------------
     INTEGER :: NEDGE ! Number of branches
     ! nodes at the start and end of each branch
     INTEGER, POINTER :: EDGENODE(:,:)
     ! spatial position of branch starting points
     ! branch direction and length
     ! branch direction is normalized by length
     DOUBLE PRECISION, POINTER :: EDGESTART(:,:), EDGEDIR(:,:), EDGELEN(:)
     ! which edges are active
     LOGICAL, POINTER :: EDGEACT(:)
     ! how many edges are active
     INTEGER :: EACTIVE
     ! stack of inactive edges
     TYPE(STACK),POINTER :: EDGESTACK
     ! index of highest active edge
     INTEGER :: EDGEHIGH
     ! index of longest edge
     INTEGER :: EIDMAX
     ! boundary edge marker
     LOGICAL, POINTER :: EDGEBND(:)
     ! past edgenode info
     INTEGER, POINTER :: EDGENODEP(:,:)
     
          
     ! ARRAYSET: Have arrays been allocated?
     ! STRUCTURESET: has node and edge information been set up?
     LOGICAL :: ARRAYSET = .FALSE., STRUCTURESET = .FALSE.
     
  END TYPE MINNET

CONTAINS
  

  SUBROUTINE INSERTINTERNODE(MINNETP,EDGE,FRAC,NID,EIDIN,WIPEPREV)
    ! insert intermediate node along an edge
    ! reconnecting everything appropriately
    ! place new node in index NID, new edge in index EIDIN
    ! new node will be along edge EDGE, a fraction FRAC of the way along it
    ! WIPEPREV: if true, then the node NID was already set up
    ! as a degree 2 node and needs to be removed from the network
    ! in this case, use the 2nd edge attached to it instead of EIDIN
    ! So when WIPEPREV is true, the value of EIDIN does not matter
    ! IMPORTANT: Must set fixed/growing node fields after function call
    
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: EDGE, NID, EIDIN
    DOUBLE PRECISION, INTENT(IN) :: FRAC
    LOGICAL, INTENT(IN) :: WIPEPREV

    INTEGER :: EID, NODE1, NODE2, NC, IND
    DOUBLE PRECISION :: LEN1, LEN2, POS(MINNETP%DIM)

    IF (FRAC.LT.-TINY(0D0).OR.FRAC.GT.1D0+TINY(1D0)) THEN
       PRINT*, 'ERROR IN INSERTINTERNODE: frac out of range', FRAC
       STOP 1
    ENDIF
    
    IF (WIPEPREV) THEN
       IF (MINNETP%NODEDEG(NID).NE.2) THEN
          PRINT*, 'ERROR IN INSERTINTERNODE: previous node does not have degree 2'
          STOP 1
       ENDIF
       PRINT*, 'ERROR IN INSERTINTERNODE: wipeprev not set up yet'
       STOP 1
    ELSE
       EID = EIDIN
    ENDIF

    ! edge where contraction point is       
    NODE1 = MINNETP%EDGENODE(EDGE,1) ! first node for this edge
    NODE2 = MINNETP%EDGENODE(EDGE,2) ! second node for this edge
       
    ! length of edge up to contraction and after it
    LEN1 = FRAC*MINNETP%EDGELEN(EDGE)
    LEN2 = MINNETP%EDGELEN(EDGE)-LEN1

    ! new node position in space       
    POS = MINNETP%EDGESTART(EDGE,:) &
            & + LEN1*MINNETP%EDGEDIR(EDGE,:)

    MINNETP%NODEPOS(NID,:) = POS
    ! if node already active, assume it was degree 1 growth/retraction->fuse event
    IF (MINNETP%NODEACT(NID)) THEN
      IF (VERBOSE) THEN
        PRINT*, 'already active'
      ENDIF
      MINNETP%NODEDEG(NID) = 3
      ! connect this node to the nodes belonging to this edge
      MINNETP%NODENODE(NID,2:3) = MINNETP%EDGENODE(EDGE,:)
      MINNETP%NODEEDGE(NID,2:3) = (/EDGE,EID/)
    ELSE
      MINNETP%NODEDEG(NID) = 2; ! extra node is inserted into edge, deg=2
      ! connect this node to the nodes belonging to this edge
      MINNETP%NODENODE(NID,1:2) = MINNETP%EDGENODE(EDGE,:)
      MINNETP%NODEEDGE(NID,1:2) = (/EDGE,EID/)
    ENDIF

    ! break up edge into 2. Current EDGE index keeps the first half
    ! EID index gets the second half
    MINNETP%EDGELEN(EDGE) = LEN1   
    MINNETP%EDGENODE(EDGE,2) = NID

    MINNETP%EDGELEN(EID) = LEN2
    MINNETP%EDGENODE(EID,:) = (/NID,NODE2/)
    MINNETP%EDGESTART(EID,:) = POS
    MINNETP%EDGEDIR(EID,:) = MINNETP%EDGEDIR(EDGE,:)
    
    IND = FINDLOC(MINNETP%NODEEDGE(NODE1,:),EDGE,DIM=1)
    MINNETP%NODENODE(NODE1,IND) = NID

    IND = FINDLOC(MINNETP%NODEEDGE(NODE2,:),EDGE,DIM=1)
    MINNETP%NODENODE(NODE2,IND) = NID
    MINNETP%NODEEDGE(NODE2,IND) = EID

    IF (WIPEPREV) THEN
    ELSEIF (MINNETP%NODEACT(NID)) THEN
      IF (VERBOSE) THEN
        PRINT*, 'already active 2'
      ENDIF
      ! node already active, no need to do any updating of nodes, but update edges
      MINNETP%EDGEACT(EID) = .TRUE.
      MINNETP%EACTIVE = MINNETP%EACTIVE + 1
      IF (EID.GT.MINNETP%EDGEHIGH) THEN
        MINNETP%EDGEHIGH = EID
      ENDIF
      MINNETP%EDGEBND(EID) = .FALSE.
    ELSE
      ! node was not already active, update active edges/nodes
      MINNETP%NODEACT(NID) = .TRUE.
      MINNETP%NACTIVE = MINNETP%NACTIVE + 1
      IF (NID.GT.MINNETP%NODEHIGH) THEN
        MINNETP%NODEHIGH = NID
      ENDIF

      MINNETP%EDGEACT(EID) = .TRUE.
      MINNETP%EACTIVE = MINNETP%EACTIVE + 1
      IF (EID.GT.MINNETP%EDGEHIGH) THEN
        MINNETP%EDGEHIGH = EID
      ENDIF
      MINNETP%EDGEBND(EID) = .FALSE.
    ENDIF

    ! PRINT*, 'inserted node', NID, 'at position', MINNETP%NODEPOS(NID,:)
    ! PRINT*, 'from splitting edge: ', EDGE
    ! PRINT*, 'between neighbors:', NODE1, NODE2


  END SUBROUTINE INSERTINTERNODE

  

  SUBROUTINE DOUBLEINSERTINTERNODE(MINNETP,EC,T1,T2,NID,NID2,EIDIN,EIDIN2)
    ! calls insertinternode twice, properly scales fractions
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: EC, NID, NID2, EIDIN, EIDIN2
    DOUBLE PRECISION, INTENT(IN) :: T1, T2

    IF (T1.GT.T2) THEN
      CALL INSERTINTERNODE(MINNETP,EC,T1,NID,EIDIN,.FALSE.)
      CALL INSERTINTERNODE(MINNETP,EC,T2/T1,NID2,EIDIN2,.FALSE.)
    ELSEIF (T2.GT.T1) THEN
      CALL INSERTINTERNODE(MINNETP,EC,T2,NID,EIDIN,.FALSE.)
      CALL INSERTINTERNODE(MINNETP,EC,T1/T2,NID2,EIDIN2,.FALSE.)
    ELSEIF (T1.EQ.T2) THEN
      PRINT*, 'ERROR IN DOUBLEINSERTINTERNODE: T1=T2'
      STOP 1
    ENDIF

    ! set these by default makes life eaiser in partutil, only place this is used
    ! node not growing
    MINNETP%NODEGROW(NID2) = 0
    ! this node is not fixed
    MINNETP%NODEFIX(NID2) = .FALSE.
    ! set previous position
    MINNETP%NODEPOSP(NID2,:) = MINNETP%NODEPOS(NID2,:)
    MINNETP%EDGENODEP(EC,:) = MINNETP%EDGENODE(EC,:)
    MINNETP%EDGENODEP(EIDIN2,:) = MINNETP%EDGENODE(EIDIN2,:)

  END SUBROUTINE DOUBLEINSERTINTERNODE


  
  SUBROUTINE SETEDGELENSSIMPLE(MINNETP)
    ! sets edge lengths based off of new node positions
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER :: EID, NODE1, NODE2

    DO EID = 1,MINNETP%EDGEHIGH
      IF (MINNETP%EDGEACT(EID)) THEN
        ! find indices of bounding nodes
        NODE1 = MINNETP%EDGENODE(EID,1)
        NODE2 = MINNETP%EDGENODE(EID,2)

        ! reset EDGESTART as new node 1 position and EDGEDIR as difference of positions
        MINNETP%EDGESTART(EID,:) = MINNETP%NODEPOS(NODE1,:)
        MINNETP%EDGEDIR(EID,:) = MINNETP%NODEPOS(NODE2,:)-MINNETP%NODEPOS(NODE1,:)

        ! EDGELEN next
        MINNETP%EDGELEN(EID) = SQRT(SUM( MINNETP%EDGEDIR(EID,:)*MINNETP%EDGEDIR(EID,:) ))

        ! go back and normalize edgedir
        MINNETP%EDGEDIR(EID,:) = MINNETP%EDGEDIR(EID,:)/MINNETP%EDGELEN(EID)

        IF (MINNETP%NODEBND(NODE1).AND.MINNETP%NODEBND(NODE2).AND.(.NOT.MINNETP%EDGEBND(EID))) THEN
          MINNETP%EDGEBND(EID) = .TRUE.
        ELSEIF ((.NOT.MINNETP%NODEBND(NODE1)).AND.(.NOT.MINNETP%NODEBND(NODE2)).AND.MINNETP%EDGEBND(EID)) THEN
          MINNETP%EDGEBND(EID) = .FALSE.
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE SETEDGELENSSIMPLE



  SUBROUTINE MOVENODE(MINNETP,NODE,NEWPOS)
    ! sets edge lengths based off of new node positions
    ! should only be used during T1 code for d3 nodes
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER :: EC, EID, NODE1, NODE2, DEG, EDGES(MAXBRANCH)
    INTEGER, INTENT(IN) :: NODE
    DOUBLE PRECISION, INTENT(IN) :: NEWPOS(MINNETP%DIM)

    ! first move the node
    MINNETP%NODEPOS(NODE,:) = NEWPOS
    
    ! reset edge directions and lengths for connected edges
    DEG = MINNETP%NODEDEG(NODE)
    EDGES = MINNETP%NODEEDGE(NODE,1:DEG)

    DO EC = 1,DEG
      EID = EDGES(EC)
      ! find indices of bounding nodes
      NODE1 = MINNETP%EDGENODE(EID,1)
      NODE2 = MINNETP%EDGENODE(EID,2)

      ! reset EDGESTART as new node 1 position and EDGEDIR as difference of positions
      MINNETP%EDGESTART(EID,:) = MINNETP%NODEPOS(NODE1,:)
      MINNETP%EDGEDIR(EID,:) = MINNETP%NODEPOS(NODE2,:)-MINNETP%NODEPOS(NODE1,:)

      ! EDGELEN next
      MINNETP%EDGELEN(EID) = SQRT(SUM( MINNETP%EDGEDIR(EID,:)*MINNETP%EDGEDIR(EID,:) ))

      ! go back and normalize edgedir
      MINNETP%EDGEDIR(EID,:) = MINNETP%EDGEDIR(EID,:)/MINNETP%EDGELEN(EID)
    ENDDO

  END SUBROUTINE MOVENODE


  
  SUBROUTINE SETEDGELENS(MINNETP,SHORTEDGES,LTOT)
    ! sets edge lengths based off of new node positions
    ! input:
    !   MINNETP -- pointer to minimal network object
    ! for now, loops through all edges and outputs:
    !   SHORTEDGES - nedge+1 array, first slot keeps track of how many edges are
    !                below minimum cutoff, then next slots give edge numbers
    !   LTOT       - total length
    USE KEYS, ONLY : DX
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(OUT) :: SHORTEDGES(MINNETP%NEDGE+1)
    DOUBLE PRECISION, INTENT(OUT) :: LTOT
    DOUBLE PRECISION :: LMAX
    INTEGER :: EID, NODE1, NODE2

    SHORTEDGES(1) = 0
    LTOT = 0D0

    ! update longest edge length
    LMAX = 0D0

    DO EID = 1,MINNETP%EDGEHIGH
      IF (MINNETP%EDGEACT(EID)) THEN
        ! find indices of bounding nodes
        NODE1 = MINNETP%EDGENODE(EID,1)
        NODE2 = MINNETP%EDGENODE(EID,2)

        ! reset EDGESTART as new node 1 position and EDGEDIR as difference of positions
        MINNETP%EDGESTART(EID,:) = MINNETP%NODEPOS(NODE1,:)
        MINNETP%EDGEDIR(EID,:) = MINNETP%NODEPOS(NODE2,:)-MINNETP%NODEPOS(NODE1,:)   

        ! EDGELEN next
        MINNETP%EDGELEN(EID) = SQRT(SUM( MINNETP%EDGEDIR(EID,:)*MINNETP%EDGEDIR(EID,:) ))

        ! go back and normalize edgedir
        MINNETP%EDGEDIR(EID,:) = MINNETP%EDGEDIR(EID,:)/MINNETP%EDGELEN(EID)

        ! add to total edge length
        LTOT = LTOT + MINNETP%EDGELEN(EID)

        ! check if longer than LMAX
        IF (MINNETP%EDGELEN(EID).GT.LMAX) THEN
          LMAX = MINNETP%EDGELEN(EID)
          MINNETP%EIDMAX = EID
        ENDIF
        
        ! update SHORTEDGES list if edgelength is less than DX
        IF (MINNETP%EDGELEN(EID).LT.DX) THEN
          SHORTEDGES(1) = SHORTEDGES(1) + 1
          SHORTEDGES(1+SHORTEDGES(1)) = EID
        ENDIF
      ENDIF

    ENDDO

  END SUBROUTINE SETEDGELENS



  SUBROUTINE UPDATEPREV(MINNETP)
    ! update edgenodep and nodeposp, just transfer over current edgenode and nodepos
    ! later, when rearrangements, or growths, or deletions occur make sure to update
    ! these values accordingly
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP

    ! that's it...
    MINNETP%EDGENODEP(:,:) = MINNETP%EDGENODE(:,:)
    MINNETP%NODEPOSP(:,:) = MINNETP%NODEPOS(:,:)

  END SUBROUTINE UPDATEPREV



  SUBROUTINE REARRANGE(MINNETP,SHORTEDGES)
    ! Given shortedges, begins to look for possible rearrangements that would shorten total length
    USE KEYS, ONLY : MAXBRANCH, DX
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: SHORTEDGES(MINNETP%NEDGE+1)

    INTEGER :: EID, NODE1, NODE2, TMPNODE, NC, EC, CT, ECC
    INTEGER :: DEG1, DEG2
    DOUBLE PRECISION :: TOTLEN, NEWLEN1, NEWLEN2, DIST

    IF (VERBOSE) THEN
      PRINT*, 'Looping through shortedges'
      PRINT*, 'shortedges', SHORTEDGES(1:(SHORTEDGES(1)+1))
    ENDIF
    DO EC = 1,SHORTEDGES(1)
      EID = SHORTEDGES(EC+1)
      ! checking to make sure previous shortedge didn't deactivate/lengthen this one
      IF ((MINNETP%EDGEACT(EID)).AND.(MINNETP%EDGELEN(EID).LT.DX)) THEN
        ! find indices of bounding nodes
        NODE1 = MINNETP%EDGENODE(EID,1)
        NODE2 = MINNETP%EDGENODE(EID,2)
        IF (VERBOSE) THEN
          PRINT*, 'NODE1', NODE1
          PRINT*, 'NODE2', NODE2
        ENDIF
   
        ! find their degree
        DEG1 = MINNETP%NODEDEG(NODE1)
        DEG2 = MINNETP%NODEDEG(NODE2)

   
        ! see if this is a Pa, T1, or C
        IF ((DEG1.EQ.2).AND.(DEG2.EQ.3)) THEN
          CALL PASSEVENT(MINNETP,NODE1,NODE2)

        ELSEIF ((DEG1.EQ.3).AND.(DEG2.EQ.2)) THEN
          CALL PASSEVENT(MINNETP,NODE2,NODE1)

        ELSEIF ((DEG1.EQ.3).AND.(DEG2.EQ.3)) THEN
          CALL T1(MINNETP,NODE1,NODE2)

        ELSEIF (DEG1.EQ.1) THEN
          CALL CATASTROPHE(MINNETP,NODE1)

        ELSEIF (DEG2.EQ.1) THEN
          CALL CATASTROPHE(MINNETP,NODE2)

        ELSEIF ((DEG2.EQ.2).AND.(DEG1.EQ.2)) THEN
          ! merge nodes
          CALL D2ABSORBD2(MINNETP,NODE1,NODE2)

        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE REARRANGE



  SUBROUTINE CATASTROPHE(MINNETP,NODE1)
    ! node1 is degree 1, check if this is a full retraction (ie. node1 is currently shrinking and near its neighbor)
    ! if so, remove the node
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER :: NODE1

    ! if currently retracting
    IF (MINNETP%NODEGROW(NODE1).EQ.2) THEN
      CALL REMOVED1NODE(MINNETP,NODE1)
    ENDIF


  END SUBROUTINE CATASTROPHE



  SUBROUTINE FINDMAXEDGE(MINNETP)
    ! simple subroutine to loop throuh all edges and find EIDMAX
    USE KEYS, ONLY : RCIRC
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    DOUBLE PRECISION :: LMAX
    INTEGER :: EID, N1, N2

    LMAX = 0D0

    DO EID=1,MINNETP%EDGEHIGH
      ! check if longer than LMAX
      IF ((MINNETP%EDGELEN(EID).GT.LMAX).AND.(MINNETP%EDGEACT(EID))) THEN
        LMAX = MINNETP%EDGELEN(EID)
        MINNETP%EIDMAX = EID
      ENDIF
    ENDDO

  END SUBROUTINE FINDMAXEDGE


  
  SUBROUTINE FUSECHECK(MINNETP,GROWNODE,FUSED)
    ! new new version: 
    !  1) check for nearby nodes, etc,etc
    !  2) check for intersections using LLSLIDEINT or LLSLIDEINT2, essentially slides the growing node
    !     along its path, while nearby edges slide from initial position to final position, interpolating
    !     all movement as linear, and checking for hits of growing node and edge
    !  3) function returns
    !       a) T, (between 0 and 1) which is how far along edge the growing node hits and
    !       b) TSTAR, time (between 0 and 1) at which this hit occurs (useful when multiple hits)
    USE KEYS, ONLY : GVEL, DELT, FUSEPROB, DX, VERBOSE, TRACKNODES, RCIRC,&
                    & PERIODIC, CENTERENCLOSE
    USE mt19937, ONLY : GRND
    USE GENUTIL, ONLY : LLINTERSECT, LLSLIDEINT, LLSLIDEINT2
    USE STACKUTIL, ONLY : POP
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: GROWNODE
    LOGICAL, INTENT(OUT) :: FUSED
    DOUBLE PRECISION :: LMAX, DIST, RN, T, T2, TSTAR, TSTAR2, DMAX, MDXND(MINNETP%DIM)
    INTEGER :: NC, EID, EIDIN, EC, DEG, N1, N2, NEWNODE, IND, N1P, N2P
    INTEGER :: CROSS, NCROSS, CROSSES(MINNETP%EDGEHIGH)
    DOUBLE PRECISION :: TS(MINNETP%EDGEHIGH), TSTARS(MINNETP%EDGEHIGH), POS(MINNETP%DIM), RPOS

    FUSED=.FALSE.

    POS = MINNETP%NODEPOS(GROWNODE,:)
    RPOS = SQRT(SUM(POS**2))

    ! max edgelength
    LMAX = MINNETP%EDGELEN(MINNETP%EIDMAX)
    ! max possible distance
    DMAX = SQRT((LMAX/2)**2 + (GVEL*DELT)**2 + DX**2)
    NCROSS = 0
    CROSSES = 0
    TS = 0D0
    TSTARS = 0D0

    ! first loop through nodes and make list of crossings
    DO NC=1,MINNETP%NODEHIGH 

      IF ((MINNETP%NODEACT(NC)).AND.(NC.NE.GROWNODE)) THEN
        ! active node, check if nodepos within certain distance
        CALL NODENODEDIST(MINNETP,NC,GROWNODE,DIST)

        IF ((DIST.LE.DMAX).AND.(NC.NE.MINNETP%NODENODE(GROWNODE,1))) THEN
          ! node within minimum distance, check each of its edges for crossings
          DEG = MINNETP%NODEDEG(NC)

          DO EC=1,DEG
            CROSS=0
            EID = MINNETP%NODEEDGE(NC,EC)
            N1 = MINNETP%EDGENODE(EID,1)
            N2 = MINNETP%EDGENODE(EID,2)
            N1P = MINNETP%EDGENODEP(EID,1)
            N2P = MINNETP%EDGENODEP(EID,2)

            ! if N1P=N2, swap N1P and N2P
            IF (N1P.EQ.N2) THEN
              N1P = N2P
              N2P = N2
            ELSEIF (N2P.EQ.N1) THEN
              N2P = N1P
              N1P = N1
            ENDIF

            ! if intersecting edge is stationary, call different function
            IF (((MINNETP%NODEFIX(N1)).AND.(MINNETP%NODEFIX(N2))).OR. &
             &((MINNETP%NODEPIN(N1)).AND.(MINNETP%NODEPIN(N2))).OR. &
             &((MINNETP%NODEFIX(N1)).AND.(MINNETP%NODEPIN(N2))).OR. &
             &((MINNETP%NODEPIN(N1)).AND.(MINNETP%NODEFIX(N2))))  THEN
              CALL LLSLIDEINT2(MINNETP%NODEPOS(N1,1),MINNETP%NODEPOS(N1,2), &
                 & MINNETP%NODEPOS(N2,1),MINNETP%NODEPOS(N2,2), &
                 & MINNETP%NODEPOSP(GROWNODE,1),MINNETP%NODEPOSP(GROWNODE,2), &
                 & MINNETP%NODEPOS(GROWNODE,1),MINNETP%NODEPOS(GROWNODE,2), &
                 & CROSS, T, TSTAR)
            ELSE
              CALL LLSLIDEINT(MINNETP%NODEPOSP(N1P,1),MINNETP%NODEPOSP(N1P,2), &
                 & MINNETP%NODEPOSP(N2P,1),MINNETP%NODEPOSP(N2P,2), &
                 & MINNETP%NODEPOSP(GROWNODE,1),MINNETP%NODEPOSP(GROWNODE,2), &
                 & MINNETP%NODEPOS(N1,1),MINNETP%NODEPOS(N1,2), &
                 & MINNETP%NODEPOS(N2,1),MINNETP%NODEPOS(N2,2), &
                 & MINNETP%NODEPOS(GROWNODE,1),MINNETP%NODEPOS(GROWNODE,2), &
                 & CROSS, T, TSTAR, T2, TSTAR2)
            ENDIF

            ! if crossed, and haven't already documented this cross before save it
            IF (CROSS.EQ.1) THEN
              IF (T.GT.0D0) THEN ! don't allow nodes to be placed on top of eachother
                IF (.NOT.ANY(CROSSES.EQ.EID)) THEN
                  NCROSS = NCROSS + 1
                  CROSSES(NCROSS) = EID
                  TS(NCROSS) = T
                  TSTARS(NCROSS) = TSTAR
                  IF (VERBOSE) THEN
                    PRINT*, GROWNODE, 'crossing event! with node ', NC, ' neighbor ', EC, 'edge,', EID
                    PRINT*, 'N1, N2, T, TSTAR', N1, N2, T, TSTAR
                  ENDIF
                ENDIF
              ENDIF

            ELSEIF (CROSS.EQ.2) THEN
              IF (T.GT.0D0) THEN ! don't allow nodes to be placed on top of eachother
                IF (.NOT.ANY(CROSSES.EQ.EID)) THEN
                  ! save first crossing
                  NCROSS = NCROSS + 1
                  CROSSES(NCROSS) = EID
                  TS(NCROSS) = T
                  TSTARS(NCROSS) = TSTAR

                  ! include second crossing as well
                  NCROSS = NCROSS + 1
                  CROSSES(NCROSS) = EID
                  TS(NCROSS) = T2
                  TSTARS(NCROSS) = TSTAR2
                  IF (VERBOSE) THEN
                    PRINT*, GROWNODE, 'double crossing event! with node ', NC, ' neighbor ', EC, 'edge,', EID
                    PRINT*, 'N1, N2, T, TSTAR, T2, TSTAR2', N1, N2, T, TSTAR, T2, TSTAR2
                  ENDIF
                ENDIF
              ENDIF
            ENDIF

          ENDDO
        ENDIF
      ENDIF
    ENDDO

    ! loop through crossing events, choose the one that is farthest along growing edge
    ! (ie. would have hit first). Sample fusion probability, if doesn't fuse, go to next
    ! farthest crossing event and sample
    DO WHILE (NCROSS.GT.0)

      ! find minimum of TSTARS, this is the 'first' intersection in this time step
      IND = MINLOC(TSTARS(1:NCROSS), DIM=1)

      ! which edge is it crossing?
      EID = CROSSES(IND)
      RN = GRND()

      ! if RN less than fuseprob or if edge is boundary edge, fuse!
      IF ((RN.LE.FUSEPROB).OR.(MINNETP%EDGEBND(EID))) THEN
        T = TS(IND)
        IF (VERBOSE) THEN
          PRINT*, 'it stuck! with edge, (BND?) ', EID, MINNETP%EDGEBND(EID)
          PRINT*, 'nodes, ', MINNETP%EDGENODE(EID,1), MINNETP%EDGENODE(EID,2)
          PRINT*, 'bnd nodes? ', MINNETP%NODEBND(MINNETP%EDGENODE(EID,1)), MINNETP%NODEBND(MINNETP%EDGENODE(EID,2))
        ENDIF

        ! with fusion, need another free edge
        CALL POP(MINNETP%EDGESTACK,EIDIN)

        ! fusion, GROWNODE already active
        CALL INSERTINTERNODE(MINNETP,EID,T,GROWNODE,EIDIN,.FALSE.)

        IF (TRACKNODES) THEN
          ! track nodes, this is when tracking starts for newly formed nodes
          MINNETP%NODEON = MINNETP%NODEON+1
          MINNETP%NODETRACK(MINNETP%NODEON) = GROWNODE
        ENDIF

        ! if edge we just fused to was a boundary edge, set new edge as boundary edge and fix node
        IF (MINNETP%EDGEBND(EID)) THEN
          MINNETP%EDGEBND(EIDIN) = .TRUE.
          ! no longer growing, update and end subroutine
          MINNETP%NODEGROW(GROWNODE) = 0
          ! now nodedir stores the direction of boundary edge (used for two purposes)
          MINNETP%NODEDIR(GROWNODE,:) = MINNETP%EDGEDIR(EID,:)
          MINNETP%NODEFIX(GROWNODE) = .FALSE.
          MINNETP%NODEBND(GROWNODE) = .TRUE.
        ELSE
          ! no longer growing, update and end subroutine
          MINNETP%NODEGROW(GROWNODE) = 0
          MINNETP%NODEDIR(GROWNODE,:) = 0D0
          MINNETP%NODEFIX(GROWNODE) = .FALSE.
        ENDIF

        MINNETP%EDGENODEP(EID,:) = MINNETP%EDGENODE(EID,:)
        MINNETP%EDGENODEP(EIDIN,:) = MINNETP%EDGENODE(EIDIN,:)

        FUSED=.TRUE.

        ! if the edge we just broke into two pieces was the max-length edge, find new max-length edge
        IF (MINNETP%EIDMAX.EQ.EID) THEN
          CALL FINDMAXEDGE(MINNETP)
          IF (VERBOSE) THEN 
            PRINT*, 'Growth fusion event split max edge. Old and new:', EID, MINNETP%EIDMAX
            PRINT*, 'Length:', MINNETP%EDGELEN(MINNETP%EIDMAX)
          ENDIF
        ENDIF
        RETURN

      ENDIF

      ! remove this event if didn't fuse, and move to next crossing
      NCROSS = NCROSS - 1
      TSTARS = PACK(TSTARS,CROSSES.NE.EID,TSTARS)
      TS = PACK(TS,CROSSES.NE.EID,TS)
      CROSSES = PACK(CROSSES,CROSSES.NE.EID,CROSSES)

    ENDDO

  END SUBROUTINE FUSECHECK


  SUBROUTINE GROWANDPIN(MINNETP)
    ! check for growth and pin/unpin events. All Poisson processes with a rate given by KEYS
    USE KEYS, ONLY : GROWTH, CATA, PIN, UNPIN, DELT, CENTERENCLOSE, RCIRC, GVEL, DX, GASTD
    USE mt19937, ONLY : GRND, RNORM
    USE GENUTIL, ONLY : BINSEARCH, PI
    USE STACKUTIL, ONLY : POP
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER :: NC, EC, I, CT, EID, NID, EIDIN, NIDGROW, EIDGROW
    DOUBLE PRECISION :: RN, RN2, RN3, RN4, PROB,  DIR(MINNETP%DIM)
    DOUBLE PRECISION :: EDGECDF(MINNETP%NEDGE), RN5, ANGLE
    INTEGER :: NGROW, MAXK, KC, NCATA, GROWING(MINNETP%NNODE)
    INTEGER :: NUNPIN, NPIN, PINNEDNODES(MINNETP%NNODE), UNPINNEDNODES(MINNETP%NNODE)
    DOUBLE PRECISION :: NGCDF(MINNETP%NNODE), CDELT, RGRO, LTOT, RPOS, LEN1, RGPOS, POS(MINNETP%DIM)
    DOUBLE PRECISION :: UPDELT, RPIN
    LOGICAL :: FUSED

    !--------------------
    ! catastrophe and attaching events
    !--------------------
    ! catastrophe is when a growing node suddenly stops growing, begins to retract
    ! attaching event is when a growing/retracting node fuses with an edge along its path

    NGROW=0
    ! loop over nodes to count number of growing

    DO NC=1,MINNETP%NODEHIGH

      IF (MINNETP%NODEGROW(NC).EQ.0) THEN
        CYCLE
      ELSEIF (MINNETP%NODEGROW(NC).EQ.1) THEN
        ! growing node, check for nearby edges to fuse

        CALL FUSECHECK(MINNETP,NC,FUSED)

        IF (FUSED) THEN
          ! it fused, don't add to list, continue looping through nodes
          CYCLE
        ELSE
          ! if didn't fuse, add to growing list
          NGROW=NGROW+1
          GROWING(NGROW) = NC
        ENDIF
      ELSE
        ! retracting node, check for nearby edges to fuse
        CALL FUSECHECK(MINNETP,NC,FUSED)

      ENDIF
    ENDDO
    IF (VERBOSE) THEN 
      PRINT*, 'count of growing nodes', NGROW
      IF (NGROW.GT.0) THEN
        PRINT*, 'list of growing nodes', GROWING(1:NGROW)
      ENDIF
    ENDIF


    IF (NGROW.GT.0) THEN

      CDELT = CATA*DELT*NGROW
      MAXK = MAX(INT(10*CDELT),10)
      MAXK = MIN(NGROW,MAXK)

      NGCDF = 0D0
      DO KC = 0,MAXK
        IF (KC.EQ.0) THEN
          NGCDF(1) = EXP( -CDELT )
        ELSE
          NGCDF(KC+1) = NGCDF(KC) + EXP( KC * LOG(CDELT) - CDELT - LOG_GAMMA(DBLE(KC+1)) )
        ENDIF
      ENDDO
      ! sample random number from 0 to 1, used to determine how many catastrophe events will occur
      RN = GRND()

      ! if random number is greater than the sampled CDF, it only happened because 
      ! there are a small number of growth events, and all of them catastrophizzle out
      IF (RN.GE.NGCDF(MAXK+1)) THEN
        PRINT*, 'ERROR in GROWANDPIN (catastrophe): random number greater than sampled CDF'
        PRINT*, RN
        PRINT*, NGCDF(MAXK+1)
        NCATA = NGROW
        PRINT*, 'number of catastrophe events, ', NCATA
      ELSE
        NCATA = BINSEARCH(NGCDF(1:(MAXK+1)),RN) - 1
        ! gives how many catastrophe events occur
        IF (VERBOSE) THEN
          PRINT*, 'number of catastrophe events, ', NCATA
        ENDIF
      ENDIF

      DO WHILE (NCATA.GT.0)
        NC = FLOOR(GRND()*NGROW)+1
        ! node ID from growing list
        NID = GROWING(NC)
        ! now the edge is retracting, detach motor, allow for unfixed motion (BD and length minimizing)
        MINNETP%NODEGROW(NID) = 2
        MINNETP%NODEFIX(NID) = .FALSE.
        IF (VERBOSE) THEN
          PRINT*, 'NID retracting', NID
        ENDIF

        GROWING = PACK(GROWING,GROWING.NE.NID,GROWING)

        NCATA = NCATA-1
        NGROW = NGROW-1
      ENDDO


    ENDIF


    !--------------------
    ! growth events
    !--------------------


    ! creates edgelength CDF, used when finding which edge growth event occurs on
    DO EC = 1,MINNETP%EDGEHIGH
      IF (EC.EQ.1) THEN ! first active edge
        ! if this is a boundary edge, give it zero weight
        IF (MINNETP%EDGEBND(EC)) THEN
          EDGECDF(EC) = 0D0
        ELSE
          EDGECDF(EC) = MINNETP%EDGELEN(EC)
        ENDIF
      ELSE
        ! if this is a boundary edge, give it zero weight
        IF (MINNETP%EDGEBND(EC)) THEN
          EDGECDF(EC) = EDGECDF(EC-1) + 0D0
        ELSE
          EDGECDF(EC) = EDGECDF(EC-1) + MINNETP%EDGELEN(EC) ! cumulative edgelength
        ENDIF
      ENDIF
    ENDDO

    ! normalize edgecdf by final value
    LTOT = EDGECDF(MINNETP%EDGEHIGH)
    EDGECDF = EDGECDF/LTOT


    ! if rate of growth events is very small, only need to consider 10 events happening per timestep
    ! the probability plummets for large k
    RGRO = GROWTH*LTOT*DELT
    MAXK = MAX(INT(10*RGRO),10)

    NGCDF = 0D0
    DO KC = 0,MAXK
      IF (KC.EQ.0) THEN
        NGCDF(1) = EXP( -RGRO )
      ELSE
        NGCDF(KC+1) = NGCDF(KC) + EXP( KC * LOG(RGRO) - RGRO - LOG_GAMMA(DBLE(KC+1)) )
      ENDIF
    ENDDO
    ! sample random number from 0 to 1, used to determine how many growth events will occur
    RN = GRND()

    IF (RN.GE.NGCDF(MAXK+1)) THEN
      PRINT*, 'ERROR in GROWANDPIN (growth): random number greater than sampled CDF'
      PRINT*, RN
      PRINT*, NGCDF(MAXK+1)
      STOP 1
    ENDIF

    ! gives how many growth events occur
    NGROW = BINSEARCH(NGCDF(1:(MAXK+1)),RN) - 1
    IF (VERBOSE) THEN
      PRINT*, 'number of growth events, ', NGROW
    ENDIF

    ! loop while growth events remain, randomly select edges weighted by length
    IF (NGROW.GT.0) THEN
      DO WHILE (NGROW.GT.0)

        ! growth event is set to occur, now sample from all edges to determine where it happens
        RN2 = GRND()
        ! EID gives edge from which growth is set to occur
        EID = BINSEARCH(EDGECDF(1:MINNETP%EDGEHIGH),RN2)
        IF (VERBOSE) THEN
          PRINT*, 'Growth event'
          PRINT*, 'EID', EID
        ENDIF

        ! where along edge?
        IF (EID.EQ.1) THEN
          RN3 = RN2/MINNETP%EDGELEN(1)*LTOT
        ELSE
          RN3 = (RN2 - EDGECDF(EID-1))/MINNETP%EDGELEN(EID)*LTOT
        ENDIF

        ! what direction to grow in?

        ! normally distributed around pi/2 or -pi/2
        RN4 = CEILING(GRND()*2) ! determine sign of angle (ie. which side of edge to grow out of)
        RN5 = RNORM()
        ANGLE = (-1)**(RN4) * (PI/2 + RN5*GASTD) + ATAN2(MINNETP%EDGEDIR(EID,2),MINNETP%EDGEDIR(EID,1))

        DIR(1) = COS(ANGLE)
        DIR(2) = SIN(ANGLE)
        
        IF (CENTERENCLOSE) THEN
          ! length of edge up to intersection
          LEN1 = RN3*MINNETP%EDGELEN(EID)
          ! new node radial position in space 
          POS = MINNETP%EDGESTART(EID,:) + LEN1*MINNETP%EDGEDIR(EID,:)
          RPOS = SQRT(SUM((POS)**2))
          ! new growth node radial position in space
          RGPOS = SQRT(SUM((POS + (GVEL*DELT+2*DX)*DIR)**2))

          ! if either is past the boundary circle, do not insert node, skip this growth event
          IF ((RPOS.GE.RCIRC).OR.(RGPOS.GE.RCIRC)) THEN
            NGROW = NGROW - 1
            CYCLE
          ENDIF
        ENDIF

        ! pop tops of NODESTACK and EDGESTACK for new node/edge
        CALL POP(MINNETP%NODESTACK,NID)
        CALL POP(MINNETP%EDGESTACK,EIDIN)

        ! if not centerenclose, or if rpos and initial growth node pos were not beyond boundary circle can try inserting 
        CALL INSERTINTERNODE(MINNETP,EID,RN3,NID,EIDIN,.FALSE.)
        MINNETP%NODEGROW(NID) = 0
        ! this node is not fixed
        MINNETP%NODEFIX(NID) = .FALSE.
        ! set previous position
        MINNETP%NODEPOSP(NID,:) = MINNETP%NODEPOS(NID,:)

        ! now add growth node, attached to new node above
        CALL POP(MINNETP%NODESTACK,NIDGROW)
        CALL POP(MINNETP%EDGESTACK,EIDGROW)

        ! insert a d1 node
        CALL INSERTD1NODE(MINNETP,NIDGROW,EIDGROW,DIR,NID)

        ! set new node as growing
        MINNETP%NODEGROW(NIDGROW) = 1
        ! set direction of new node growth
        MINNETP%NODEDIR(NIDGROW,:) = DIR
 
        ! recalculate EDGECDF, since new edges were added
        ! don't bother if it's the last step
        IF (NGROW.GT.1) THEN
          DO EC = 1,MINNETP%EDGEHIGH
            IF (EC.EQ.1) THEN ! first active edge
              ! if this is a boundary edge, give it zero weight
              IF (MINNETP%EDGEBND(EC)) THEN
                EDGECDF(EC) = 0D0
              ELSE
                EDGECDF(EC) = MINNETP%EDGELEN(EC)
              ENDIF
            ELSE
              ! if this is a boundary edge, give it zero weight
              IF (MINNETP%EDGEBND(EC)) THEN
                EDGECDF(EC) = EDGECDF(EC-1) + 0D0
              ELSE
                EDGECDF(EC) = EDGECDF(EC-1) + MINNETP%EDGELEN(EC) ! cumulative edgelength
              ENDIF
            ENDIF
          ENDDO

          LTOT = EDGECDF(MINNETP%EDGEHIGH)
          EDGECDF = EDGECDF/LTOT
        ENDIF

        ! finished growth event, decrease NGROW (remaining growth events)
        NGROW = NGROW - 1
      ENDDO
    ENDIF

    !-----------------------------!
    !--------- Unpinning ---------!
    !-----------------------------!

    NPIN = 0

    ! find pinned nodes

    DO NC = 1,MINNETP%NODEHIGH
      IF (.NOT.MINNETP%NODEPIN(NC)) THEN
        CYCLE
      ELSEIF (MINNETP%NODEPIN(NC)) THEN
        NPIN = NPIN + 1
        PINNEDNODES(NPIN) = NC
      ENDIF
    ENDDO

    IF (VERBOSE) THEN
      PRINT*, 'number of pinned nodes: ', NPIN
    ENDIF

    ! see if any nodes will unpin from Poisson process of unpinning
    ! default case depends on number of pinned nodes
    UPDELT = UNPIN*DELT*NPIN
    MAXK = MAX(INT(10*UPDELT),10)
    MAXK = MIN(NPIN,MAXK)

    NGCDF = 0D0
    DO KC = 0,MAXK
      IF (KC.EQ.0) THEN
        NGCDF(1) = EXP(-UPDELT)
      ELSE
        NGCDF(KC+1) = NGCDF(KC) + EXP( KC * LOG(UPDELT) - UPDELT - LOG_GAMMA(DBLE(KC+1)) )
      ENDIF
    ENDDO
    ! sample random number from 0 to 1, used to determine how many unpinning events will occur
    RN = GRND()

    ! if random number is greater than the sampled CDF, it only happened because 
    ! there are a small number of pinned nodes and they all unpin
    IF (RN.GE.NGCDF(MAXK+1)) THEN
      PRINT*, 'ERROR in GROWANDPIN (unpinning): random number greater than sampled CDF'
      PRINT*, RN
      PRINT*, NGCDF(MAXK+1)
      NUNPIN = NPIN
      PRINT*, 'number of unpinning events, ', NUNPIN
    ELSE
      NUNPIN = BINSEARCH(NGCDF(1:(MAXK+1)),RN) - 1
      ! gives how many unpinning events occur
      IF (VERBOSE) THEN
        PRINT*, 'number of unpinning events, ', NUNPIN
      ENDIF
    ENDIF

    ! now actually unpin NUNPIN randomly chosen pinned nodes
    DO WHILE (NUNPIN.GT.0)
      NC = FLOOR(GRND()*NPIN)+1
      ! node ID from pinned list
      NID = PINNEDNODES(NC)

      MINNETP%NODEPIN(NID) = .FALSE.
      
      IF (VERBOSE) THEN
        PRINT*, 'NID unpinned', NID
      ENDIF

      PINNEDNODES = PACK(PINNEDNODES,PINNEDNODES.NE.NID,PINNEDNODES)

      NUNPIN = NUNPIN-1
      NPIN = NPIN-1
    ENDDO

    IF (VERBOSE) THEN
      PRINT*, 'number of pinned nodes: ', NPIN
    ENDIF

    !-----------------------------!
    !---------- Pinning ----------!
    !-----------------------------!
 
    ! count all non-boundary, non-growing, non-fixed, non-pinned nodes
    IF (PIN.GT.0D0) THEN

    NUNPIN=0
    UNPINNEDNODES=0

    DO NC = 1,MINNETP%NODEHIGH
      IF (.NOT.(MINNETP%NODEBND(NC).OR.(MINNETP%NODEGROW(NC).GT.0).OR.MINNETP%NODEFIX(NC)&
        &.OR.MINNETP%NODEPIN(NC))) THEN
        NUNPIN=NUNPIN+1
        UNPINNEDNODES(NUNPIN)=NC
      ENDIF
    ENDDO

    IF (NUNPIN.EQ.0) THEN
      ! no nodes to unpin
    ELSE
      ! same as unpinning, constant rate per second, independent of how many nodes there are
      RPIN = PIN*DELT
      MAXK = MAX(INT(10*RPIN),10)
 
      NGCDF = 0D0
      DO KC = 0,MAXK
        IF (KC.EQ.0) THEN
          NGCDF(1) = EXP(-RPIN)
        ELSE
          NGCDF(KC+1) = NGCDF(KC) + EXP( KC * LOG(RPIN) - RPIN - LOG_GAMMA(DBLE(KC+1)) )
        ENDIF
      ENDDO
      ! sample random number from 0 to 1, used to determine how many pinning events will occur
      RN = GRND()
 
      IF (RN.GE.NGCDF(MAXK+1)) THEN
        PRINT*, 'warning in GROWANDPIN (pinning): random number greater than sampled CDF'
        PRINT*, RN
        PRINT*, NGCDF(MAXK+1)
        NPIN = NUNPIN
        PRINT*, 'number of pinning events, ', NPIN
      ELSE
        NPIN = BINSEARCH(NGCDF(1:(MAXK+1)),RN) - 1

        ! correct if somehow greater than the number of unpinned nodes
        IF (NPIN.GT.NUNPIN) THEN
          IF (VERBOSE) THEN
            PRINT*, 'number of pinning events prior to correction, ', NPIN
          ENDIF
          NPIN = NUNPIN
        ENDIF

        ! gives how many pinning events occur
        IF (VERBOSE) THEN
          PRINT*, 'number of pinning events, ', NPIN
        ENDIF
      ENDIF
 
      ! now actually pin NPIN randomly chosen unpinned nodes
      DO WHILE (NPIN.GT.0)
        NC = FLOOR(GRND()*NUNPIN)+1
        ! node ID from pinned list
        NID = UNPINNEDNODES(NC)
        MINNETP%NODEPIN(NID) = .TRUE.
        
        IF (VERBOSE) THEN
          PRINT*, 'NID pinned', NID
        ENDIF

        UNPINNEDNODES = PACK(UNPINNEDNODES,UNPINNEDNODES.NE.NID,UNPINNEDNODES)
        NPIN = NPIN-1
        NUNPIN = NUNPIN -1
      ENDDO
    ENDIF

    ENDIF

    IF (VERBOSE) THEN
      PRINT*, 'number of unpinned nodes: ', NUNPIN
    ENDIF


  END SUBROUTINE GROWANDPIN



  SUBROUTINE INSERTD1NODE(MINNETP,NIDGROW,EIDGROW,DIR,NID)
    ! insert d1 node at (GVEL*DT+2*DX)*DIR + NODEPOS(NID), using NIDGROW and EIDGROW slots for new node/edge
    USE KEYS, ONLY : DX, GVEL, DELT
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NIDGROW, EIDGROW, NID
    DOUBLE PRECISION :: DIR(MINNETP%DIM)
    INTEGER :: DEG

    ! setup new node
    ! MINNETP%NODEPOS(NIDGROW,:) = MINNETP%NODEPOS(NID,:) + 10*DX*DIR
    MINNETP%NODEPOS(NIDGROW,:) = MINNETP%NODEPOS(NID,:) + (GVEL*DELT+2*DX)*DIR
    MINNETP%NODENODE(NIDGROW,1) = NID
    MINNETP%NODEDEG(NIDGROW) = 1
    MINNETP%NODEEDGE(NIDGROW,1) = EIDGROW
    MINNETP%NODEFIX(NIDGROW) = .TRUE.
    MINNETP%NODEACT(NIDGROW) = .TRUE.
    MINNETP%NACTIVE = MINNETP%NACTIVE + 1
    IF (NIDGROW.GT.MINNETP%NODEHIGH) THEN
      MINNETP%NODEHIGH = NIDGROW
    ENDIF
    ! MINNETP%NODEHIGH = FINDLOC(MINNETP%NODEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)
    ! previous position, set it so that it is 2DX from edge
    MINNETP%NODEPOSP(NIDGROW,:) = MINNETP%NODEPOS(NID,:) + (GVEL*DELT)*DIR
    MINNETP%NODEBND(NIDGROW) = .FALSE.

    ! setup new edge, starting from NID, ending at NIDGROW
    MINNETP%EDGENODE(EIDGROW,1) = NID
    MINNETP%EDGENODE(EIDGROW,2) = NIDGROW
    MINNETP%EDGENODEP(EIDGROW,:) = MINNETP%EDGENODE(EIDGROW,:)
    MINNETP%EDGESTART(EIDGROW,:) = MINNETP%NODEPOS(NID,:)
    MINNETP%EDGEDIR(EIDGROW,:) = MINNETP%NODEPOS(NIDGROW,:) - MINNETP%NODEPOS(NID,:)
    MINNETP%EDGELEN(EIDGROW) = SQRT(SUM( MINNETP%EDGEDIR(EIDGROW,:)*MINNETP%EDGEDIR(EIDGROW,:) ))
    ! go back and normalize edgedir
    MINNETP%EDGEDIR(EIDGROW,:) = MINNETP%EDGEDIR(EIDGROW,:)/MINNETP%EDGELEN(EIDGROW)
    MINNETP%EDGEACT(EIDGROW) = .TRUE.
    MINNETP%EACTIVE = MINNETP%EACTIVE + 1
    IF (EIDGROW.GT.MINNETP%EDGEHIGH) THEN
      MINNETP%EDGEHIGH = EIDGROW
    ENDIF
    ! MINNETP%EDGEHIGH = FINDLOC(MINNETP%EDGEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)
    MINNETP%EDGEBND(EIDGROW) = .FALSE.

    ! adjust already existing node
    DEG = MINNETP%NODEDEG(NID)+1
    MINNETP%NODEDEG(NID) = DEG
    MINNETP%NODENODE(NID,DEG) = NIDGROW
    MINNETP%NODEEDGE(NID,DEG) = EIDGROW


  END SUBROUTINE INSERTD1NODE



  SUBROUTINE WIPEEDGE(MINNETP,EID)
    ! remove edge EID, make inactive and update all related fields
    USE STACKUTIL, ONLY: PUSH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: EID

    ! remove edge
    MINNETP%EDGENODE(EID,1:2) = 0
    MINNETP%EDGESTART(EID,:) = 0D0
    MINNETP%EDGEDIR(EID,:) = 0D0
    MINNETP%EDGELEN(EID) = 0D0
    MINNETP%EDGEACT(EID) = .FALSE.
    MINNETP%EACTIVE = MINNETP%EACTIVE - 1
    MINNETP%EDGEHIGH = FINDLOC(MINNETP%EDGEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)
    MINNETP%EDGEBND(EID) = .FALSE.

    CALL PUSH(MINNETP%EDGESTACK,EID)

    IF (VERBOSE) THEN
      PRINT*, 'Deactivated edge, ', EID
    ENDIF


  END SUBROUTINE WIPEEDGE



  SUBROUTINE WIPENODE(MINNETP,NID)
    ! remove node NID, make inactive and update all related fields
    USE STACKUTIL, ONLY: PUSH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NID

    ! remove node
    MINNETP%NODEPOS(NID,:) = 0D0
    MINNETP%NODENODE(NID,1) = 0
    MINNETP%NODEDEG(NID) = 0
    MINNETP%NODEEDGE(NID,1) = 0
    MINNETP%NODEFIX(NID) = .FALSE.
    MINNETP%NODEPIN(NID) = .FALSE.
    MINNETP%NODEACT(NID) = .FALSE.
    MINNETP%NACTIVE = MINNETP%NACTIVE - 1
    MINNETP%NODEGROW(NID) = 0
    MINNETP%NODEDIR(NID,:) = 0D0
    MINNETP%NODEHIGH = FINDLOC(MINNETP%NODEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)
    MINNETP%NODEBND(NID) = .FALSE.
    
    CALL PUSH(MINNETP%NODESTACK,NID)

    IF (VERBOSE) THEN
      PRINT*, 'Deactivated node, ', NID
    ENDIF


  END SUBROUTINE WIPENODE
    


  SUBROUTINE REMOVED1NODE(MINNETP,NID)
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NID
    INTEGER :: IND, I, EID, NEIGH
    LOGICAL :: REMOVED

    EID = MINNETP%NODEEDGE(NID,1)
    NEIGH = MINNETP%NODENODE(NID,1)

    ! ! track node, probably don't need to worry about node tracking here
    ! ! since only ever tracking d3/d2 nodes now!
    ! IF (.NOT.MINNETP%NODEDONE(NID)) THEN
    !   MINNETP%NODETRACK(NID) = NEIGH
    !   MINNETP%NODEDONE(NID) = .TRUE.
    ! ENDIF

    CALL WIPENODE(MINNETP,NID)
    CALL WIPEEDGE(MINNETP,EID)

    ! adjust already existing node
    IND = FINDLOC(MINNETP%NODEEDGE(NEIGH,:),DIM=1,VALUE=EID)
    ! removing NID from NEIGH connections
    REMOVED=.FALSE.
    DO I = 1,MINNETP%NODEDEG(NEIGH)
      IF (IND.EQ.I) THEN
        REMOVED=.TRUE.
        CYCLE
      ENDIF
      IF (REMOVED) THEN
        MINNETP%NODENODE(NEIGH,I-1) = MINNETP%NODENODE(NEIGH,I)
        MINNETP%NODEEDGE(NEIGH,I-1) = MINNETP%NODEEDGE(NEIGH,I)
      ENDIF
    ENDDO
    ! after do loop, want to remove unused nodenode and nodeedge values, 
    ! ie. entries past the nodedeg. just to keep these arrays clean
    MINNETP%NODENODE(NEIGH,MINNETP%NODEDEG(NEIGH)) = 0
    MINNETP%NODEEDGE(NEIGH,MINNETP%NODEDEG(NEIGH)) = 0
    ! update node degree, lost a connection
    MINNETP%NODEDEG(NEIGH) = MINNETP%NODEDEG(NEIGH) - 1

  END SUBROUTINE REMOVED1NODE



  SUBROUTINE T1(MINNETP,NODE1,NODE2)
    ! T1 rearrangemnt
    USE KEYS, ONLY : MAXBRANCH, DX
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE1, NODE2
    INTEGER :: NODES1(MAXBRANCH), NODES2(MAXBRANCH)
    INTEGER :: TMPNODE, NC, EC, CT, ECC, TMPNODES2(MAXBRANCH), TMPNODES1(MAXBRANCH)
    INTEGER :: NODE11, NODE12, NODE21, NODE22, INVEDGES(2*MAXBRANCH-1)
    DOUBLE PRECISION :: TOTLEN, NEWLEN1, NEWLEN2, NEWLEN3, NEWLEN4, NEWLEN5
    DOUBLE PRECISION :: DIST, TMPNODEPOS(MINNETP%DIM)
    INTEGER :: NN1(3), NN2(3), CT1, CT2, EIDS(2)
    INTEGER :: EDGE21, EDGE11
    LOGICAL :: N1BOOL, N2BOOL


    ! first check if this involves a degenerate loop that has shrunk to below DX
    INVEDGES = 0
    EIDS = 0
    NN1 = MINNETP%NODENODE(NODE1,1:3)
    NN2 = MINNETP%NODENODE(NODE2,1:3)
    NODES1=0
    NODES2=0

    CT = 0
    ! these keep track of how many times each node appears in the other's NODENODE
    CT1 = 0
    CT2 = 0
    ! find their neighbors, and keep track of "involved edges" INVEDGES
    DO NC = 1,3
      CT = CT + 1
      TMPNODE = NN1(NC)
      IF (TMPNODE.EQ.NODE2) THEN
        CT1 = CT1+1
        ! track which edges are between node 1 and 2 (especially if there are two!)
        EIDS(CT1) = MINNETP%NODEEDGE(NODE1,NC)
        CYCLE
      ENDIF
      CT = CT - 1
      INVEDGES(NC-CT) = MINNETP%NODEEDGE(NODE1,NC)
      NODES1(NC-CT) = NC
    ENDDO
    CT = 0
    DO NC = 1,3
      CT = CT + 1
      TMPNODE = NN2(NC)
      IF (TMPNODE.EQ.NODE1) THEN
        CT2 = CT2+1
        CYCLE
      ENDIF
      CT = CT - 1
      NODES2(NC-CT) = NC
      INVEDGES(3 - CT1 + NC-CT) = MINNETP%NODEEDGE(NODE2,NC)
    ENDDO

    ! check if ct1=ct2
    IF (CT1.NE.CT2) THEN
      PRINT*, 'ERROR IN T1: discrepancy in nodenode for nodes ', NODE1, NODE2
      PRINT*, NN1
      PRINT*, NN2
      STOP 1
    ENDIF
    
    ! if there is a degenerate edge, and it is the short edge, remove the loop and done with this call to T1
    IF (CT1.EQ.2) THEN
      IF (VERBOSE) THEN
        PRINT*, INVEDGES
        PRINT*, EIDS
        PRINT*, NODES1,NODES2
      ENDIF
      CALL REMOVELOOP(MINNETP,NODE1,NODE2,NODES1(1),NODES2(1),EIDS,INVEDGES)
      RETURN
    ENDIF

    ! ---------------------
    ! T1 involving boundary
    ! ---------------------
    IF (ANY(MINNETP%EDGEBND(INVEDGES(1:4)))) THEN

      ! check two possible rearrangements: NODE1 becomes doubly connected off of boundary
      ! or NODE2 is doubly connected off of boundary, other node taking boundary edges
      ! the nodes stay where they are, and edge connecting them remains

      ! check if node 2 is connected to boundary edge
      IF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE2,NODES2(1)))) THEN
        EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(1))
        TMPNODES2 = NODES2
      ELSEIF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE2,NODES2(2)))) THEN
        EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(2))
        TMPNODES2(2) = NODES2(1)
        TMPNODES2(1) = NODES2(2)
      ELSE
        ! attempting T1 between a node on the boundary and a node not on the boundary
        IF (VERBOSE) THEN
          PRINT*, 'T1 between bnd and not bnd'
        ENDIF
        CALL BNDT1(MINNETP,NODE2,NODE1,NODES2,NODES1,INVEDGES)

        ! whether or not a T1 was carried out, return
        RETURN
      ENDIF

      IF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE1,NODES1(1)))) THEN
        EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
        TMPNODES1 = NODES1
      ELSEIF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE1,NODES1(2)))) THEN
        EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(2))
        TMPNODES1(2) = NODES1(1)
        TMPNODES1(1) = NODES1(2)
      ELSE
        ! attempting T1 between a node on the boundary and a node not on the boundary
        IF (VERBOSE) THEN
          PRINT*, 'T1 between bnd and not bnd'
        ENDIF
        CALL BNDT1(MINNETP,NODE1,NODE2,NODES1,NODES2,INVEDGES)

        ! whether or not a T1 was carried out, return
        RETURN
      ENDIF

      IF (VERBOSE) THEN
        PRINT*, 'T1 between two boundary nodes'
      ENDIF

      ! now we have boundary edges stored in TMPNODES1(1) and TMPNODES2(1)
      ! could be the case that the other edges are also boundary edges...
      ! next check possible rearrangements

      TOTLEN = 0D0
      DO ECC = 1,4
        TOTLEN = TOTLEN + MINNETP%EDGELEN(INVEDGES(ECC))
      ENDDO

      ! since there is a chance of projecting nodes here, include edge between nodes 1 and 2
      CALL NODENODEDIST(MINNETP,NODE1,NODE2,DIST)
      TOTLEN = TOTLEN + DIST

      NODE11 = MINNETP%NODENODE(NODE1,TMPNODES1(1))
      NODE12 = MINNETP%NODENODE(NODE1,TMPNODES1(2))
      NODE21 = MINNETP%NODENODE(NODE2,TMPNODES2(1))
      NODE22 = MINNETP%NODENODE(NODE2,TMPNODES2(2))

      ! first combo:
      ! 11 <-> 22
      IF (MINNETP%NODEFIX(NODE1).OR.MINNETP%NODEPIN(NODE1)) THEN
        NEWLEN1 = TOTLEN+1D3 ! this transition is not allowed if node 1 is a corner/fixed node
      ELSE
        NEWLEN1 = SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE22,:))**2))
        NEWLEN1 = NEWLEN1 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE12,:))**2))
        NEWLEN1 = NEWLEN1 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE21,:))**2))
        NEWLEN1 = NEWLEN1 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE11,:))**2)) + DIST
      ENDIF


      ! second combo:
      ! 11 <-> 21
      IF (MINNETP%NODEFIX(NODE1).OR.MINNETP%NODEPIN(NODE1)) THEN
        ! if node 1 is fixed, then project node 2 position onto line connecting nodes 1 and 11
        ! so node slides along boundary
        CALL GETPROJ(MINNETP,NODE2,EDGE11,TMPNODEPOS)

        IF (VERBOSE) THEN
          PRINT*, 'projected node 2 onto 11. Old and new pos:, ', MINNETP%NODEPOS(NODE2,:), TMPNODEPOS
        ENDIF

        NEWLEN2 = SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE21,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE12,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((TMPNODEPOS - MINNETP%NODEPOS(NODE11,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((TMPNODEPOS - MINNETP%NODEPOS(NODE22,:))**2))
        ! distance between node 1 and 2 
        NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - TMPNODEPOS)**2))
      ELSEIF (MINNETP%NODEFIX(NODE2).OR.MINNETP%NODEPIN(NODE2)) THEN
        ! if node 2 is fixed, then project node 1 position onto line connecting nodes 2 and 21
        ! so node slides along boundary
        CALL GETPROJ(MINNETP,NODE1,EDGE21,TMPNODEPOS)

        IF (VERBOSE) THEN
          PRINT*, 'projected node 1 onto 21. Old and new pos:, ', MINNETP%NODEPOS(NODE1,:), TMPNODEPOS
        ENDIF

        NEWLEN2 = SQRT(SUM((TMPNODEPOS - MINNETP%NODEPOS(NODE21,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((TMPNODEPOS - MINNETP%NODEPOS(NODE12,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE11,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE22,:))**2))
        ! distance between node 1 and 2 
        NEWLEN2 = NEWLEN2 + SQRT(SUM((TMPNODEPOS - MINNETP%NODEPOS(NODE2,:))**2))

      ELSE
        ! no corner nodes, carry on as usual
        NEWLEN2 = SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE21,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE12,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE11,:))**2))
        NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE22,:))**2)) + DIST
        ! also adding on DIST

      ENDIF


      ! third combo:
      ! 12 <-> 21
      IF (MINNETP%NODEFIX(NODE2).OR.MINNETP%NODEPIN(NODE2)) THEN
        NEWLEN3 = TOTLEN+1D3 ! this transition is not allowed if node 2 is a corner/fixed node
      ELSE
        NEWLEN3 = SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE11,:))**2))
        NEWLEN3 = NEWLEN3 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE21,:))**2))
        NEWLEN3 = NEWLEN3 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE12,:))**2))
        NEWLEN3 = NEWLEN3 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE22,:))**2)) + DIST
        ! also adding on DIST
      ENDIF

      ! fourth combo:
      ! 11 <-> 21, 12 <-> 22
      IF ((MINNETP%NODEFIX(NODE1)).OR.(MINNETP%NODEFIX(NODE2)) &
        & .OR.(MINNETP%NODEPIN(NODE1)).OR.(MINNETP%NODEPIN(NODE2))) THEN
        NEWLEN4 = TOTLEN+1D3 ! this transition is not allowed if either node is a corner/fixed node
      ELSE
        NEWLEN4 = SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE21,:))**2))
        NEWLEN4 = NEWLEN4 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE22,:))**2))
        NEWLEN4 = NEWLEN4 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE11,:))**2))
        NEWLEN4 = NEWLEN4 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE12,:))**2)) + DIST
        ! also adding on DIST
      ENDIF

      ! fifth combo:
      ! 12 <-> 22
      NEWLEN5 = SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE11,:))**2))
      NEWLEN5 = NEWLEN5 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE22,:))**2))
      NEWLEN5 = NEWLEN5 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE21,:))**2))
      NEWLEN5 = NEWLEN5 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE12,:))**2)) + DIST
      ! also adding on DIST

      IF (VERBOSE) THEN
        PRINT*, 'total length before rearrangement', TOTLEN
        PRINT*, 'new length after rearrangement 1', NEWLEN1
        PRINT*, 'new length after rearrangement 2', NEWLEN2
        PRINT*, 'new length after rearrangement 3', NEWLEN3
        PRINT*, 'new length after rearrangement 4', NEWLEN4
        PRINT*, 'new length after rearrangement 5', NEWLEN5
      ENDIF


      IF ((NEWLEN1.LT.TOTLEN).AND.(NEWLEN1.LE.NEWLEN2).AND.(NEWLEN1.LE.NEWLEN3)&
        & .AND.(NEWLEN1.LE.NEWLEN4).AND.(NEWLEN1.LE.NEWLEN5)) THEN
        ! combo 1
        TMPNODE = TMPNODES2(1)
        TMPNODES2(1) = TMPNODES2(2)
        TMPNODES2(2) = TMPNODE

        CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,TMPNODES1,TMPNODES2)

        ! node 1 may no longer be a boundary node
        ! check if it is still connected to a boundary node (in addition to node 2)
        N1BOOL = MINNETP%NODEBND(MINNETP%NODENODE(NODE1,TMPNODES1(1)))
        N2BOOL = MINNETP%NODEBND(MINNETP%NODENODE(NODE1,TMPNODES1(2)))
        ! if neither of its other connections are boundary nodes, no longer on the boundary!
        IF (.NOT.(N1BOOL.OR.N2BOOL)) THEN
          MINNETP%NODEBND(NODE1) = .FALSE.
          MINNETP%NODEDIR(NODE1,:) = 0D0
          ! edge connecting the two nodes is no longer boundary edge
          MINNETP%EDGEBND(EIDS(1)) = .FALSE.
        ENDIF

      ELSEIF ((NEWLEN2.LT.TOTLEN).AND.(NEWLEN2.LE.NEWLEN3).AND.(NEWLEN2.LE.NEWLEN4).AND.(NEWLEN2.LE.NEWLEN5)) THEN
        ! combo 2

        CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,TMPNODES1,TMPNODES2)
        ! node 1 and 2 retain their boundary status, since they just swapped boundary edges

      ELSEIF ((NEWLEN3.LT.TOTLEN).AND.(NEWLEN3.LE.NEWLEN4).AND.(NEWLEN3.LE.NEWLEN5)) THEN
        ! combo 3

        TMPNODE = TMPNODES1(1)
        TMPNODES1(1) = TMPNODES1(2)
        TMPNODES1(2) = TMPNODE

        CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,TMPNODES1,TMPNODES2)

        ! node 2 may no longer be a boundary node
        ! check if it is still connected to a boundary node (in addition to node 2)
        N1BOOL = MINNETP%NODEBND(MINNETP%NODENODE(NODE2,TMPNODES2(1)))
        N2BOOL = MINNETP%NODEBND(MINNETP%NODENODE(NODE2,TMPNODES2(2)))
        ! if neither of its other connections are boundary nodes, no longer on the boundary!
        IF (.NOT.(N1BOOL.OR.N2BOOL)) THEN
          MINNETP%NODEBND(NODE2) = .FALSE.
          MINNETP%NODEDIR(NODE2,:) = 0D0
          ! edge connecting the two nodes is no longer boundary edge
          MINNETP%EDGEBND(EIDS(1)) = .FALSE.
        ENDIF

      ELSEIF ((NEWLEN4.LT.TOTLEN).AND.(NEWLEN4.LE.NEWLEN5)) THEN
        ! this is not allowed
      ELSEIF (NEWLEN5.LT.TOTLEN) THEN

        TMPNODE = TMPNODES1(1)
        TMPNODES1(1) = TMPNODES1(2)
        TMPNODES1(2) = TMPNODE
        TMPNODE = TMPNODES2(1)
        TMPNODES2(1) = TMPNODES2(2)
        TMPNODES2(2) = TMPNODE

        CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,TMPNODES1,TMPNODES2)

      ENDIF

      RETURN
    ENDIF
    ! ---------------------
    ! end boundary code
    ! ---------------------

    TOTLEN = 0D0
    DO ECC = 1,4
      TOTLEN = TOTLEN + MINNETP%EDGELEN(INVEDGES(ECC))
    ENDDO
    
    ! now see if T1 rearrangement lowers total length
    NEWLEN1 = 0D0
    NEWLEN2 = 0D0
    NEWLEN3 = 0D0
    NEWLEN4 = 0D0
    NEWLEN5 = 0D0
    NODE11 = MINNETP%NODENODE(NODE1,NODES1(1))
    NODE12 = MINNETP%NODENODE(NODE1,NODES1(2))
    NODE21 = MINNETP%NODENODE(NODE2,NODES2(1))
    NODE22 = MINNETP%NODENODE(NODE2,NODES2(2))


    ! add up lengths, four possible combinations:

    ! combo 1, 11 <-> 22
    CALL NODENODEDIST(MINNETP,NODE1,NODE22,DIST)
    NEWLEN1 = NEWLEN1 + DIST
    CALL NODENODEDIST(MINNETP,NODE1,NODE12,DIST)
    NEWLEN1 = NEWLEN1 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    NEWLEN1 = NEWLEN1 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE21,DIST)
    NEWLEN1 = NEWLEN1 + DIST
 
    ! combo 2, 11 <-> 21
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    NEWLEN2 = NEWLEN2 + DIST
    CALL NODENODEDIST(MINNETP,NODE1,NODE12,DIST)
    NEWLEN2 = NEWLEN2 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    NEWLEN2 = NEWLEN2 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE22,DIST)
    NEWLEN2 = NEWLEN2 + DIST
 
    ! combo 3, 12 <-> 21
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    NEWLEN3 = NEWLEN3 + DIST
    CALL NODENODEDIST(MINNETP,NODE1,NODE11,DIST)
    NEWLEN3 = NEWLEN3 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE12,DIST)
    NEWLEN3 = NEWLEN3 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE22,DIST)
    NEWLEN3 = NEWLEN3 + DIST
 
    ! combo 4, 11 <-> 21, 12 <-> 22
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    NEWLEN4 = NEWLEN4 + DIST
    CALL NODENODEDIST(MINNETP,NODE1,NODE22,DIST)
    NEWLEN4 = NEWLEN4 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    NEWLEN4 = NEWLEN4 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE12,DIST)
    NEWLEN4 = NEWLEN4 + DIST

    ! combo 5, 12 <-> 22
    CALL NODENODEDIST(MINNETP,NODE1,NODE11,DIST)
    NEWLEN5 = NEWLEN5 + DIST
    CALL NODENODEDIST(MINNETP,NODE1,NODE22,DIST)
    NEWLEN5 = NEWLEN5 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE21,DIST)
    NEWLEN5 = NEWLEN5 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,NODE12,DIST)
    NEWLEN5 = NEWLEN5 + DIST

    IF (VERBOSE) THEN
      PRINT*, 'total length before rearrangement', TOTLEN
      PRINT*, 'new length after rearrangement 1', NEWLEN1
      PRINT*, 'new length after rearrangement 2', NEWLEN2
      PRINT*, 'new length after rearrangement 3', NEWLEN3
      PRINT*, 'new length after rearrangement 4', NEWLEN4
      PRINT*, 'new length after rearrangement 5', NEWLEN5
    ENDIF

    IF ((NEWLEN1.LT.TOTLEN).AND.(NEWLEN1.LE.NEWLEN2).AND.(NEWLEN1.LE.NEWLEN3)&
      & .AND.(NEWLEN1.LE.NEWLEN4).AND.(NEWLEN1.LE.NEWLEN5)) THEN
      ! combo 1

      TMPNODES2(1) = NODES2(2)
      TMPNODES2(2) = NODES2(1)

      CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,NODES1,TMPNODES2)

    ELSEIF ((NEWLEN2.LT.TOTLEN).AND.(NEWLEN2.LE.NEWLEN3).AND.(NEWLEN2.LE.NEWLEN4).AND.(NEWLEN2.LE.NEWLEN5)) THEN
      ! combo 2

      CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,NODES1,NODES2)

    ELSEIF ((NEWLEN3.LT.TOTLEN).AND.(NEWLEN3.LE.NEWLEN4).AND.(NEWLEN3.LE.NEWLEN5)) THEN
      ! combo 3

      TMPNODES1(1) = NODES1(2)
      TMPNODES1(2) = NODES1(1)

      CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,TMPNODES1,NODES2)

    ELSEIF ((NEWLEN4.LT.TOTLEN).AND.(NEWLEN4.LE.NEWLEN5)) THEN
      ! combo 4
 
      ! this is the special case of both neighbors being swapped
      ! not possible to form a loop in this case (except for in the intermediate step
      ! but these two calls to swapneigh happen instantaneously
      TMPNODES1(1) = NODES1(2)
      TMPNODES1(2) = NODES2(1)
      TMPNODES2(1) = NODES2(2)
      TMPNODES2(2) = NODES1(1)
 
      CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,NODES1,NODES2)
      CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,TMPNODES1,TMPNODES2)

    ELSEIF (NEWLEN5.LT.TOTLEN) THEN

      TMPNODES1(1) = NODES1(2)
      TMPNODES1(2) = NODES2(1)
      TMPNODES2(1) = NODES2(2)
      TMPNODES2(2) = NODES1(1)

      CALL SWAPNEIGHT1(MINNETP,NODE1,NODE2,TMPNODES1,TMPNODES2)

    ENDIF


  END SUBROUTINE T1



  SUBROUTINE BNDT1(MINNETP,NODE1,NODE2,NODES1,NODES2,INVEDGES)
    ! T1 rearrangemnt between a boundary node and a non-boundary node
    ! NODE1 is non-boundary node, NODE2 boundary node (either fixed or sliding)
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE1, NODE2
    INTEGER, INTENT(IN) :: NODES1(MAXBRANCH), NODES2(MAXBRANCH), INVEDGES(2*MAXBRANCH-1)
    INTEGER :: EC, TMPNODES2(MAXBRANCH), TMPNODES1(MAXBRANCH)
    INTEGER :: NODE11, NODE12, NODE21, NODE22
    DOUBLE PRECISION :: TOTLEN, TMPNODEPOS1(MINNETP%DIM), TMPNODEPOS2(MINNETP%DIM)
    DOUBLE PRECISION :: NEWLEN1, NEWLEN2, NEWLEN3, NEWLEN4, NEWLEN5
    DOUBLE PRECISION :: DIST
    ! INTEGER :: NN1(3), NN2(3), CT1, CT2, EIDS(2)
    INTEGER :: EDGE21, EDGE22

    EDGE21=MINNETP%NODEEDGE(NODE2,NODES2(1))
    EDGE22=MINNETP%NODEEDGE(NODE2,NODES2(2))

    TOTLEN = 0D0
    DO EC = 1,4
      TOTLEN = TOTLEN + MINNETP%EDGELEN(INVEDGES(EC))
    ENDDO
    ! since there is a chance of projecting nodes here, include edge between nodes 1 and 2
    CALL NODENODEDIST(MINNETP,NODE1,NODE2,DIST)
    TOTLEN = TOTLEN + DIST

    NODE11 = MINNETP%NODENODE(NODE1,NODES1(1))
    NODE12 = MINNETP%NODENODE(NODE1,NODES1(2))
    NODE21 = MINNETP%NODENODE(NODE2,NODES2(1))
    NODE22 = MINNETP%NODENODE(NODE2,NODES2(2))

    ! first combo:
    ! 11 <-> 22, project node1 onto 2-22 line
    CALL GETPROJ(MINNETP,NODE1,EDGE22,TMPNODEPOS2)

    NEWLEN1 = SQRT(SUM((TMPNODEPOS2 - MINNETP%NODEPOS(NODE22,:))**2))
    NEWLEN1 = NEWLEN1 + SQRT(SUM((TMPNODEPOS2 - MINNETP%NODEPOS(NODE12,:))**2))
    NEWLEN1 = NEWLEN1 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE21,:))**2))
    NEWLEN1 = NEWLEN1 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE11,:))**2))
    NEWLEN1 = NEWLEN1 + SQRT(SUM((TMPNODEPOS2 - MINNETP%NODEPOS(NODE2,:))**2))

    ! second combo:
    ! 11 <-> 21, project node1 onto 2-21 line
    CALL GETPROJ(MINNETP,NODE1,EDGE21,TMPNODEPOS1)

    NEWLEN2 = SQRT(SUM((TMPNODEPOS1 - MINNETP%NODEPOS(NODE21,:))**2))
    NEWLEN2 = NEWLEN2 + SQRT(SUM((TMPNODEPOS1 - MINNETP%NODEPOS(NODE12,:))**2))
    NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE11,:))**2))
    NEWLEN2 = NEWLEN2 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE22,:))**2))
    NEWLEN2 = NEWLEN2 + SQRT(SUM((TMPNODEPOS1 - MINNETP%NODEPOS(NODE2,:))**2))

    ! third combo:
    ! 12 <-> 21, project node1 onto 2-21 line
    NEWLEN3 = SQRT(SUM((TMPNODEPOS1 - MINNETP%NODEPOS(NODE11,:))**2))
    NEWLEN3 = NEWLEN3 + SQRT(SUM((TMPNODEPOS1 - MINNETP%NODEPOS(NODE21,:))**2))
    NEWLEN3 = NEWLEN3 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE12,:))**2))
    NEWLEN3 = NEWLEN3 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE22,:))**2))
    NEWLEN3 = NEWLEN3 + SQRT(SUM((TMPNODEPOS1 - MINNETP%NODEPOS(NODE2,:))**2))

    ! fourth combo:
    ! 11 <-> 21, 12 <-> 22, not allowed, node1 would be projected onto node2

    ! fifth combo:
    ! 12 <-> 22, project node1 onto 2-22 line
    NEWLEN5 = SQRT(SUM((TMPNODEPOS2 - MINNETP%NODEPOS(NODE11,:))**2))
    NEWLEN5 = NEWLEN5 + SQRT(SUM((TMPNODEPOS2 - MINNETP%NODEPOS(NODE22,:))**2))
    NEWLEN5 = NEWLEN5 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE21,:))**2))
    NEWLEN5 = NEWLEN5 + SQRT(SUM((MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE12,:))**2))
    NEWLEN5 = NEWLEN5 + SQRT(SUM((TMPNODEPOS2 - MINNETP%NODEPOS(NODE2,:))**2))

    IF (VERBOSE) THEN
      PRINT*, 'total length before rearrangement', TOTLEN
      PRINT*, 'new length after rearrangement 1', NEWLEN1
      PRINT*, 'new length after rearrangement 2', NEWLEN2
      PRINT*, 'new length after rearrangement 3', NEWLEN3
      PRINT*, 'new length after rearrangement 5', NEWLEN5
    ENDIF

    IF ((NEWLEN1.LT.TOTLEN).AND.(NEWLEN1.LE.NEWLEN2).AND.(NEWLEN1.LE.NEWLEN3)&
      & .AND.(NEWLEN1.LE.NEWLEN5)) THEN
      ! combo 1

      TMPNODES2(1) = NODES2(2)
      TMPNODES2(2) = NODES2(1)

      CALL SWAPNEIGHT1BND(MINNETP,NODE1,NODE2,NODES1,TMPNODES2,TMPNODEPOS2)


    ELSEIF ((NEWLEN2.LT.TOTLEN).AND.(NEWLEN2.LE.NEWLEN3).AND.(NEWLEN2.LE.NEWLEN5)) THEN
      ! combo 2

      CALL SWAPNEIGHT1BND(MINNETP,NODE1,NODE2,NODES1,NODES2,TMPNODEPOS1)

    ELSEIF ((NEWLEN3.LT.TOTLEN).AND.(NEWLEN3.LE.NEWLEN5)) THEN
      ! combo 3

      TMPNODES1(1) = NODES1(2)
      TMPNODES1(2) = NODES1(1)

      CALL SWAPNEIGHT1BND(MINNETP,NODE1,NODE2,TMPNODES1,NODES2,TMPNODEPOS1)

    ELSEIF (NEWLEN5.LT.TOTLEN) THEN

      TMPNODES1(1) = NODES1(2)
      TMPNODES1(2) = NODES2(1)
      TMPNODES2(1) = NODES2(2)
      TMPNODES2(2) = NODES1(1)

      CALL SWAPNEIGHT1BND(MINNETP,NODE1,NODE2,TMPNODES1,TMPNODES2,TMPNODEPOS2)

    ENDIF

  END SUBROUTINE BNDT1



  SUBROUTINE  REMOVELOOP23(MINNETP,NODE1,NODE2,NODES2,EIDS,INVEDGES)
    ! remove the two input edges, shrinking down to a node
    ! this is for the case of a degree 2 and a degree 3 node, NODE2 is deg3, save that one
    USE KEYS, ONLY : MAXBRANCH, TRACKNODES
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE1, NODE2, NODES2, EIDS(2), INVEDGES(2*MAXBRANCH-1)
    INTEGER :: NC, NN2(MAXBRANCH), NE2(MAXBRANCH), IND(1), EN11, EN12, EN21, EN22

    ! save nn2
    NN2 = MINNETP%NODENODE(NODE2,:)
    NE2 = MINNETP%NODEEDGE(NODE2,:)

    ! set new nodenode for node2
    MINNETP%NODENODE(NODE2,:) = 0
    MINNETP%NODENODE(NODE2,1) = NN2(NODES2)
    ! node degree
    MINNETP%NODEDEG(NODE2) = 1
    ! set new nodeedge for node2
    MINNETP%NODEEDGE(NODE2,:) = 0
    MINNETP%NODEEDGE(NODE2,1) = INVEDGES(1)
    ! only average positions if nodes weren't fixed
    IF (MINNETP%NODEFIX(NODE1)) THEN
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE1,:)
      MINNETP%NODEFIX(NODE2) = .TRUE.
    ELSEIF (MINNETP%NODEFIX(NODE2)) THEN
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE2,:)
    ELSEIF (MINNETP%NODEPIN(NODE1)) THEN ! if node1 was pinned
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE1,:)
      MINNETP%NODEPIN(NODE2) = .TRUE.
    ELSEIF (MINNETP%NODEPIN(NODE2)) THEN ! if node2 was pinned
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE2,:)
    ELSE
      ! new position is average of previous positions
      MINNETP%NODEPOS(NODE2,:) = (MINNETP%NODEPOS(NODE1,:) + MINNETP%NODEPOS(NODE2,:))/2
    ENDIF

    ! update edges
    ! invedge(1) has no changes to edgenode, just change the start, dir, and length
    EN11 = MINNETP%EDGENODE(INVEDGES(1),1)
    EN12 = MINNETP%EDGENODE(INVEDGES(1),2)

    IF (EN11.EQ.NODE2) THEN
      MINNETP%EDGESTART(INVEDGES(1),:) = MINNETP%NODEPOS(NODE2,:)
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%NODEPOS(EN12,:)-MINNETP%NODEPOS(NODE2,:)
      MINNETP%EDGELEN(INVEDGES(1)) = SQRT(SUM( MINNETP%EDGEDIR(INVEDGES(1),:)**2 ))
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%EDGEDIR(INVEDGES(1),:)/MINNETP%EDGELEN(INVEDGES(1))
    ELSE
      MINNETP%EDGESTART(INVEDGES(1),:) = MINNETP%NODEPOS(EN11,:)
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%NODEPOS(NODE2,:)-MINNETP%NODEPOS(EN11,:)
      MINNETP%EDGELEN(INVEDGES(1)) = SQRT(SUM( MINNETP%EDGEDIR(INVEDGES(1),:)**2 ))
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%EDGEDIR(INVEDGES(1),:)/MINNETP%EDGELEN(INVEDGES(1))
    ENDIF

    ! check if new max edge length
    IF (MINNETP%EDGELEN(INVEDGES(1)).GT.MINNETP%EDGELEN(MINNETP%EIDMAX)) THEN
      MINNETP%EIDMAX = INVEDGES(1)
    ENDIF

    ! track node
    IF (TRACKNODES) THEN
      ! find all occurrences of NODE1 and update to NODE2
      DO NC=1,MINNETP%NODEON
        IF (MINNETP%NODETRACK(NC).EQ.NODE1) THEN
          MINNETP%NODETRACK(NC) = NODE2
        ENDIF
      ENDDO
    ENDIF

    ! wipe NODE1
    CALL WIPENODE(MINNETP,NODE1)
    ! wipe two short edges
    CALL WIPEEDGE(MINNETP,EIDS(1))
    CALL WIPEEDGE(MINNETP,EIDS(2))

  END SUBROUTINE REMOVELOOP23



  SUBROUTINE  REMOVELOOP(MINNETP,NODE1,NODE2,NODES1,NODES2,EIDS,INVEDGES)
    ! remove the two input edges, shrinking down to a node
    USE KEYS, ONLY : MAXBRANCH, TRACKNODES
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE1, NODE2, NODES1, NODES2, EIDS(2), INVEDGES(2*MAXBRANCH-1)
    INTEGER :: NC, NN1(MAXBRANCH), NE1(MAXBRANCH), IND(1), EN11, EN12, EN21, EN22

    ! save nn1
    NN1 = MINNETP%NODENODE(NODE1,:)
    NE1 = MINNETP%NODEEDGE(NODE1,:)

    ! set new nodenode for node1
    MINNETP%NODENODE(NODE1,:) = 0
    MINNETP%NODENODE(NODE1,1) = NN1(NODES1)
    MINNETP%NODENODE(NODE1,2) = MINNETP%NODENODE(NODE2,NODES2)
    ! node degree
    MINNETP%NODEDEG(NODE1) = 2
    ! set new nodeedge for node1
    MINNETP%NODEEDGE(NODE1,:) = 0
    MINNETP%NODEEDGE(NODE1,1) = INVEDGES(1)
    MINNETP%NODEEDGE(NODE1,2) = INVEDGES(2)
    ! only average positions if nodes weren't fixed
    IF (MINNETP%NODEFIX(NODE1).OR.(MINNETP%NODEPIN(NODE1))) THEN
      MINNETP%NODEPOS(NODE1,:) = MINNETP%NODEPOS(NODE1,:)
    ELSEIF (MINNETP%NODEFIX(NODE2)) THEN
      MINNETP%NODEPOS(NODE1,:) = MINNETP%NODEPOS(NODE2,:)
      MINNETP%NODEFIX(NODE1) = .TRUE.
    ELSEIF (MINNETP%NODEPIN(NODE2)) THEN
      MINNETP%NODEPOS(NODE1,:) = MINNETP%NODEPOS(NODE2,:)
      MINNETP%NODEPIN(NODE1) = .TRUE.
    ELSE
      ! new position is average of previous positions
      MINNETP%NODEPOS(NODE1,:) = (MINNETP%NODEPOS(NODE1,:) + MINNETP%NODEPOS(NODE2,:))/2
    ENDIF

    ! update nodenode of neighbors
    ! this gives index of NODE2 in it's neighbors nodenode
    IND = FINDLOC(MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2),:),NODE2,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2),IND) = NODE1

    ! update edges
    ! invedge(1) has no changes to edgenode, just change the start, dir, and length
    ! invedge(2), must change edgenode in addition to spatial properties
    EN11 = MINNETP%EDGENODE(INVEDGES(1),1)
    EN12 = MINNETP%EDGENODE(INVEDGES(1),2)

    IF (EN11.EQ.NODE1) THEN
      MINNETP%EDGESTART(INVEDGES(1),:) = MINNETP%NODEPOS(NODE1,:)
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%NODEPOS(EN12,:)-MINNETP%NODEPOS(NODE1,:)
      MINNETP%EDGELEN(INVEDGES(1)) = SQRT(SUM( MINNETP%EDGEDIR(INVEDGES(1),:)**2 ))
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%EDGEDIR(INVEDGES(1),:)/MINNETP%EDGELEN(INVEDGES(1))
    ELSE
      MINNETP%EDGESTART(INVEDGES(1),:) = MINNETP%NODEPOS(EN11,:)
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%NODEPOS(NODE1,:)-MINNETP%NODEPOS(EN11,:)
      MINNETP%EDGELEN(INVEDGES(1)) = SQRT(SUM( MINNETP%EDGEDIR(INVEDGES(1),:)**2 ))
      MINNETP%EDGEDIR(INVEDGES(1),:) = MINNETP%EDGEDIR(INVEDGES(1),:)/MINNETP%EDGELEN(INVEDGES(1))
    ENDIF

    ! check if new max edge length
    IF (MINNETP%EDGELEN(INVEDGES(1)).GT.MINNETP%EDGELEN(MINNETP%EIDMAX)) THEN
      MINNETP%EIDMAX = INVEDGES(1)
    ENDIF

    ! invedges(2) updates
    EN21 = MINNETP%EDGENODE(INVEDGES(2),1)
    EN22 = MINNETP%EDGENODE(INVEDGES(2),2)

    IF (EN21.EQ.NODE2) THEN
      MINNETP%EDGESTART(INVEDGES(2),:) = MINNETP%NODEPOS(NODE1,:)
      MINNETP%EDGEDIR(INVEDGES(2),:) = MINNETP%NODEPOS(EN22,:)-MINNETP%NODEPOS(NODE1,:)
      MINNETP%EDGELEN(INVEDGES(2)) = SQRT(SUM( MINNETP%EDGEDIR(INVEDGES(2),:)**2 ))
      MINNETP%EDGEDIR(INVEDGES(2),:) = MINNETP%EDGEDIR(INVEDGES(2),:)/MINNETP%EDGELEN(INVEDGES(2))
      MINNETP%EDGENODE(INVEDGES(2),1) = NODE1
    ELSE
      MINNETP%EDGESTART(INVEDGES(2),:) = MINNETP%NODEPOS(EN21,:)
      MINNETP%EDGEDIR(INVEDGES(2),:) = MINNETP%NODEPOS(NODE1,:)-MINNETP%NODEPOS(EN21,:)
      MINNETP%EDGELEN(INVEDGES(2)) = SQRT(SUM( MINNETP%EDGEDIR(INVEDGES(2),:)**2 ))
      MINNETP%EDGEDIR(INVEDGES(2),:) = MINNETP%EDGEDIR(INVEDGES(2),:)/MINNETP%EDGELEN(INVEDGES(2))
      MINNETP%EDGENODE(INVEDGES(2),2) = NODE1
    ENDIF

    ! check if new max edge length
    IF (MINNETP%EDGELEN(INVEDGES(2)).GT.MINNETP%EDGELEN(MINNETP%EIDMAX)) THEN
      MINNETP%EIDMAX = INVEDGES(2)
    ENDIF

    ! track node
    IF (TRACKNODES) THEN
      ! find all occurrences of NODE2 and update to NODE1
      DO NC=1,MINNETP%NODEON
        IF (MINNETP%NODETRACK(NC).EQ.NODE2) THEN
          MINNETP%NODETRACK(NC) = NODE1
        ENDIF
      ENDDO
    ENDIF

    ! wipe NODE2
    CALL WIPENODE(MINNETP,NODE2)
    ! wipe two short edges
    CALL WIPEEDGE(MINNETP,EIDS(1))
    CALL WIPEEDGE(MINNETP,EIDS(2))

  END SUBROUTINE REMOVELOOP



  SUBROUTINE SWAPNEIGHT1BND(MINNETP,NODE1,NODE2,NODES1,NODES2,TMPNODEPOS)
    ! T1 between node1 and node2, swap specified neighbors
    ! swaps NODES1(1) with NODES2(1)
    ! where NODES1 contains index of NODENODE(NODE1,:) 
    ! example: node1 = 2, node2 = 5. nodes1 = 2 3, nodes2 = 2 3

    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) ::  NODE1, NODE2, NODES1(MAXBRANCH), NODES2(MAXBRANCH)
    DOUBLE PRECISION, INTENT(IN) :: TMPNODEPOS(MINNETP%DIM)
    INTEGER :: NODENODE1(MAXBRANCH), NODEEDGE1(MAXBRANCH), NN1NEIGH(MAXBRANCH), NN2NEIGH(MAXBRANCH)
    INTEGER :: EDGE11, EDGE21, NODE11, NODE21, CT, IND, EID
    DOUBLE PRECISION :: DIST

    ! move node
    CALL MOVENODE(MINNETP,NODE1,TMPNODEPOS)
  
    ! edge11 is edge from node 1 that is changing hands, edge21 is edge from node 2 that is changing hands
    EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
    EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(1))

    IF (MINNETP%EDGENODE(EDGE11,1).EQ.NODE1) THEN
      NODE11 = MINNETP%EDGENODE(EDGE11,2)
    ELSE
      NODE11 = MINNETP%EDGENODE(EDGE11,1)
    ENDIF
    IF (MINNETP%EDGENODE(EDGE21,1).EQ.NODE2) THEN
      NODE21 = MINNETP%EDGENODE(EDGE21,2)
    ELSE
      NODE21 = MINNETP%EDGENODE(EDGE21,1)
    ENDIF

    MINNETP%EDGENODE(EDGE11,1) = NODE2
    MINNETP%EDGENODE(EDGE11,2) = NODE11
    MINNETP%EDGESTART(EDGE11,:) = MINNETP%NODEPOS(NODE2,:)
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%NODEPOS(NODE11,:) - MINNETP%NODEPOS(NODE2,:)
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    MINNETP%EDGELEN(EDGE11) = DIST
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%EDGEDIR(EDGE11,:)/MINNETP%EDGELEN(EDGE11)

    MINNETP%EDGENODE(EDGE21,1) = NODE1
    MINNETP%EDGENODE(EDGE21,2) = NODE21
    MINNETP%EDGESTART(EDGE21,:) = MINNETP%NODEPOS(NODE1,:)
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%NODEPOS(NODE21,:) - MINNETP%NODEPOS(NODE1,:)
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    MINNETP%EDGELEN(EDGE21) = DIST
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%EDGEDIR(EDGE21,:)/MINNETP%EDGELEN(EDGE21)

    ! reset node values
    NODENODE1 = MINNETP%NODENODE(NODE1,:)
    NODEEDGE1 = MINNETP%NODEEDGE(NODE1,:)

    ! set new nodenode values for neighbors of nodes 1 and 2
    NN1NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),:)
    NN2NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),:)

    ! find what was connected to node1, attach to node2, using nodeedges to find correct index
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE1,NODES1(1)),:),EDGE11,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),IND) = NODE2

    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE2,NODES2(1)),:),EDGE21,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),IND) = NODE1

    ! setting new neighbor for node 1
    MINNETP%NODENODE(NODE1,NODES1(1)) = MINNETP%NODENODE(NODE2,NODES2(1))
    MINNETP%NODEEDGE(NODE1,NODES1(1)) = MINNETP%NODEEDGE(NODE2,NODES2(1))

    ! setting new neighbor for node 2
    MINNETP%NODENODE(NODE2,NODES2(1)) = NODENODE1(NODES1(1))
    MINNETP%NODEEDGE(NODE2,NODES2(1)) = NODEEDGE1(NODES1(1))

    ! node 1 is now bnd node
    MINNETP%NODEBND(NODE1) = .TRUE.
    ! edge connecting node 1 and 2 is bnd edge
    IND = FINDLOC(MINNETP%NODENODE(NODE1,:),NODE2,DIM=1)
    EID = MINNETP%NODEEDGE(NODE1,IND)
    MINNETP%EDGEBND(EID) = .TRUE.

    ! if node 11 was a boundary node then make the edge between 11 and node2 a boundary edge
    ! IF (MINNETP%NODEBND(NODE11).AND.NODE1.EQ.416) THEN
    IF (MINNETP%NODEBND(NODE11)) THEN
      MINNETP%EDGEBND(EDGE11) = .TRUE.
    ENDIF

    ! set nodedir (since node1 is now sliding)
    MINNETP%NODEDIR(NODE1,:) = MINNETP%EDGEDIR(EDGE21,:)


  END SUBROUTINE SWAPNEIGHT1BND



  SUBROUTINE SWAPNEIGHT1(MINNETP,NODE1,NODE2,NODES1,NODES2)
    ! T1 between node1 and node2, swap specified neighbors
    ! swaps NODES1(1) with NODES2(1)
    ! where NODES1 contains index of NODENODE(NODE1,:) 
    ! example: node1 = 2, node2 = 5. nodes1 = 2 3, nodes2 = 2 3
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) ::  NODE1, NODE2, NODES1(MAXBRANCH), NODES2(MAXBRANCH)
    INTEGER :: NODENODE1(MAXBRANCH), NODEEDGE1(MAXBRANCH), NN1NEIGH(MAXBRANCH), NN2NEIGH(MAXBRANCH)
    INTEGER :: EDGE11, EDGE21, NODE11, NODE21, CT, IND
    DOUBLE PRECISION :: DIST
  
    ! reset edge dir and edge len for effected edges
    ! edge11 is edge from node 1 that is changing hands, edge21 is edge from node 2 that is changing hands
    EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
    EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(1))

    IF (MINNETP%EDGENODE(EDGE11,1).EQ.NODE1) THEN
      NODE11 = MINNETP%EDGENODE(EDGE11,2)
    ELSE
      NODE11 = MINNETP%EDGENODE(EDGE11,1)
    ENDIF
    IF (MINNETP%EDGENODE(EDGE21,1).EQ.NODE2) THEN
      NODE21 = MINNETP%EDGENODE(EDGE21,2)
    ELSE
      NODE21 = MINNETP%EDGENODE(EDGE21,1)
    ENDIF

    MINNETP%EDGENODE(EDGE11,1) = NODE2
    MINNETP%EDGENODE(EDGE11,2) = NODE11
    MINNETP%EDGESTART(EDGE11,:) = MINNETP%NODEPOS(NODE2,:)
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%NODEPOS(NODE11,:) - MINNETP%NODEPOS(NODE2,:)
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    MINNETP%EDGELEN(EDGE11) = DIST
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%EDGEDIR(EDGE11,:)/MINNETP%EDGELEN(EDGE11)

    MINNETP%EDGENODE(EDGE21,1) = NODE1
    MINNETP%EDGENODE(EDGE21,2) = NODE21
    MINNETP%EDGESTART(EDGE21,:) = MINNETP%NODEPOS(NODE1,:)
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%NODEPOS(NODE21,:) - MINNETP%NODEPOS(NODE1,:)
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    MINNETP%EDGELEN(EDGE21) = DIST
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%EDGEDIR(EDGE21,:)/MINNETP%EDGELEN(EDGE21)

    ! reset node values
    NODENODE1 = MINNETP%NODENODE(NODE1,:)
    NODEEDGE1 = MINNETP%NODEEDGE(NODE1,:)

    ! set new nodenode values for neighbors of nodes 1 and 2
    NN1NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),:)
    NN2NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),:)

    ! find what was connected to node1, attach to node2, using nodeedges to find correct index
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE1,NODES1(1)),:),EDGE11,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),IND) = NODE2

    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE2,NODES2(1)),:),EDGE21,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),IND) = NODE1

    ! setting new neighbor for node 1
    MINNETP%NODENODE(NODE1,NODES1(1)) = MINNETP%NODENODE(NODE2,NODES2(1))
    MINNETP%NODEEDGE(NODE1,NODES1(1)) = MINNETP%NODEEDGE(NODE2,NODES2(1))

    ! setting new neighbor for node 2
    MINNETP%NODENODE(NODE2,NODES2(1)) = NODENODE1(NODES1(1))
    MINNETP%NODEEDGE(NODE2,NODES2(1)) = NODEEDGE1(NODES1(1))

  END SUBROUTINE SWAPNEIGHT1


  SUBROUTINE REMOVENNENTRY(MINNETP,NODE1,NODE2,NODE11,NODES2)
    ! commonly performed removal of a node node entry
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) ::  NODE1, NODE2, NODE11, NODES2(MAXBRANCH)
    INTEGER :: NN1NEIGH(MAXBRANCH)
    INTEGER :: CT, I, IND
    LOGICAL :: REMOVED

    ! set new nodenode values for neighbors of nodes 1 and 2
    NN1NEIGH = MINNETP%NODENODE(NODE11,:)

    ! have to remove an entry from nodenode(node11,:) because it will lose a connection after this loop is fused
    DO CT = 1,MAXBRANCH
      IF (NN1NEIGH(CT).EQ.NODE1) THEN
        IND = CT
        EXIT
      ENDIF
    ENDDO

    ! loop through nodenode and nodedeg and shift accordingly
    REMOVED = .FALSE.
    DO I = 1,MINNETP%NODEDEG(NODE11)
      IF (IND.EQ.I) THEN
        REMOVED=.TRUE.
        CYCLE
      ENDIF
      IF (REMOVED) THEN
        ! removed something lower, so shift everything down one index!
        MINNETP%NODENODE(NODE11,I-1) = MINNETP%NODENODE(NODE11,I)
        MINNETP%NODEEDGE(NODE11,I-1) = MINNETP%NODEEDGE(NODE11,I)
      ENDIF
    ENDDO
    ! after do loop, want to remove unused nodenode and nodeedge values, 
    ! ie. entries past the nodedeg. just to keep these arrays clean
    MINNETP%NODENODE(NODE11,MINNETP%NODEDEG(NODE11)) = 0
    MINNETP%NODEEDGE(NODE11,MINNETP%NODEDEG(NODE11)) = 0
    ! update node degree, lost a connection
    MINNETP%NODEDEG(NODE11) = MINNETP%NODEDEG(NODE11) - 1
 
 
    ! removing neighbor 1 for node 2
    REMOVED=.FALSE.
    DO I = 1,MINNETP%NODEDEG(NODE2)
      IF (NODES2(1).EQ.I) THEN
        REMOVED=.TRUE.
        CYCLE
      ENDIF
      IF (REMOVED) THEN
        MINNETP%NODENODE(NODE2,I-1) = MINNETP%NODENODE(NODE2,I)
        MINNETP%NODEEDGE(NODE2,I-1) = MINNETP%NODEEDGE(NODE2,I)
      ENDIF
    ENDDO
    ! after do loop, want to remove unused nodenode and nodeedge values, 
    ! ie. entries past the nodedeg. just to keep these arrays clean
    MINNETP%NODENODE(NODE2,MINNETP%NODEDEG(NODE2)) = 0
    MINNETP%NODEEDGE(NODE2,MINNETP%NODEDEG(NODE2)) = 0
    ! update node degree, lost a connection
    MINNETP%NODEDEG(NODE2) = MINNETP%NODEDEG(NODE2) - 1


  END SUBROUTINE REMOVENNENTRY



  SUBROUTINE D3ABSORBD2(MINNETP,NODE1,NODE2,NODES1)
    ! node1 = d2, node2 = d3. node 2 absorb node 1 and its connection
    !     *
    !      \
    !       \ 2    1          NODES1
    !        * -- * -------- *
    !       /
    !      /
    !     *
    USE KEYS, ONLY : TRACKNODES
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE1, NODE2, NODES1(1)
    INTEGER :: NODE11, EDGE11, IND, EID, NC
    DOUBLE PRECISION :: DIST

    ! what edge is connected to node1?
    EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
    ! NODES1(1) is either 1 or 2
    ! so EID is edge connecting node 1 and node 2
    EID = MINNETP%NODEEDGE(NODE1,3-NODES1(1))

    IF (MINNETP%EDGENODE(EDGE11,1).EQ.NODE1) THEN
      NODE11 = MINNETP%EDGENODE(EDGE11,2)
    ELSE
      NODE11 = MINNETP%EDGENODE(EDGE11,1)
    ENDIF

    ! new nodepos
    IF ((MINNETP%NODEFIX(NODE2)).OR.MINNETP%NODEBND(NODE2)) THEN
      ! don't need to update position
    ELSEIF (MINNETP%NODEFIX(NODE1)) THEN
      MINNETP%NODEFIX(NODE2) = .TRUE.
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE1,:)
    ELSEIF (MINNETP%NODEPIN(NODE1)) THEN
      MINNETP%NODEPIN(NODE2) = .TRUE.
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE1,:)
    ELSE
      ! average positions
      MINNETP%NODEPOS(NODE2,:) = (MINNETP%NODEPOS(NODE1,:) + MINNETP%NODEPOS(NODE2,:))/2
    ENDIF

    MINNETP%EDGENODE(EDGE11,1) = NODE2
    MINNETP%EDGENODE(EDGE11,2) = NODE11
    MINNETP%EDGESTART(EDGE11,:) = MINNETP%NODEPOS(NODE2,:)
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%NODEPOS(NODE11,:) - MINNETP%NODEPOS(NODE2,:)
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    MINNETP%EDGELEN(EDGE11) = DIST
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%EDGEDIR(EDGE11,:)/MINNETP%EDGELEN(EDGE11)

    ! find what was connected to node1, attach to node2, using nodeedges to find correct index
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE1,NODES1(1)),:),EDGE11,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),IND) = NODE2

    ! update node2 nn,ne
    IND = FINDLOC(MINNETP%NODEEDGE(NODE2,:),EID,DIM=1)
    MINNETP%NODENODE(NODE2,IND) = NODE11
    MINNETP%NODEEDGE(NODE2,IND) = EDGE11

    ! track node
    IF (TRACKNODES) THEN
      ! find all occurrences of NODE1 and update to NODE2
      DO NC=1,MINNETP%NODEON
        IF (MINNETP%NODETRACK(NC).EQ.NODE1) THEN
          MINNETP%NODETRACK(NC) = NODE2
        ENDIF
      ENDDO
    ENDIF

    ! deactivate node1
    CALL WIPENODE(MINNETP,NODE1)

    ! deactivate eid
    CALL WIPEEDGE(MINNETP,EID)


  END SUBROUTINE D3ABSORBD2


  SUBROUTINE D2ABSORBD2(MINNETP,NODE1,NODE2)
    ! node1 = d2, node2 = d2. node 2 absorb node 1 and its connection
    USE KEYS, ONLY : MAXBRANCH, TRACKNODES
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE1, NODE2
    INTEGER :: NODES1(MAXBRANCH), NODES2(MAXBRANCH), EDGE11, IND
    INTEGER :: EID, TMPNODE, NC, EC, CT, ECC, TMPNODES2(MAXBRANCH), NODE11, NODE21, NODE22, INVEDGES(2*MAXBRANCH-1)
    DOUBLE PRECISION :: TOTLEN, NEWLEN1, NEWLEN2, DIST
    INTEGER :: NN1(2), NN2(2), EIDS(2), CT1, CT2

    INVEDGES = 0
    EIDS = 0
    ! first check if this involves a degenerate loop that has shrunk to below DX
    NN1 = MINNETP%NODENODE(NODE1,1:2)
    NN2 = MINNETP%NODENODE(NODE2,1:2)
    NODES1=0
    NODES2=0

    CT = 0
    ! these keep track of how many times each node appears in the other's NODENODE
    CT1 = 0
    CT2 = 0
    ! find their neighbors, and keep track of "involved edges" INVEDGES
    DO NC = 1,2
      CT = CT + 1
      TMPNODE = NN1(NC)
      IF (TMPNODE.EQ.NODE2) THEN
        CT1 = CT1+1
        ! track which edges are between node 1 and 2 (especially if there are two!)
        EIDS(CT1) = MINNETP%NODEEDGE(NODE1,NC)
        CYCLE
      ENDIF
      CT = CT - 1
      INVEDGES(NC-CT) = MINNETP%NODEEDGE(NODE1,NC)
      NODES1(NC-CT) = NC
    ENDDO
    CT = 0
    DO NC = 1,2
      CT = CT + 1
      TMPNODE = NN2(NC)
      IF (TMPNODE.EQ.NODE1) THEN
        CT2 = CT2+1
        CYCLE
      ENDIF
      CT = CT - 1
      NODES2(NC-CT) = NC
      INVEDGES(2 - CT1 + NC-CT) = MINNETP%NODEEDGE(NODE2,NC)
    ENDDO

    IF (CT1.NE.CT2) THEN
      PRINT*, 'ERROR IN D2ABSORBD2: discrepancy in nodenode for nodes ', NODE1, NODE2
      PRINT*, NN1
      PRINT*, NN2
      STOP 1
    ENDIF

    ! if there is a degenerate edge, something is wrong
    IF (CT1.EQ.2) THEN
      ! check if ct1=ct2
      PRINT*, 'ERROR IN D2ABSORBD2: separate network with 2 d2 nodes, shrank to a point'
      STOP 1
    ENDIF

    ! new nodepos
    IF (MINNETP%NODEFIX(NODE2).OR.MINNETP%NODEPIN(NODE2)) THEN
      ! don't need to update position
    ELSEIF (MINNETP%NODEFIX(NODE1)) THEN
      MINNETP%NODEFIX(NODE2) = .TRUE.
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE1,:)
    ELSEIF (MINNETP%NODEPIN(NODE1)) THEN
      MINNETP%NODEPIN(NODE2) = .TRUE.
      MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE1,:)
    ELSE
      ! average positions
      MINNETP%NODEPOS(NODE2,:) = (MINNETP%NODEPOS(NODE1,:) + MINNETP%NODEPOS(NODE2,:))/2
    ENDIF
   
    ! what edge is connected to node1?
    EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
    ! NODES1(1) is either 1 or 2
    ! so EID is edge connecting node 1 and node 2
    EID = MINNETP%NODEEDGE(NODE1,3-NODES1(1))

    IF (MINNETP%EDGENODE(EDGE11,1).EQ.NODE1) THEN
      NODE11 = MINNETP%EDGENODE(EDGE11,2)
    ELSE
      NODE11 = MINNETP%EDGENODE(EDGE11,1)
    ENDIF

    MINNETP%EDGENODE(EDGE11,1) = NODE2
    MINNETP%EDGENODE(EDGE11,2) = NODE11
    MINNETP%EDGESTART(EDGE11,:) = MINNETP%NODEPOS(NODE2,:)
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%NODEPOS(NODE11,:) - MINNETP%NODEPOS(NODE2,:)
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    MINNETP%EDGELEN(EDGE11) = DIST
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%EDGEDIR(EDGE11,:)/MINNETP%EDGELEN(EDGE11)

    ! find what was connected to node1, attach to node2, using nodeedges to find correct index
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE1,NODES1(1)),:),EDGE11,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),IND) = NODE2

    ! update node2 nn,ne
    IND = FINDLOC(MINNETP%NODEEDGE(NODE2,:),EID,DIM=1)
    MINNETP%NODENODE(NODE2,IND) = NODE11
    MINNETP%NODEEDGE(NODE2,IND) = EDGE11

    ! track node
    IF (TRACKNODES) THEN
      ! find all occurrences of NODE1 and update to NODE2
      DO NC=1,MINNETP%NODEON
        IF (MINNETP%NODETRACK(NC).EQ.NODE1) THEN
          MINNETP%NODETRACK(NC) = NODE2
        ENDIF
      ENDDO
    ENDIF

    ! deactivate node1
    CALL WIPENODE(MINNETP,NODE1)

    ! deactivate eid
    CALL WIPEEDGE(MINNETP,EID)


  END SUBROUTINE D2ABSORBD2



  SUBROUTINE PASSEVENT(MINNETP,NODE1,NODE2)
    ! passing event
    ! try alternate arrangement of connections
    ! always set up so that first node input is deg2, second is deg3. total degree will always be 5, with 3 invedges and the 1 in
    ! between two nodes of interest
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE1, NODE2
    INTEGER :: NODES1(MAXBRANCH), NODES2(MAXBRANCH)
    INTEGER :: EID, TMPNODE, NC, EC, CT, ECC, TMPNODES2(MAXBRANCH), NODE11, NODE21, NODE22, INVEDGES(2*MAXBRANCH-1)
    DOUBLE PRECISION :: TOTLEN, NEWLEN1, NEWLEN2, DIST
    INTEGER :: NN1(2), NN2(3), EIDS(2), CT1, CT2, EDGE21, EDGE11
    DOUBLE PRECISION :: NEWNP2(MINNETP%DIM), V1(MINNETP%DIM), V2(MINNETP%DIM)

    INVEDGES = 0
    EIDS = 0
    ! first check if this involves a degenerate loop that has shrunk to below DX
    NN1 = MINNETP%NODENODE(NODE1,1:2)
    NN2 = MINNETP%NODENODE(NODE2,1:3)
    NODES1=0
    NODES2=0

    CT = 0
    ! these keep track of how many times each node appears in the other's NODENODE
    CT1 = 0
    CT2 = 0
    ! find their neighbors, and keep track of "involved edges" INVEDGES
    DO NC = 1,2
      CT = CT + 1
      TMPNODE = NN1(NC)
      IF (TMPNODE.EQ.NODE2) THEN
        CT1 = CT1+1
        ! track which edges are between node 1 and 2 (especially if there are two!)
        EIDS(CT1) = MINNETP%NODEEDGE(NODE1,NC)
        CYCLE
      ENDIF
      CT = CT - 1
      INVEDGES(NC-CT) = MINNETP%NODEEDGE(NODE1,NC)
      NODES1(NC-CT) = NC
    ENDDO
    CT = 0
    DO NC = 1,3
      CT = CT + 1
      TMPNODE = NN2(NC)
      IF (TMPNODE.EQ.NODE1) THEN
        CT2 = CT2+1
        CYCLE
      ENDIF
      CT = CT - 1
      NODES2(NC-CT) = NC
      INVEDGES(2 - CT1 + NC-CT) = MINNETP%NODEEDGE(NODE2,NC)
    ENDDO

    IF (CT1.NE.CT2) THEN
      PRINT*, 'ERROR IN PASSEVENT: discrepancy in nodenode for nodes ', NODE1, NODE2
      PRINT*, NN1
      PRINT*, NN2
      STOP 1
    ENDIF


    ! if there is a degenerate edge, and it is the short edge, remove the loop and done with this call
    IF (CT1.EQ.2) THEN
      ! check if ct1=ct2

      IF (VERBOSE) THEN
        PRINT*, 'passevent loop fusing'
        PRINT*, INVEDGES
        PRINT*, EIDS
        PRINT*, NODES1,NODES2
      ENDIF
      ! have to rearrange invedges since call to removeloop uses first input node as future node,
      CALL REMOVELOOP23(MINNETP,NODE1,NODE2,NODES2(1),EIDS,INVEDGES)
      RETURN
    ENDIF

    ! already ruled out double connected case, now checking for boundary edges!
    IF (ANY(MINNETP%EDGEBND(INVEDGES(1:3)))) THEN
      IF (VERBOSE) THEN
        PRINT*, 'found a PA on boundary'
      ENDIF
      ! have to be careful here, must move the d3 node to be in correct place along next boundary edge

      ! find edge21, defined to be edge from node 2 that is boundary edge, but not the one between node1 and node2
      IF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE2,NODES2(1)))) THEN
        EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(1))
        TMPNODES2 = NODES2
      ELSEIF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE2,NODES2(2)))) THEN
        EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(2))
        TMPNODES2(2) = NODES2(1)
        TMPNODES2(1) = NODES2(2)
      ELSE
        PRINT*, 'ERROR IN PASSEVENT, not enough boundary edges, might be just a d3 node passing close to d2 bnd'
        ! STOP 1
        RETURN
      ENDIF

      ! check if d2 node is fixed, if it is, it is a fixed boundary node, and you can continue, if not, absorb it
      IF (.NOT.(MINNETP%NODEFIX(NODE1))) THEN
        CALL D3ABSORBD2(MINNETP,NODE1,NODE2,NODES1(1))
        RETURN
      ENDIF

      TOTLEN = 0D0
      DO ECC = 1,3
        TOTLEN = TOTLEN + MINNETP%EDGELEN(INVEDGES(ECC))
      ENDDO

      ! counting edgelen between nodes since it could change if node2 is projected onto the boundary
      TOTLEN = TOTLEN + MINNETP%EDGELEN(EIDS(1))


      ! start with position of NODE2, must be placed along next boundary edge
      ! dot product to project NODE2 onto edge
      V1 = MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE1,:)
      EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
      V2 = MINNETP%EDGEDIR(EDGE11,:)
      NEWNP2 = MINNETP%NODEPOS(NODE1,:) + DOT_PRODUCT(V1,V2)*V2

      ! calculate lengths
      NODE11 = MINNETP%NODENODE(NODE1,NODES1(1))
      NODE21 = MINNETP%NODENODE(NODE2,TMPNODES2(1))
      NODE22 = MINNETP%NODENODE(NODE2,TMPNODES2(2))

      ! now see if passing event lowers total length

      ! edge between n1 n2
      NEWLEN1 = SQRT(SUM((NEWNP2 - MINNETP%NODEPOS(NODE1,:))**2))
      NEWLEN1 = NEWLEN1 + SQRT(SUM((NEWNP2 - MINNETP%NODEPOS(NODE11,:))**2))
      NEWLEN1 = NEWLEN1 + SQRT(SUM((NEWNP2 - MINNETP%NODEPOS(NODE22,:))**2))
      NEWLEN1 = NEWLEN1 + SQRT(SUM((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE21,:))**2))

      IF (VERBOSE) THEN
        PRINT*, 'boundary passing oldlength', TOTLEN
        PRINT*, 'boundary passing newlength', NEWLEN1
      ENDIF

      IF (NEWLEN1.LE.TOTLEN) THEN
        CALL SWAPNEIGHPABND(MINNETP,NODE1,NODE2,NODES1,TMPNODES2)
      ENDIF

      RETURN
    ENDIF

    ! if the d2 node is not fixed, then have d3 node absorb d2 node
    ! if d3 is fixed, should absorb since meaningless d2 node
    ! if d3 is not fixed, but d2 isn't either, also absorb
      ! also confirm that this is not a pinned node, in which case don't absorb it
    IF ((.NOT.MINNETP%NODEFIX(NODE1)).AND.(.NOT.MINNETP%NODEPIN(NODE1))) THEN
      CALL D3ABSORBD2(MINNETP,NODE1,NODE2,NODES1(1))
      RETURN
    ENDIF

    TOTLEN = 0D0
    DO ECC = 1,3
      TOTLEN = TOTLEN + MINNETP%EDGELEN(INVEDGES(ECC))
    ENDDO

    NODE11 = MINNETP%NODENODE(NODE1,NODES1(1))
    NODE21 = MINNETP%NODENODE(NODE2,NODES2(1))
    NODE22 = MINNETP%NODENODE(NODE2,NODES2(2))

    ! now see if passing event lowers total length
    NEWLEN1 = 0D0
    NEWLEN2 = 0D0

    ! add up lengths, two possible combinations:

    ! combo 1
    CALL NODENODEDIST(MINNETP,NODE2,MINNETP%NODENODE(NODE1,NODES1(1)),DIST)
    NEWLEN1 = NEWLEN1 + DIST
    CALL NODENODEDIST(MINNETP,NODE1,MINNETP%NODENODE(NODE2,NODES2(1)),DIST)
    NEWLEN1 = NEWLEN1 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,MINNETP%NODENODE(NODE2,NODES2(2)),DIST)
    NEWLEN1 = NEWLEN1 + DIST
    ! don't have to include edge in between two nodes of interest, the short one, because it is in both estimates

    ! combo 2
    CALL NODENODEDIST(MINNETP,NODE2,MINNETP%NODENODE(NODE1,NODES1(1)),DIST)
    NEWLEN2 = NEWLEN2 + DIST
    CALL NODENODEDIST(MINNETP,NODE1,MINNETP%NODENODE(NODE2,NODES2(2)),DIST)
    NEWLEN2 = NEWLEN2 + DIST
    CALL NODENODEDIST(MINNETP,NODE2,MINNETP%NODENODE(NODE2,NODES2(1)),DIST)
    NEWLEN2 = NEWLEN2 + DIST
    ! don't have to include edge in between two nodes of interest, the short one, because it is in both estimates

    IF (VERBOSE) THEN
      PRINT*, 'total length before rearrangement', TOTLEN
      PRINT*, 'new length after rearrangement 1', NEWLEN1
      PRINT*, 'new length after rearrangement 2', NEWLEN2
    ENDIF

    IF ((NEWLEN1.LT.TOTLEN).AND.(NEWLEN1.LE.NEWLEN2)) THEN
      ! combo 1

      CALL SWAPNEIGHPA(MINNETP,NODE1,NODE2,NODES1,NODES2)

    ELSEIF (NEWLEN2.LT.TOTLEN) THEN
      ! combo 2
        
      TMPNODES2(1) = NODES2(2)
      TMPNODES2(2) = NODES2(1)

      CALL SWAPNEIGHPA(MINNETP,NODE1,NODE2,NODES1,TMPNODES2)
    ENDIF

  END SUBROUTINE PASSEVENT


  SUBROUTINE PAFUSELOOP(MINNETP,NODE1,NODE2,NODES1,NODES2)
    ! Passing event between node1 and node2, but delete degenerate edges so that node2 becomes degree 2
    ! NOTE: assumes NODE1 is the degree 2 node, NODE2 degree 3, make sure called correctly
    ! swaps NODES1(1) with NODES2(1)
    ! assumes node11=node22
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) ::  NODE1, NODE2, NODES1(MAXBRANCH), NODES2(MAXBRANCH)
    INTEGER :: NODENODE1(MAXBRANCH), NODEEDGE1(MAXBRANCH), NN1NEIGH(MAXBRANCH), NN2NEIGH(MAXBRANCH)
    INTEGER :: EDGE11, EDGE21, NODE11, NODE21, CT, I, IND
    DOUBLE PRECISION :: DIST
    LOGICAL :: REMOVED
  
    ! reset edge dir and edge len for effected edges
    ! edge11 is edge from node 1 that is changing hands, edge21 is edge from node 2 that is changing hands
    EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
    EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(1))

    NODE21 = MINNETP%NODENODE(NODE2,NODES2(1))
 
    ! edge11 is deemed inactive and it no longer connects anything
    CALL WIPEEDGE(MINNETP,EDGE11)
 
    MINNETP%EDGENODE(EDGE21,1) = NODE1
    MINNETP%EDGENODE(EDGE21,2) = NODE21
    MINNETP%EDGESTART(EDGE21,:) = MINNETP%NODEPOS(NODE1,:)
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%NODEPOS(NODE21,:) - MINNETP%NODEPOS(NODE1,:)
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    MINNETP%EDGELEN(EDGE21) = DIST
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%EDGEDIR(EDGE21,:)/MINNETP%EDGELEN(EDGE21)

    ! reset node values
    NODENODE1 = MINNETP%NODENODE(NODE1,:)
    NODEEDGE1 = MINNETP%NODEEDGE(NODE1,:)

    ! commonly used values
    NODE11 = MINNETP%NODENODE(NODE1,NODES1(1))
    NODE21 = MINNETP%NODENODE(NODE2,NODES2(1))
 
    ! set new nodenode values for neighbors of nodes 1 and 2
    NN2NEIGH = MINNETP%NODENODE(NODE21,:)

    CALL REMOVENNENTRY(MINNETP,NODE1,NODE2,NODE11,NODES2)

    ! update nodenode for node21, connect to node1
    DO CT = 1,MAXBRANCH
      IF (NN2NEIGH(CT).EQ.NODE2) THEN
        IND = CT
        EXIT
      ENDIF
    ENDDO
    MINNETP%NODENODE(NODE21,IND) = NODE1
 
    ! setting new neighbor for node 1
    MINNETP%NODENODE(NODE1,NODES1(1)) = NODE21
    MINNETP%NODEEDGE(NODE1,NODES1(1)) = EDGE21

  END SUBROUTINE PAFUSELOOP


  SUBROUTINE SWAPNEIGHPABND(MINNETP,NODE1,NODE2,NODES1,TMPNODES2)
    ! Passing event between node1 and node2, swap their neighbors
    ! NOTE: assumes NODE1 is the degree 2 node, NODE2 degree 3, make sure called correctly
    ! This is the special case of NODE1 is a fixed node along the boundary, NODE2 is a
    ! sliding node on a boundary edge
    ! must swap neighbors carefully, and adjust position of NODE2 to be on neighboring edge
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) ::  NODE1, NODE2, NODES1(MAXBRANCH), TMPNODES2(MAXBRANCH)
    INTEGER :: NODENODE1(MAXBRANCH), NODEEDGE1(MAXBRANCH), NN1NEIGH(MAXBRANCH), NN2NEIGH(MAXBRANCH)
    INTEGER :: EDGE11, EDGE21, NODE11, NODE21, CT, IND, EID, NODES2(MAXBRANCH), NID
    DOUBLE PRECISION :: DIST
    DOUBLE PRECISION :: V1(MINNETP%DIM), V2(MINNETP%DIM)

    ! find edge21, defined to be edge from node 2 that is boundary edge, but not the one between node1 and node2
    IF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE2,TMPNODES2(1)))) THEN
      EDGE21 = MINNETP%NODEEDGE(NODE2,TMPNODES2(1))
      NODES2 = TMPNODES2
    ELSEIF (MINNETP%EDGEBND(MINNETP%NODEEDGE(NODE2,TMPNODES2(2)))) THEN
      EDGE21 = MINNETP%NODEEDGE(NODE2,TMPNODES2(2))
      NODES2(1) = TMPNODES2(2)
      NODES2(2) = TMPNODES2(1)
    ELSE
      PRINT*, 'ERROR IN SWAPNEIGHPABND, not enough boundary edges'
    ENDIF

    ! start with position of NODE2, must be placed along next boundary edge
    ! dot product to project NODE2 onto edge
    V1 = MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NODE1,:)
    EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
    V2 = MINNETP%EDGEDIR(EDGE11,:)
    MINNETP%NODEPOS(NODE2,:) = MINNETP%NODEPOS(NODE1,:) + DOT_PRODUCT(V1,V2)*V2

    ! fix edgelen, edgedir, etc of node we just shifted!
    EID = MINNETP%NODEEDGE(NODE2,NODES2(2))
    ! if node2 is the start of edge, then:
    IF (MINNETP%EDGENODE(EID,1).EQ.NODE2) THEN
      NID = MINNETP%EDGENODE(EID,2)
      MINNETP%EDGESTART(EID,:) = MINNETP%NODEPOS(NODE2,:)
      MINNETP%EDGEDIR(EID,:) = MINNETP%NODEPOS(NID,:) - MINNETP%NODEPOS(NODE2,:)
      CALL NODENODEDIST(MINNETP,NODE2,NID,DIST)
      MINNETP%EDGELEN(EID) = DIST
      MINNETP%EDGEDIR(EID,:) = MINNETP%EDGEDIR(EID,:)/MINNETP%EDGELEN(EID)
    ELSE
      NID = MINNETP%EDGENODE(EID,1)
      MINNETP%EDGESTART(EID,:) = MINNETP%NODEPOS(NID,:)
      MINNETP%EDGEDIR(EID,:) = MINNETP%NODEPOS(NODE2,:) - MINNETP%NODEPOS(NID,:)
      CALL NODENODEDIST(MINNETP,NID,NODE2,DIST)
      MINNETP%EDGELEN(EID) = DIST
      MINNETP%EDGEDIR(EID,:) = MINNETP%EDGEDIR(EID,:)/MINNETP%EDGELEN(EID)
    ENDIF

    ! reset edge dir and edge len for effected edges
    IF (MINNETP%EDGENODE(EDGE11,1).EQ.NODE1) THEN
      NODE11 = MINNETP%EDGENODE(EDGE11,2)
    ELSE
      NODE11 = MINNETP%EDGENODE(EDGE11,1)
    ENDIF
    IF (MINNETP%EDGENODE(EDGE21,1).EQ.NODE2) THEN
      NODE21 = MINNETP%EDGENODE(EDGE21,2)
    ELSE
      NODE21 = MINNETP%EDGENODE(EDGE21,1)
    ENDIF


    ! Do Not update edgedir here, may lead to compounding errors over a long time
    MINNETP%EDGENODE(EDGE11,1) = NODE2
    MINNETP%EDGENODE(EDGE11,2) = NODE11
    MINNETP%EDGESTART(EDGE11,:) = MINNETP%NODEPOS(NODE2,:)
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    MINNETP%EDGELEN(EDGE11) = DIST

    MINNETP%EDGENODE(EDGE21,1) = NODE1
    MINNETP%EDGENODE(EDGE21,2) = NODE21
    MINNETP%EDGESTART(EDGE21,:) = MINNETP%NODEPOS(NODE1,:)
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    MINNETP%EDGELEN(EDGE21) = DIST
    ! NODES1 or 2 contain indices of nodenode(NODE1,:) that aren't the other swapping node. i.e.
    ! if nodes 41 and 42 are really close and they are degree 2 and 3, node1=41, node2=42, and 
    ! we could have nodes1=1, nodes2 = 2 3, showing which of nodenode the other neighbors are. 
    ! In this function we swap nodes2(1) and nodes1(1)
    ! Must update their nodenode as well

    NODENODE1 = MINNETP%NODENODE(NODE1,:)
    NODEEDGE1 = MINNETP%NODEEDGE(NODE1,:)

    ! set new nodenode values for neighbors of nodes 1 and 2
    NN1NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),:)
    NN2NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),:)

    ! find what was connected to node1, attach to node2, using nodeedges to find correct index
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE1,NODES1(1)),:),EDGE11,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),IND) = NODE2
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE2,NODES2(1)),:),EDGE21,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),IND) = NODE1

    ! setting new neighbor for node 1
    MINNETP%NODENODE(NODE1,NODES1(1)) = MINNETP%NODENODE(NODE2,NODES2(1))
    MINNETP%NODEEDGE(NODE1,NODES1(1)) = MINNETP%NODEEDGE(NODE2,NODES2(1))

    ! setting new neighbor for node 2
    MINNETP%NODENODE(NODE2,NODES2(1)) = NODENODE1(NODES1(1))
    MINNETP%NODEEDGE(NODE2,NODES2(1)) = NODEEDGE1(NODES1(1))

    MINNETP%NODEDIR(NODE2,:) = MINNETP%EDGEDIR(EDGE11,:)

  END SUBROUTINE SWAPNEIGHPABND


    
  SUBROUTINE SWAPNEIGHPA(MINNETP,NODE1,NODE2,NODES1,NODES2)
    ! Passing event between node1 and node2, swap their neighbors
    ! NOTE: assumes NODE1 is the degree 2 node, NODE2 degree 3, make sure called correctly
    ! swaps NODES1(1) with NODES2(1)
    USE KEYS, ONLY : MAXBRANCH
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) ::  NODE1, NODE2, NODES1(MAXBRANCH), NODES2(MAXBRANCH)
    INTEGER :: NODENODE1(MAXBRANCH), NODEEDGE1(MAXBRANCH), NN1NEIGH(MAXBRANCH), NN2NEIGH(MAXBRANCH)
    INTEGER :: EDGE11, EDGE21, NODE11, NODE21, CT, IND
    DOUBLE PRECISION :: DIST
  
    ! reset edge dir and edge len for effected edges
    ! edge11 is edge from node 1 that is changing hands, edge21 is edge from node 2 that is changing hands
    EDGE11 = MINNETP%NODEEDGE(NODE1,NODES1(1))
    EDGE21 = MINNETP%NODEEDGE(NODE2,NODES2(1))

    IF (MINNETP%EDGENODE(EDGE11,1).EQ.NODE1) THEN
      NODE11 = MINNETP%EDGENODE(EDGE11,2)
    ELSE
      NODE11 = MINNETP%EDGENODE(EDGE11,1)
    ENDIF
    IF (MINNETP%EDGENODE(EDGE21,1).EQ.NODE2) THEN
      NODE21 = MINNETP%EDGENODE(EDGE21,2)
    ELSE
      NODE21 = MINNETP%EDGENODE(EDGE21,1)
    ENDIF

    MINNETP%EDGENODE(EDGE11,1) = NODE2
    MINNETP%EDGENODE(EDGE11,2) = NODE11
    MINNETP%EDGESTART(EDGE11,:) = MINNETP%NODEPOS(NODE2,:)
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%NODEPOS(NODE11,:) - MINNETP%NODEPOS(NODE2,:)
    CALL NODENODEDIST(MINNETP,NODE2,NODE11,DIST)
    MINNETP%EDGELEN(EDGE11) = DIST
    MINNETP%EDGEDIR(EDGE11,:) = MINNETP%EDGEDIR(EDGE11,:)/MINNETP%EDGELEN(EDGE11)

    MINNETP%EDGENODE(EDGE21,1) = NODE1
    MINNETP%EDGENODE(EDGE21,2) = NODE21
    MINNETP%EDGESTART(EDGE21,:) = MINNETP%NODEPOS(NODE1,:)
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%NODEPOS(NODE21,:) - MINNETP%NODEPOS(NODE1,:)
    CALL NODENODEDIST(MINNETP,NODE1,NODE21,DIST)
    MINNETP%EDGELEN(EDGE21) = DIST
    MINNETP%EDGEDIR(EDGE21,:) = MINNETP%EDGEDIR(EDGE21,:)/MINNETP%EDGELEN(EDGE21)

    ! NODES1 or 2 contain indices of nodenode(NODE1,:) that aren't the other swapping node. i.e.
    ! if nodes 41 and 42 are really close and they are degree 2 and 3, node1=41, node2=42, and 
    ! we could have nodes1=1, nodes2 = 2 3, showing which of nodenode the other neighbors are. 
    ! In this function we swap nodes2(1) and nodes1(1)
    ! Must update their nodenode as well

    NODENODE1 = MINNETP%NODENODE(NODE1,:)
    NODEEDGE1 = MINNETP%NODEEDGE(NODE1,:)

    ! set new nodenode values for neighbors of nodes 1 and 2
    NN1NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),:)
    NN2NEIGH = MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),:)

    ! find what was connected to node1, attach to node2, using nodeedges to find correct index
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE1,NODES1(1)),:),EDGE11,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE1,NODES1(1)),IND) = NODE2
    IND = FINDLOC(MINNETP%NODEEDGE(MINNETP%NODENODE(NODE2,NODES2(1)),:),EDGE21,DIM=1)
    MINNETP%NODENODE(MINNETP%NODENODE(NODE2,NODES2(1)),IND) = NODE1

    ! setting new neighbor for node 1
    MINNETP%NODENODE(NODE1,NODES1(1)) = MINNETP%NODENODE(NODE2,NODES2(1))
    MINNETP%NODEEDGE(NODE1,NODES1(1)) = MINNETP%NODEEDGE(NODE2,NODES2(1))

    ! setting new neighbor for node 2
    MINNETP%NODENODE(NODE2,NODES2(1)) = NODENODE1(NODES1(1))
    MINNETP%NODEEDGE(NODE2,NODES2(1)) = NODEEDGE1(NODES1(1))

  END SUBROUTINE SWAPNEIGHPA

    
  SUBROUTINE NETWORKCHECK(MINNETP, FLAG)
    ! network self consistency checks
    ! returns PASSFAIL = FALSE if failed, TRUE if passed
    ! if failed, prints relevant information

    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    LOGICAL, INTENT(OUT) :: FLAG

  END SUBROUTINE NETWORKCHECK


  SUBROUTINE GETPROJ(MINNETP,NODE,EDGE,NODEPOS)
    ! projection of node onto given edge
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NODE, EDGE
    DOUBLE PRECISION, INTENT(OUT) :: NODEPOS(MINNETP%DIM)
    DOUBLE PRECISION :: POS(MINNETP%DIM), EDIR(MINNETP%DIM), TMP(MINNETP%DIM)

    POS = MINNETP%NODEPOS(NODE,:)
    EDIR = MINNETP%EDGEDIR(EDGE,:)

    TMP = MINNETP%NODEPOS(NODE,:)-MINNETP%EDGESTART(EDGE,:)

    NODEPOS = MINNETP%EDGESTART(EDGE,:) + DOT_PRODUCT(TMP,MINNETP%EDGEDIR(EDGE,:)) &
      & *MINNETP%EDGEDIR(EDGE,:)
    
    ! NODEPOS = (POS(1) - MINNETP%EDGESTART(EDGE,1))*EDIR(1) - &
      ! & (POS(2) - MINNETP%EDGESTART(EDGE,2))*EDIR(2)
    ! NODEPOS = NODEPOS*MINNETP%EDGEDIR(EDGE,:)

  END SUBROUTINE GETPROJ


  SUBROUTINE NODENODEDIST(MINNETP,NODE1,NODE2,DIST)

    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) ::  NODE1, NODE2
    DOUBLE PRECISION, INTENT(OUT) :: DIST
    DOUBLE PRECISION :: DIFFSQ(MINNETP%DIM)

    DIFFSQ = ((MINNETP%NODEPOS(NODE1,:) - MINNETP%NODEPOS(NODE2,:))**2)
    DIST = SQRT(SUM(DIFFSQ))

  END SUBROUTINE NODENODEDIST


  SUBROUTINE SETUPPINNODES(MINNETP)
    ! set up pinned nodes on a network
    ! either from a given list, or set randomly
    USE mt19937, ONLY : GRND
    USE KEYS, ONLY : NPIN, PINNODES, RANDPINNODES
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER :: A, CC, TRY, NID
    INTEGER, PARAMETER :: MAXNTRY= 1000
  
    IF (.NOT.MINNETP%ARRAYSET) THEN
       PRINT*, 'ERROR IN SETPINNODES: network not set up yet'
       STOP 1    
    ENDIF
  
    IF (RANDPINNODES) THEN
       PRINT*, 'setting pinned nodes'
       ! set random fixed nodes
       MINNETP%NODEPIN(:) = .FALSE.
       ! Sample fix nodes for this field
       DO CC = 1,NPIN
          DO TRY = 1,MAXNTRY ! try picking a new node
             NID = FLOOR(GRND()*MINNETP%NODEHIGH)+1
             ! don't try to pin non-existent nodes, pinned, boundary, fixed or growing nodes
             IF (NID.GT.MINNETP%NODEHIGH) THEN
               CYCLE
             ELSEIF (MINNETP%NODEPIN(NID)) THEN
               CYCLE
             ELSEIF (MINNETP%NODEBND(NID)) THEN
               CYCLE
             ELSEIF (MINNETP%NODEFIX(NID)) THEN
               CYCLE
             ELSEIF (MINNETP%NODEGROW(NID).GT.0) THEN
               CYCLE
             ELSE
               EXIT
             ENDIF
          ENDDO
          IF (TRY.GE.MAXNTRY) THEN
             PRINT*, 'warning when setting pinned nodes: continuing simulation with fewer pinned nodes than desired'
             PRINT*, 'active nodes: ', MINNETP%NODEHIGH
             PRINT*, 'Did not perform nth pin, n=', CC
          ELSE
            PRINT*, 'PINNED NODE:', NID
            MINNETP%NODEPIN(NID) = .TRUE.
          ENDIF
       ENDDO
    ELSE ! set the desired fixed nodes
       MINNETP%NODEPIN = .FALSE.
       DO A = 1,NPIN
          IF (PINNODES(A).GT.MINNETP%NODEHIGH.OR.PINNODES(A).LT.1) THEN
             PRINT*, 'ERROR IN SETUPPINNODES: bad pinned node', A, PINNODES(A)
             STOP 1
          ENDIF
          MINNETP%NODEPIN(PINNODES(A)) = .TRUE.             
       ENDDO
    ENDIF

  END SUBROUTINE SETUPPINNODES

  
  SUBROUTINE SETUPFIXNODES(MINNETP)
    ! set up fixed nodes on a network
    ! either from a given list, or set randomly
    USE mt19937, ONLY : GRND
    USE KEYS, ONLY : NFIX, FIXNODES, RANDFIXNODES
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    !INTEGER, INTENT(IN) :: NFIX(:)
    !INTEGER, INTENT(IN), OPTIONAL :: FIXNODES(:,:)
    INTEGER :: A, CC, TRY, NID
    INTEGER, PARAMETER :: MAXNTRY= 100
  
    IF (.NOT.MINNETP%ARRAYSET) THEN
       PRINT*, 'ERROR IN SETFIXNODES: network not set up yet'
       STOP 1    
    ENDIF
  
    IF (RANDFIXNODES) THEN
       ! set random fixed nodes
       MINNETP%NODEFIX(:) = .FALSE.
       ! Sample fix nodes for this field
       DO CC = 1,NFIX
          DO TRY = 1,MAXNTRY ! try picking a new node
             NID = FLOOR(GRND()*MINNETP%NODEHIGH)+1
             IF (NID.GT.MINNETP%NODEHIGH) THEN
               CYCLE
             ELSEIF (MINNETP%NODEFIX(NID)) THEN
               CYCLE
             ELSE
               EXIT
             ENDIF
          ENDDO
          IF (TRY.GE.MAXNTRY) THEN
             PRINT*, 'ERROR IN SETTING FIXED NODES:   failed to find a new node', CC, NID
             STOP 1
          ENDIF

          PRINT*, 'FIXED NODE:', NID
          MINNETP%NODEFIX(NID) = .TRUE.
       ENDDO
    ELSE ! set the desired fixed nodes
       MINNETP%NODEFIX = .FALSE.
       DO A = 1,NFIX
          IF (FIXNODES(A).GT.MINNETP%NODEHIGH.OR.FIXNODES(A).LT.1) THEN
             PRINT*, 'ERROR IN SETUPFIXNODES: bad fixed node', A, FIXNODES(A)
             STOP 1
          ENDIF
          MINNETP%NODEFIX(FIXNODES(A)) = .TRUE.             
       ENDDO
    ENDIF

  END SUBROUTINE SETUPFIXNODES
  

  SUBROUTINE GETMINNETFORCE(MINNETP,FORCES)
    ! loop through active nodes and get forces
    USE KEYS, ONLY : MAXBRANCH, MOB
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    DOUBLE PRECISION, INTENT(OUT) :: FORCES(MINNETP%NODEHIGH,MINNETP%DIM)
    DOUBLE PRECISION, DIMENSION(MINNETP%NODEHIGH,MINNETP%DIM) :: FTENS
    INTEGER :: S, I, DEG, EDGES(MAXBRANCH), NEIGH(MAXBRANCH), EID
    DOUBLE PRECISION :: LENS(MAXBRANCH)
    
    FORCES = 0D0
    FTENS = 0D0

    ! go through nodes
    IF (VERBOSE) THEN
      PRINT*, 'NODEHIGH', MINNETP%NODEHIGH
    ENDIF

    DO S = 1,MINNETP%NODEHIGH
      IF (MINNETP%NODEACT(S)) THEN

        DEG = MINNETP%NODEDEG(S) ! degree of node
        NEIGH(1:DEG) = MINNETP%NODENODE(S,1:DEG) ! neighboring nodes
        EDGES(1:DEG) = MINNETP%NODEEDGE(S,1:DEG) ! branches from node
        LENS(1:DEG) = MINNETP%EDGELEN(EDGES(1:DEG)) ! lengths of branches from node

        ! loop through neighbors
        DO I = 1,DEG

          IF (LENS(I).EQ.0D0) THEN
            PRINT*, 'sum of lengths less than zero!'
            PRINT*, 'lengths', LENS(1:DEG)
            PRINT*, 'node, degree, neighbors, and edges ', S, DEG, NEIGH, EDGES
          ENDIF
          FTENS(S,:) = FTENS(S,:) + (MINNETP%NODEPOS(S,:)-MINNETP%NODEPOS(NEIGH(I),:))/LENS(I)
        ENDDO
      ENDIF
    ENDDO

    FORCES = -MOB*FTENS;

    
  END SUBROUTINE GETMINNETFORCE


  SUBROUTINE MINNETFROMFILE(MINNETP,NETFILE)
    USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
    USE KEYS, ONLY : MAXBRANCH, MINNETDIM, FIXNODEFROMNETFILE, FIXNODES, NFIX,&
                     & CENTERENCLOSE, RCIRC, CONTINUERUN, NPIN
    USE GENUTIL, ONLY : NORMALIZE
    USE STACKUTIL, ONLY: PUSH
    ! Set up (allocate) a network structure, reading in connectivity and
    ! geometry from an input file    
    
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    CHARACTER(LEN=*), INTENT(IN) :: NETFILE
    LOGICAL :: LDUM, CASESET
    LOGICAL :: FILEEND=.FALSE., DIMSET = .FALSE.
    CHARACTER(LEN=100) :: WORD, LBL, STR, BND
    INTEGER :: NITEMS, NNODE, NEDGE, NODE1, NODE2, DIM, EID, NID, NE
    INTEGER :: LC, I, CC
    DOUBLE PRECISION :: LEN, LENS(MAXBRANCH), LENSAVE(MAXBRANCH), LMAX
    DOUBLE PRECISION, ALLOCATABLE :: DIR(:)
    INTEGER :: TMPARRAY(MAXBRANCH)
    LOGICAL, ALLOCATABLE :: EDGELENSET(:)
    LOGICAL, ALLOCATABLE :: ISFIXED(:)
    INTEGER :: NC, RC, EC, CT
    
    INTEGER, PARAMETER :: NF = 55 ! input file unit number

    ! deallocate any previously set arrays
    IF (MINNETP%ARRAYSET) CALL CLEANUPMINNET(MINNETP)

    ! if detect output is from previous run of code, change this to true to make sure everything
    ! is set up properly!
    CONTINUERUN = .FALSE.
    
    ! go through file and count nodes and branches (total and from each node)
    PRINT*, 'Reading network structure file: ', NETFILE
    INQUIRE(FILE=NETFILE,EXIST=LDUM)
    IF (.NOT.LDUM) THEN
      PRINT*, 'ERROR in MINNETFROMFILE: network file ', TRIM(ADJUSTL(NETFILE)), ' does not exist.'
      STOP 1
    ENDIF
    OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')

    ! Count number of nodes and edges, set spatial dimension
    DIMSET = .FALSE.
    NNODE = 0
    NEDGE = 0

    ! predefine network dimension
    IF (MINNETDIM.GT.0) THEN
       DIMSET = .TRUE.
       DIM = MINNETDIM
    END IF    
    
    DO 
       CALL READLINE(NF,FILEEND,NITEMS)
       
       IF (FILEEND.and.nitems.eq.0) EXIT
       ! skip empty lines
       IF (NITEMS.EQ.0) CYCLE
       ! Read in the keyword for this line
       CALL READA(WORD,CASESET=1)

       IF (WORD.EQ.'NODE') THEN
          NNODE = NNODE+1
          
          ! IF (MINNETDIM.LE.0) DIM = NITEMS-2
          IF (MINNETDIM.LE.0) DIM = 2 ! default choose dimension 2, code can't run with d=3
          IF (NITEMS.GT.10) THEN
            CONTINUERUN = .TRUE.
          ENDIF

       ELSEIF (WORD.EQ.'EDGE') THEN
          NEDGE = NEDGE + 1
       ENDIF
    ENDDO
    CLOSE(NF)

    PRINT*, 'Number of nodes, edges, dimension: ', NNODE, NEDGE, DIM
    ! allocate arrays         
    CALL SETUPMINNET(MINNETP,NNODE,NEDGE,DIM,MAXBRANCH)
        
    ALLOCATE(EDGELENSET(MINNETP%NEDGE), ISFIXED(MINNETP%NNODE))
    ALLOCATE(DIR(MINNETP%DIM))
    EDGELENSET=.FALSE. ! which edges are set directly from file

    PRINT*, 'Reading in node positions and edge connectivity...'
    OPEN(UNIT=NF, FILE=NETFILE, STATUS='OLD')
    ! Get node and edge information
    IF (FIXNODEFROMNETFILE) ISFIXED = .FALSE.
    DO 
      CALL READLINE(NF,FILEEND,NITEMS)
      IF (FILEEND.and.nitems.eq.0) EXIT
      ! skip empty lines
      IF (NITEMS.EQ.0) CYCLE
      ! Read in the keyword for this line
      CALL READA(WORD,CASESET=1)
      ! Skip any comment lines
      IF (WORD(1:1).EQ.'#') CYCLE

      IF (WORD.EQ.'NODE') THEN
        CALL READI(NID)
        IF ((NID.LT.1.OR.NID.GT.NNODE).AND.(.NOT.CONTINUERUN)) THEN
          PRINT*, 'ERROR IN MINNETFROMFILE: invalid node index', NNODE, NID
          STOP 1
        ENDIF
        CALL READF(MINNETP%NODEPOS(NID,1))
        CALL READF(MINNETP%NODEPOS(NID,2))
        MINNETP%NODEACT(NID) = .TRUE.

        IF (NITEMS.EQ.11) THEN ! Specific case of previously run network, import all fields
          CALL READO(MINNETP%NODEFIX(NID))
          CALL READI(MINNETP%NODEGROW(NID))
          CALL READF(MINNETP%NODEDIR(NID,1))
          CALL READF(MINNETP%NODEDIR(NID,2))
          CALL READF(MINNETP%NODEPOSP(NID,1))
          CALL READF(MINNETP%NODEPOSP(NID,2))
          CALL READO(MINNETP%NODEBND(NID))
          ! CONTINUERUN = .TRUE.
        ELSEIF (NITEMS.EQ.12) THEN ! Specific case of previously run network, now with pinning
          CALL READO(MINNETP%NODEFIX(NID))
          CALL READI(MINNETP%NODEGROW(NID))
          CALL READF(MINNETP%NODEDIR(NID,1))
          CALL READF(MINNETP%NODEDIR(NID,2))
          CALL READF(MINNETP%NODEPOSP(NID,1))
          CALL READF(MINNETP%NODEPOSP(NID,2))
          CALL READO(MINNETP%NODEBND(NID))
          CALL READO(MINNETP%NODEPIN(NID))
          ! CONTINUERUN = .TRUE.
        ELSE
          ! loop through remaining items looking for fixed nodes
          DO I = 1,NITEMS-2-MINNETP%DIM
            CALL READA(LBL)
            IF (LBL(1:1).EQ.'F'.AND.FIXNODEFROMNETFILE) THEN ! set this as a fixed node
              ISFIXED(NID) = .TRUE.
            ENDIF
          ENDDO
        ENDIF

        ! IF (MINNETP%DIM.GT.2) THEN
        !   CALL READF(MINNETP%NODEPOS(NID,3))!this line had to be added for 3d
        ! ENDIF

      ELSEIF (WORD.EQ.'EDGE') THEN
        CALL READI(EID) ! edge id
        CALL READI(NODE1)
        CALL READI(NODE2)
        MINNETP%EDGEACT(EID) = .TRUE.

        IF (NITEMS.EQ.8) THEN ! previously run network, import all information!
          CALL READF(MINNETP%EDGELEN(EID))
          EDGELENSET(EID) = .TRUE.
          CALL READO(MINNETP%EDGEBND(EID))
          CALL READI(MINNETP%EDGENODEP(EID,1))
          CALL READI(MINNETP%EDGENODEP(EID,2))
          ! CONTINUERUN = .TRUE.
        ELSEIF (NITEMS.GT.4) THEN ! read in edge length directly
          CALL READF(MINNETP%EDGELEN(EID))
          EDGELENSET(EID) = .TRUE.
        ENDIF
        
        IF (CONTINUERUN) THEN 
          IF (NODE1.LT.1.OR.NODE2.LT.1.OR.EID.LT.1) THEN
            PRINT*, 'ERROR IN MINNET FROM FILE: &
                 & bad edge or node indices while reading edge.', &
                 ' EID, NODE1, NODE2, NNODE, NEDGE:', &
                 & EID, NODE1, NODE2, NNODE, NEDGE
            STOP 1
          ENDIF
        ! not continuuing run, check for correct node/edge indices
        ELSE 
          IF (NODE1.LT.1.OR.NODE2.LT.1.OR.EID.LT.1&
               & .OR.NODE1.GT.NNODE.OR.NODE2.GT.NNODE.OR.EID.GT.NEDGE) THEN
            PRINT*, 'ERROR IN MINNET FROM FILE: &
                 & bad edge or node indices while reading edge.', &
                 ' EID, NODE1, NODE2, NNODE, NEDGE:', &
                 & EID, NODE1, NODE2, NNODE, NEDGE
            STOP 1
          ENDIF
        ENDIF

        ! nodes connected to this edge           
        MINNETP%EDGENODE(EID,:) = (/NODE1,NODE2/)
      ENDIF        
    END DO
    CLOSE(NF)

    IF (FIXNODEFROMNETFILE) THEN
       ! update fixed nodes from network file
       
      NFIX = 0
      FIXNODES = 0
      DO NC = 1,NNODE
        IF (ISFIXED(NC)) THEN
          NFIX = NFIX + 1
          FIXNODES(NFIX) = NC
        ENDIF
      ENDDO
      PRINT*, 'Fixed nodes for field 1 only:', NFIX, FIXNODES(1:NFIX)
    END IF
    
    PRINT*, 'Setting up connectivity and edge lengths...'

    ! for finding the longest edge
    LMAX = 0D0

    ! old way 
    ! DO EID = 1,NEDGE
    DO EID = 1,MINNETP%NEDGE
      IF (MINNETP%EDGEACT(EID)) THEN

        NODE1 = MINNETP%EDGENODE(EID,1)
        NODE2 = MINNETP%EDGENODE(EID,2)
 
        ! edge start, length, and (normalized) direction
        MINNETP%EDGESTART(EID,:) = MINNETP%NODEPOS(NODE1,:)
        DIR(1:MINNETP%DIM) = MINNETP%NODEPOS(NODE2,:)-MINNETP%NODEPOS(NODE1,:)   

        IF (.NOT.EDGELENSET(EID))  THEN ! edge length was not set in network file
          MINNETP%EDGELEN(EID) = SQRT(SUM( DIR*DIR ))
        ENDIF
        MINNETP%EDGEDIR(EID,:) = DIR/MINNETP%EDGELEN(EID)
 
        ! increment degrees of the nodes
        MINNETP%NODEDEG(NODE1) = MINNETP%NODEDEG(NODE1)+1
        MINNETP%NODEDEG(NODE2) = MINNETP%NODEDEG(NODE2)+1
 
        IF ( MAX(MINNETP%NODEDEG(NODE1),MINNETP%NODEDEG(NODE2)).GT.MAXBRANCH) THEN
          PRINT*, 'ERROR IN MINNETFROMFILE: node degree exceeds maximum.',&
               & NODE1, NODE2, MAXBRANCH, MINNETP%NODEDEG(NODE1), MINNETP%NODEDEG(NODE2)
          STOP 1
        ENDIF
 
        ! keep track of longest edge
        ! check if longer than LMAX
        IF (MINNETP%EDGELEN(EID).GT.LMAX) THEN
          LMAX = MINNETP%EDGELEN(EID)
          MINNETP%EIDMAX = EID
        ENDIF
 
        ! edges connected to each node
        MINNETP%NODEEDGE(NODE1,MINNETP%NODEDEG(NODE1)) = EID
        MINNETP%NODEEDGE(NODE2,MINNETP%NODEDEG(NODE2)) = EID
 
        ! nodes connected to each node
        MINNETP%NODENODE(NODE1,MINNETP%NODEDEG(NODE1)) = NODE2
        MINNETP%NODENODE(NODE2,MINNETP%NODEDEG(NODE2)) = NODE1              
      ENDIF
    ENDDO

    ! NODEON initially zero, in this loop we add any non-fixed nodes from the initial network
    MINNETP%NODETRACK=0
 
    ! for each node, sort edges by length
    CT=0
    DO NID = 1,MINNETP%NNODE
      IF (MINNETP%NODEACT(NID)) THEN
        NE = MINNETP%NODEDEG(NID)
        LENS(1:NE) = MINNETP%EDGELEN(MINNETP%NODEEDGE(NID,1:NE))
        LENSAVE = LENS
        TMPARRAY = MINNETP%NODEEDGE(NID,:)        
        CALL SORT2INT(NE,LENS,TMPARRAY)
        MINNETP%NODEEDGE(NID,1:NE) = TMPARRAY(1:NE)

        TMPARRAY = MINNETP%NODENODE(NID,:)
        CALL SORT2INT(NE,LENSAVE,TMPARRAY)
        MINNETP%NODENODE(NID,1:NE) = TMPARRAY(1:NE)

        ! if not a permanently fixed boundary node, add to list of tracked nodes
        ! allow tracking of sliding boundary nodes, but don't track the fixed boundary nodes
        IF (.NOT.MINNETP%NODEFIX(NID)) THEN
          MINNETP%NODEON = MINNETP%NODEON+1
          MINNETP%NODETRACK(MINNETP%NODEON) = NID
        ENDIF

      ENDIF
    ENDDO


    IF (CONTINUERUN) THEN
      ! set node and edge high
      MINNETP%NODEHIGH = FINDLOC(MINNETP%NODEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)
      MINNETP%EDGEHIGH = FINDLOC(MINNETP%EDGEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)

      ! Fixed nodes are already set up?

      ! active/inactive setup
      MINNETP%NACTIVE = NNODE
      MINNETP%EACTIVE = NEDGE

      ! push unused nodes/edges to the stacks!
      ! count down from highest node
      DO NC=1,MINNETP%NNODE
        NID = MINNETP%NNODE - NC + 1
        IF (.NOT.MINNETP%NODEACT(NID)) THEN
          CALL PUSH(MINNETP%NODESTACK,NID)
        ENDIF
      ENDDO

      DO EC=1,MINNETP%NEDGE
        EID = MINNETP%NEDGE - EC + 1
        IF (.NOT.MINNETP%EDGEACT(EID)) THEN
          CALL PUSH(MINNETP%EDGESTACK,EID)
        ENDIF
      ENDDO

      ! set pinned nodes only if NPIN is set by user
      ! NPIN being set means either RANDPINNODES was provided
      ! or they gave a list of PINNODES to set manually
      IF (NPIN.GT.0) CALL SETUPPINNODES(MINNETP)

      ! previous arrays already set up
      CENTERENCLOSE = .TRUE.

    ELSE

      ! set node and edge high
      MINNETP%NODEHIGH = NNODE
      MINNETP%EDGEHIGH = NEDGE

      ! set fixed nodes
      CALL SETUPFIXNODES(MINNETP)
      
      ! set pinned nodes
      CALL SETUPPINNODES(MINNETP)

      ! set pinned nodes

      ! set up active/inactive arrays/stacks
      ! first NNODE nodes are initially active
      MINNETP%NODEACT(1:NNODE) = .TRUE.
      MINNETP%EDGEACT(1:NEDGE) = .TRUE.
      MINNETP%NACTIVE = NNODE
      MINNETP%EACTIVE = NEDGE

      ! set up previous arrays
      MINNETP%NODEPOSP(:,:) = MINNETP%NODEPOS(:,:)
      MINNETP%EDGENODEP(:,:) = MINNETP%EDGENODE(:,:)
    ENDIF

    DEALLOCATE(EDGELENSET,ISFIXED)
    DEALLOCATE(DIR)

    MINNETP%STRUCTURESET = .TRUE.

  END SUBROUTINE MINNETFROMFILE


  SUBROUTINE CENTENC(MINNETP)
    USE STACKUTIL, ONLY : POP
    USE KEYS, ONLY : RMULT, RCIRC, CONTINUERUN, NBNDEDGES, RNODES
    USE GENUTIL, ONLY : PI, LLINTERSECT
    ! center and enclose the network in a circle of stationary nodes all connected to eachother
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    DOUBLE PRECISION :: MEANX, MEANY, RS(MINNETP%NNODE), RMAX, RPOS
    DOUBLE PRECISION :: X, Y, T, U, DIR(2), NEWRCIRC
    INTEGER :: NEDGE, EC, EID, NID, OLDNID, NODE1, N1, N2, NC
    LOGICAL :: CROSS

    ! if continuing a previous run, simply find the minimum radial distance to a boundary edge
    ! and that sets RCIRC. Nothing else is needed
    IF (CONTINUERUN) THEN
      RS = SQRT(SUM(MINNETP%NODEPOS**2,DIM=2))
      RMAX = MAXVAL(RS,MASK=MINNETP%NODEACT)
      RCIRC=1D10
      DO EC=1,MINNETP%EDGEHIGH
        IF (MINNETP%EDGEACT(EC).AND.MINNETP%EDGEBND(EC)) THEN
          ! use edgedir? find perp comp and then see if perp component hits edge segment if starting at origin?
          ! if not then use the endpoint with the smallest radial distance

          ! DIR is perpendicular to edgedir
          DIR(1) = -MINNETP%EDGEDIR(EC,2)*RMAX
          DIR(2) = MINNETP%EDGEDIR(EC,1)*RMAX
          X = 0D0
          Y = 0D0
          CROSS = .FALSE.
          T=0D0
          U=0D0
          N1 = MINNETP%EDGENODE(EC,1)
          N2 = MINNETP%EDGENODE(EC,2)

          CALL LLINTERSECT(X,Y,DIR(1),DIR(2),MINNETP%NODEPOS(N1,1),MINNETP%NODEPOS(N1,2),&
               & MINNETP%NODEPOS(N2,1),MINNETP%NODEPOS(N2,2),CROSS,T,U)

          IF (CROSS) THEN
            ! crossed first direction attempt
            ! set RCIRC
            NEWRCIRC = T*RMAX
            IF (NEWRCIRC.LT.RCIRC) THEN
              RCIRC = NEWRCIRC
            ENDIF

          ELSE
            ! try negating the vector
            CALL LLINTERSECT(X,Y,-DIR(1),-DIR(2),MINNETP%NODEPOS(N1,1),MINNETP%NODEPOS(N1,2),&
                 & MINNETP%NODEPOS(N2,1),MINNETP%NODEPOS(N2,2),CROSS,T,U)
            IF (CROSS) THEN
              ! crossed first direction attempt
              ! set RCIRC
              NEWRCIRC = T*RMAX
              IF (NEWRCIRC.LT.RCIRC) THEN
                RCIRC = NEWRCIRC
              ENDIF

            ENDIF
          ENDIF
        ENDIF
      ENDDO

      ! find RNODES as well
      ! this is the radial position of the fixed boundary nodes, which is larger than rcirc
      DO NC=1,MINNETP%NODEHIGH
        IF (MINNETP%NODEACT(NC).AND.MINNETP%NODEBND(NC).AND.MINNETP%NODEFIX(NC)) THEN
          RPOS = SQRT(SUM(MINNETP%NODEPOS(NC,:)**2))
          IF (VERBOSE) THEN
            PRINT*, 'rposition of active, fixed boundary node: ', RPOS
          ENDIF
          RNODES=RPOS
          EXIT
        ENDIF
      ENDDO

      ! having looped through edges and found the closest boundary edge to center, return
      RETURN

    ENDIF

    ! centering
    MEANX = SUM(MINNETP%NODEPOS(:,1),DIM=1,MASK=MINNETP%NODEACT)/MINNETP%NACTIVE
    MEANY = SUM(MINNETP%NODEPOS(:,2),DIM=1,MASK=MINNETP%NODEACT)/MINNETP%NACTIVE

    MINNETP%NODEPOS(:,1) = MINNETP%NODEPOS(:,1) - MEANX
    MINNETP%NODEPOS(:,2) = MINNETP%NODEPOS(:,2) - MEANY
    MEANX = SUM(MINNETP%NODEPOS(:,1),DIM=1,MASK=MINNETP%NODEACT)/MINNETP%NACTIVE
    MEANY = SUM(MINNETP%NODEPOS(:,2),DIM=1,MASK=MINNETP%NODEACT)/MINNETP%NACTIVE

    ! enclosing
    ! find max radial distance of all nodes
    RS = SQRT(SUM(MINNETP%NODEPOS**2,DIM=2))
    RMAX = MAXVAL(RS,MASK=MINNETP%NODEACT)

    RCIRC = RMULT*RMAX
    ! if user input number of boundary edges use that
    IF (NBNDEDGES.GT.0) THEN
      NEDGE = NBNDEDGES
    ELSE
    ! else create boundary circle with edges of edge length at most (max length/2)
      NEDGE = CEILING(4*PI*RCIRC/MINNETP%EDGELEN(MINNETP%EIDMAX))
    ENDIF

    ! add first node
    CALL POP(MINNETP%NODESTACK,NID)
    CALL ADDNODE(MINNETP,RCIRC,0D0,NID)
    OLDNID = NID

    ! save first node of circle to connect to final node
    NODE1 = MINNETP%NACTIVE

    DO EC=1,NEDGE-1
      X = RCIRC*COS(2*PI*EC/NEDGE)
      Y = RCIRC*SIN(2*PI*EC/NEDGE)

      CALL POP(MINNETP%NODESTACK,NID)
      CALL POP(MINNETP%EDGESTACK,EID)
      ! connect it to last added node (just nactive)
      CALL ADDNODE(MINNETP,X,Y,NID)
      CALL CONNECTNODE(MINNETP,NID,OLDNID,EID)
      
      ! set this as boundary edge
      MINNETP%EDGEBND(EID) = .TRUE.

      ! update oldnid
      OLDNID = NID
    ENDDO

    ! connect node 1 to last node
    CALL POP(MINNETP%EDGESTACK,EID)
    CALL CONNECTNODE(MINNETP,NID,NODE1,EID)
    ! set this as boundary edge
    MINNETP%EDGEBND(EID) = .TRUE.

    ! update RCIRC to be minimal distance to boundary edges which is smaller than original RCIRC
    ! due to approximating circle as polygon
    RNODES=RCIRC
    RCIRC = RCIRC*COS(PI/NEDGE)


  END SUBROUTINE CENTENC



  SUBROUTINE CONNECTNODE(MINNETP,NID,NIDIN,EID)
    ! connect node NID with node NIDIN, use edge EID
    USE STACKUTIL, ONLY : POP
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NIDIN, NID, EID

    ! update node NID
    MINNETP%NODEDEG(NID) = MINNETP%NODEDEG(NID) + 1
    MINNETP%NODENODE(NID,MINNETP%NODEDEG(NID)) = NIDIN
    MINNETP%NODEEDGE(NID,MINNETP%NODEDEG(NID)) = EID

    ! update new edge, EID
    MINNETP%EDGENODE(EID,1:2) = (/NIDIN,NID/)
    MINNETP%EDGESTART(EID,:) = MINNETP%NODEPOS(NIDIN,:)
    MINNETP%EDGEDIR(EID,:) = MINNETP%NODEPOS(NID,:) - MINNETP%NODEPOS(NIDIN,:)
    MINNETP%EDGELEN(EID) = SQRT(SUM(MINNETP%EDGEDIR(EID,:)**2))
    MINNETP%EDGEDIR(EID,:) = MINNETP%EDGEDIR(EID,:)/MINNETP%EDGELEN(EID)
    MINNETP%EDGEACT(EID) = .TRUE.
    MINNETP%EACTIVE = MINNETP%EACTIVE + 1
    IF (EID.GT.MINNETP%EDGEHIGH) THEN
      MINNETP%EDGEHIGH = EID
    ENDIF
    ! MINNETP%EDGEHIGH = FINDLOC(MINNETP%EDGEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)
    MINNETP%EDGEBND(EID) = .FALSE.

    ! update already existing node, NIDIN
    MINNETP%NODEDEG(NIDIN) = MINNETP%NODEDEG(NIDIN) + 1
    MINNETP%NODENODE(NIDIN,MINNETP%NODEDEG(NIDIN)) = NID
    MINNETP%NODEEDGE(NIDIN,MINNETP%NODEDEG(NIDIN)) = EID

    IF (VERBOSE) THEN
      PRINT*, 'connected node', NID, 'to node', NIDIN, 'with edge', EID
    ENDIF


  END SUBROUTINE CONNECTNODE



  SUBROUTINE ADDNODE(MINNETP,X,Y,NID)
    ! add node NID at position x,y
    USE STACKUTIL, ONLY : POP
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    DOUBLE PRECISION, INTENT(IN) :: X, Y
    INTEGER, INTENT(IN) :: NID

    MINNETP%NODEPOS(NID,:) = (/X,Y/)
    MINNETP%NODENODE(NID,1) = 0
    MINNETP%NODEDEG(NID) = 0
    MINNETP%NODEEDGE(NID,1) = 0
    MINNETP%NODEFIX(NID) = .TRUE.
    MINNETP%NODEACT(NID) = .TRUE.
    MINNETP%NACTIVE = MINNETP%NACTIVE + 1
    MINNETP%NODEGROW(NID) = 0
    MINNETP%NODEDIR(NID,:) = 0D0
    IF (NID.GT.MINNETP%NODEHIGH) THEN
      MINNETP%NODEHIGH = NID
    ENDIF
    ! MINNETP%NODEHIGH = FINDLOC(MINNETP%NODEACT,DIM=1,VALUE=.TRUE.,BACK=.TRUE.)
    MINNETP%NODEBND(NID) = .TRUE.

    IF (VERBOSE) THEN
      PRINT*, 'Added node, ', NID
    ENDIF


  END SUBROUTINE ADDNODE



  SUBROUTINE OUTPUTMINNET(MINNETP,NETOUTFILE,APPEND)
    ! output a network file describing network structure
    ! if NODEVAL is present: write an additional value for each node
    ! if EDGEVAL is present: write an additional value for each edge
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    CHARACTER (LEN=*), INTENT(IN) :: NETOUTFILE
    LOGICAL, INTENT(IN) :: APPEND
!    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: NODEVAL(MINNETP%NNODE), EDGEVAL(MINNETP%NEDGE)
    INTEGER, PARAMETER :: FU=51
    INTEGER :: NC, EC, LC
    CHARACTER(LEN=100) :: NODEFMTSTR, NODEFMTSTRRESV, FMTN, FMTE
    CHARACTER(LEN=6) :: NUMSTRING

    ! maximum number of nodes/edges is 999,999 (only exporting 6 positions in node number/edge number)

    WRITE(FMTN,'(A,I1,A,I1,A,I1,A)') &
    & '(A,1X,I6,1X,',MINNETP%DIM,'F20.10,1X,L1,1X,I1,1X,',MINNETP%DIM,'F20.10,1X,',MINNETP%DIM,'F20.10,1X,L1,1X,L1)'
    
    WRITE(FMTE,'(A)') '(A,1X,I6,1X,I6,1X,I6,1X,F20.10,1X,L1,1X,I6,1X,I6)'

    IF (APPEND) THEN
      OPEN(UNIT=FU,FILE=NETOUTFILE,STATUS='UNKNOWN',ACCESS='APPEND')
    ELSE
      OPEN(UNIT=FU,FILE=NETOUTFILE,STATUS='UNKNOWN')
    ENDIF
    WRITE(FU,'(A)')'# network structure output by fortran code'
    
    DO NC = 1,MINNETP%NODEHIGH
      IF (MINNETP%NODEACT(NC)) THEN
        WRITE(FU,FMTN) 'NODE ', NC, MINNETP%NODEPOS(NC,:), MINNETP%NODEFIX(NC), MINNETP%NODEGROW(NC),&
                      & MINNETP%NODEDIR(NC,:), MINNETP%NODEPOSP(NC,:), MINNETP%NODEBND(NC), MINNETP%NODEPIN(NC)
      ENDIF
    ENDDO

    DO EC = 1,MINNETP%EDGEHIGH
      IF (MINNETP%EDGEACT(EC)) THEN
        WRITE(FU,FMTE) 'EDGE ', EC, MINNETP%EDGENODE(EC,:), MINNETP%EDGELEN(EC), MINNETP%EDGEBND(EC),&
                      & MINNETP%EDGENODEP(EC,:)
        IF (MINNETP%EDGENODE(EC,1).EQ.MINNETP%EDGENODE(EC,2)) THEN
          PRINT*, 'ERROR IN OUTPUTMINNET: loop edge found.'
          STOP 1
        ENDIF
      ENDIF
    ENDDO  
    
    WRITE(FU,*) ''
    CLOSE(FU)
    PRINT*, 'State of Random Seed'
    
  END SUBROUTINE OUTPUTMINNET


  SUBROUTINE OUTPUTSNAPSHOT(MINNETP,FILENAME,APPEND)
    ! dump out a snapshot of the node positions followed by edgenodes
    ! if append=true, append to an existing file
    TYPE(MINNET), POINTER :: MINNETP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL, INTENT(IN) :: APPEND
    CHARACTER*100 :: FMT, FMTN
    INTEGER :: B
 
    IF (APPEND) THEN
      OPEN(UNIT=99,FILE=FILENAME,POSITION='APPEND')
    ELSE
      OPEN(UNIT=99,FILE=FILENAME,POSITION='REWIND')
    ENDIF
 
    ! write information line
    WRITE(99,'(A,1X,I12)') 'N',MINNETP%NACTIVE
    WRITE(99,'(A,1X,I12)') 'E',MINNETP%EACTIVE
 
    FMT = '(A,1X,12G20.10)'
    FMTN = '(A,1X,I12,12G20.10)'
 
    ! write out positions of network nodes
    DO B = 1,MINNETP%NODEHIGH
      ! only print active nodes
      IF (MINNETP%NODEACT(B)) THEN
        ! print node number and then nodepos, since not all node numbers will be filled
        WRITE(99,FMTN) 'P', B, MINNETP%NODEPOS(B,:)
        ! WRITE(99,FMTN) 'P', B, MINNETP%NODEPOS(B,:), MINNETP%NODEPIN(B)
      ENDIF
    ENDDO

    DO B = 1,MINNETP%EDGEHIGH
      ! only print active edges
      IF (MINNETP%EDGEACT(B)) THEN
        WRITE(99,FMT) 'E', MINNETP%EDGENODE(B,:)
      ENDIF
    ENDDO
  END SUBROUTINE OUTPUTSNAPSHOT


  SUBROUTINE OUTPUTPIN(MINNETP,FILENAME,APPEND)
    ! dump pinned node info
    ! if append=true, append to an existing file
    TYPE(MINNET), POINTER :: MINNETP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL, INTENT(IN) :: APPEND
    CHARACTER*100 :: FMT
    INTEGER :: B, PINNEDNODES(MINNETP%NNODE), NPIN
    
 
    IF (APPEND) THEN
      OPEN(UNIT=97,FILE=FILENAME,POSITION='APPEND')
    ELSE
      OPEN(UNIT=97,FILE=FILENAME,POSITION='REWIND')
    ENDIF

    NPIN=0
    ! count pinned nodes...
    DO B = 1,MINNETP%NNODE
      IF (MINNETP%NODEPIN(B)) THEN
        NPIN=NPIN+1
        PINNEDNODES(NPIN) = B
      ENDIF
    ENDDO
 
    ! write information line
    WRITE(97,'(A,1X,I12)') 'N',NPIN
    FMT = '(I12,1X,I12)'
 
    ! loop over all nodes, output pinned ones
    ! write out positions of network nodes
    IF (NPIN.GT.0) THEN
      DO B = 1,NPIN
        WRITE(97,FMT) B, PINNEDNODES(B)
      ENDDO
    ENDIF

  END SUBROUTINE OUTPUTPIN


  SUBROUTINE OUTPUTBOUNDS(MINNETP,FILENAME,APPEND)
    ! dump boundary node info
    ! if append=true, append to an existing file
    TYPE(MINNET), POINTER :: MINNETP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL, INTENT(IN) :: APPEND
    CHARACTER*100 :: FMT
    INTEGER :: B, BNDNODES(MINNETP%NNODE), NBND
    
 
    IF (APPEND) THEN
      OPEN(UNIT=96,FILE=FILENAME,POSITION='APPEND')
    ELSE
      OPEN(UNIT=96,FILE=FILENAME,POSITION='REWIND')
    ENDIF

    NBND=0
    ! count pinned nodes...
    DO B = 1,MINNETP%NNODE
      IF (MINNETP%NODEBND(B)) THEN
        NBND=NBND+1
        BNDNODES(NBND) = B
      ENDIF
    ENDDO
 
    ! write information line
    WRITE(96,'(A,1X,I12)') 'N',NBND
    FMT = '(I12,1X,I12)'
 
    ! loop over all nodes, output pinned ones
    ! write out positions of network nodes
    IF (NBND.GT.0) THEN
      DO B = 1,NBND
        WRITE(96,FMT) B, BNDNODES(B)
      ENDDO
    ENDIF

  END SUBROUTINE OUTPUTBOUNDS



  SUBROUTINE OUTPUTEBOUNDS(MINNETP,FILENAME,APPEND)
    ! dump boundary node info
    ! if append=true, append to an existing file
    TYPE(MINNET), POINTER :: MINNETP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL, INTENT(IN) :: APPEND
    CHARACTER*100 :: FMT
    INTEGER :: B, BNDEDGES(MINNETP%NNODE), NBND, NINACTEDGES(MINNETP%NEDGE)
    
 
    IF (APPEND) THEN
      OPEN(UNIT=95,FILE=FILENAME,POSITION='APPEND')
    ELSE
      OPEN(UNIT=95,FILE=FILENAME,POSITION='REWIND')
    ENDIF

    NBND=0
    ! have to do convoluted way of outputting the edge information so it is compatiable
    ! with the plotting in matlab. So instead of outputting edge number as in Fortran, 
    ! output (edge number - number of inactive edges preceding that edge)... bad but 
    ! should work
    ! count boundary edges
    NINACTEDGES=0

    IF (.NOT.MINNETP%EDGEACT(1)) THEN
      NINACTEDGES(1) = 1
    ENDIF
    IF (MINNETP%EDGEBND(1)) THEN
      NBND=NBND+1
      BNDEDGES(NBND) = B
    ENDIF

    DO B = 2,MINNETP%NEDGE
      IF (.NOT.MINNETP%EDGEACT(B)) THEN
        NINACTEDGES(B) = NINACTEDGES(B-1) + 1
      ELSE
        NINACTEDGES(B) = NINACTEDGES(B-1)
      ENDIF
      IF (MINNETP%EDGEBND(B)) THEN
        NBND=NBND+1
        BNDEDGES(NBND) = B - NINACTEDGES(B)
      ENDIF
    ENDDO
 
    ! write information line
    WRITE(95,'(A,1X,I12)') 'N',NBND
    FMT = '(I12,1X,I12)'
 
    ! loop over all edges, output boundary edges
    IF (NBND.GT.0) THEN
      DO B = 1,NBND
        WRITE(95,FMT) B, BNDEDGES(B)
      ENDDO
    ENDIF

  END SUBROUTINE OUTPUTEBOUNDS

 
  
  SUBROUTINE OUTPUTTRACK(MINNETP,FILENAME,APPEND)
    ! dump track node info
    ! if append=true, append to an existing file
    TYPE(MINNET), POINTER :: MINNETP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL, INTENT(IN) :: APPEND
    CHARACTER*100 :: FMT
    INTEGER :: B
 
    IF (APPEND) THEN
      OPEN(UNIT=98,FILE=FILENAME,POSITION='APPEND')
    ELSE
      OPEN(UNIT=98,FILE=FILENAME,POSITION='REWIND')
    ENDIF
 
    ! write information line
    WRITE(98,'(A,1X,I12)') 'N',MINNETP%NODEON
    FMT = '(I12,1X,I12)'
 
    IF (MINNETP%NODEON.GT.0) THEN
      ! write out positions of network nodes
      DO B = 1,MINNETP%NODEON
        WRITE(98,FMT) B, MINNETP%NODETRACK(B)
      ENDDO
    ENDIF

  END SUBROUTINE OUTPUTTRACK
 
  
  SUBROUTINE SETUPMINNET(MINNETP,NNODE,NEDGE,DIM,MAXBRANCH)
    ! set up a network by allocating arrays
    ! MINNETP: Pointer to a network
    ! NNODE: number of nodes
    ! NEDGE: number of branches
    ! DIM: spatial dimensionality where network resides
    ! MAXBRANCH: maximum branches attached to each node
    ! CONTINUERUN: if true, then continuing a previous run, must be careful with stacks
    
    USE KEYS, ONLY : EXTRANODEFACT, EXTRAEDGEFACT, CONTINUERUN, EXTRATRACKFACT
    USE STACKUTIL, ONLY : INITIALIZESTACK
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP
    INTEGER, INTENT(IN) :: NNODE, NEDGE, DIM, MAXBRANCH
    INTEGER :: LN((EXTRANODEFACT-1)*NNODE), LE((EXTRAEDGEFACT-1)*NEDGE), I

    ! all network variables initialized with space for extra nodes/edges
    MINNETP%NNODE = NNODE*EXTRANODEFACT
    MINNETP%NEDGE = NEDGE*EXTRAEDGEFACT
    MINNETP%DIM = DIM
    MINNETP%NODEON=0
    
    ! allocate node data
    ALLOCATE(MINNETP%NODENODE(NNODE*EXTRANODEFACT,MAXBRANCH), MINNETP%NODEEDGE(NNODE*EXTRANODEFACT,MAXBRANCH))
    ALLOCATE(MINNETP%NODEPOS(NNODE*EXTRANODEFACT,DIM), MINNETP%NODEDEG(NNODE*EXTRANODEFACT),&
         & MINNETP%NODEACT(NNODE*EXTRANODEFACT), MINNETP%NODEGROW(NNODE*EXTRANODEFACT),&
         & MINNETP%NODEDIR(NNODE*EXTRANODEFACT,DIM), MINNETP%NODEPOSP(NNODE*EXTRANODEFACT,DIM),&
         & MINNETP%NODEBND(NNODE*EXTRANODEFACT),MINNETP%NODETRACK(NNODE*EXTRATRACKFACT),&
         & MINNETP%NODEPIN(NNODE*EXTRANODEFACT))
    
    ! allocate branch data
    ALLOCATE(MINNETP%EDGENODE(NEDGE*EXTRAEDGEFACT,2), MINNETP%EDGESTART(NEDGE*EXTRAEDGEFACT,DIM),&
         & MINNETP%EDGEDIR(NEDGE*EXTRAEDGEFACT,DIM), MINNETP%EDGELEN(NEDGE*EXTRAEDGEFACT),&
         & MINNETP%EDGEACT(EXTRAEDGEFACT*NEDGE),MINNETP%EDGEBND(EXTRAEDGEFACT*NEDGE),&
         & MINNETP%EDGENODEP(NEDGE*EXTRAEDGEFACT,2))
 
    MINNETP%ARRAYSET = .TRUE.
    MINNETP%NODEDEG = 0
    MINNETP%NODEEDGE = 0; MINNETP%EDGENODE = 0; MINNETP%NODENODE = 0
    MINNETP%NODEACT = .FALSE.
    MINNETP%NODEBND = .FALSE.
    MINNETP%EDGEACT = .FALSE.
    MINNETP%EDGEBND = .FALSE.
    MINNETP%NODEPIN = .FALSE.
    MINNETP%NODEGROW = 0

    ! allocate the pointers for stacks
    ALLOCATE(MINNETP%NODESTACK, MINNETP%EDGESTACK)

    IF (CONTINUERUN) THEN
      ! set up the stacks, no prefilling, will add unused edges/nodes in MINNETFROMFILE
      ! once we know which ones aren't being used
      CALL INITIALIZESTACK(MINNETP%NODESTACK,NNODE*EXTRANODEFACT)
      CALL INITIALIZESTACK(MINNETP%EDGESTACK,NEDGE*EXTRAEDGEFACT)

    ELSE
      ! set up the stacks, prefilled with unused nodes/edges 
      LN = [(NNODE*EXTRANODEFACT + NNODE + 1 - I, I=(NNODE+1), NNODE*EXTRANODEFACT)]
      CALL INITIALIZESTACK(MINNETP%NODESTACK,NNODE*EXTRANODEFACT,LN)
      LE = [(NEDGE*EXTRAEDGEFACT + NEDGE + 1 - I, I=(NEDGE+1), NEDGE*EXTRAEDGEFACT)]
      CALL INITIALIZESTACK(MINNETP%EDGESTACK,NEDGE*EXTRAEDGEFACT,LE)
    ENDIF
    
    ! fixed nodes
    ALLOCATE(MINNETP%NODEFIX(MINNETP%NNODE*EXTRANODEFACT))

    ! growth nodes
    ALLOCATE(MINNETP%NODEGROW(MINNETP%NNODE*EXTRANODEFACT), MINNETP%NODEDIR(MINNETP%NNODE*EXTRANODEFACT,DIM))

  END SUBROUTINE SETUPMINNET



  SUBROUTINE CLEANUPMINNET(MINNETP)
    ! deallocate arrays for the network structure
    USE STACKUTIL, ONLY : DEALLOCATESTACK
    IMPLICIT NONE
    TYPE(MINNET), POINTER :: MINNETP

    
    IF (MINNETP%ARRAYSET) THEN
       DEALLOCATE(MINNETP%NODENODE, MINNETP%NODEEDGE, MINNETP%NODEPOS, &
            & MINNETP%NODEDEG, MINNETP%NODEFIX, MINNETP%NODEACT, &
            & MINNETP%NODEGROW, MINNETP%NODEDIR,MINNETP%NODEBND, &
            & MINNETP%NODEPIN)
       DEALLOCATE(MINNETP%EDGENODE, MINNETP%EDGESTART, MINNETP%EDGEDIR, &
            & MINNETP%EDGELEN, MINNETP%EDGEACT)

       ! deallocate stack variables
       CALL DEALLOCATESTACK(MINNETP%NODESTACK)
       CALL DEALLOCATESTACK(MINNETP%EDGESTACK)
       DEALLOCATE(MINNETP%NODESTACK, MINNETP%EDGESTACK)
    ENDIF
    
    MINNETP%ARRAYSET = .FALSE.
    MINNETP%STRUCTURESET = .FALSE.
  END SUBROUTINE CLEANUPMINNET



END MODULE MINNETUTIL
