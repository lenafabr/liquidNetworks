MODULE STACKUTIL
  ! Utilities for creating a stack, popping and pushing values

  IMPLICIT NONE

  TYPE STACK

    ! maximum size of stack
    INTEGER :: MAXSIZE 
    ! current size of stack
    INTEGER :: TOP
    ! current items
    INTEGER, POINTER :: ITEMS(:)
    ! stack is allocated?
    LOGICAL :: SET = .FALSE.
    
  END TYPE STACK

CONTAINS



  SUBROUTINE INITIALIZESTACK(STACKP,MAXSIZE,ITEMS)
    ! set up new stack, with specified maxsize
    ! optional inputs of INSIZE and ITEMS to initialize with a large initial input
    ! INSIZE specifies size of ITEMS.
    ! ITEMS is an optional initial array of integers to fill stack with

    IMPLICIT NONE
    TYPE(STACK), POINTER :: STACKP
    INTEGER, INTENT(IN) :: MAXSIZE
    INTEGER, INTENT(IN), OPTIONAL :: ITEMS(:)
    INTEGER :: INSIZE
    
    ! deallocate any previously set stack
    IF (STACKP%SET) THEN
      DEALLOCATE(STACKP%ITEMS)
      STACKP%SET = .FALSE.
    ENDIF

    ! allocate array for stack
    ALLOCATE(STACKP%ITEMS(MAXSIZE))
    STACKP%SET = .TRUE.
    STACKP%MAXSIZE = MAXSIZE

    ! using input ITEMS, set for INSIZE elements of stack
    IF (PRESENT(ITEMS)) THEN
      INSIZE = SIZE(ITEMS)
      STACKP%TOP = INSIZE
      STACKP%ITEMS(1:INSIZE) = ITEMS
    ! no optional arguments provided, initialize empty stack
    ELSE
      ! set up empty stack
      STACKP%TOP = 0
    ENDIF

  END SUBROUTINE INITIALIZESTACK


  SUBROUTINE PUSH(STACKP,ITEM)
    IMPLICIT NONE
    TYPE(STACK), POINTER :: STACKP
    INTEGER, INTENT(IN) :: ITEM

    IF (STACKP%TOP.EQ.STACKP%MAXSIZE) THEN
      PRINT*, 'ERROR, stack overflow'
      STOP 1
    ELSE
      STACKP%TOP = STACKP%TOP + 1
      STACKP%ITEMS(STACKP%TOP) = ITEM
    ENDIF

  END SUBROUTINE PUSH


  SUBROUTINE POP(STACKP,ITEM)
    IMPLICIT NONE
    TYPE(STACK), POINTER :: STACKP
    INTEGER, INTENT(OUT) :: ITEM

    IF (STACKP%TOP.EQ.0) THEN
      PRINT*, 'ERROR, stack underflow'
      STOP 1
    ELSE
      ITEM = STACKP%ITEMS(STACKP%TOP)
      STACKP%TOP = STACKP%TOP - 1
    ENDIF

  END SUBROUTINE POP
  

  SUBROUTINE DEALLOCATESTACK(STACKP)
    ! deallocate stack object
    IMPLICIT NONE
    TYPE(STACK), POINTER :: STACKP

    IF (STACKP%SET) THEN
      DEALLOCATE(STACKP%ITEMS)
      STACKP%SET = .FALSE.
    ENDIF


  END SUBROUTINE DEALLOCATESTACK



END MODULE STACKUTIL
