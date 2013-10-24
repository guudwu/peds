! Test rkf45-solver

! Last modified on 2013-10-23 by GuuD WU.

! Read ODE parameters from data-files,
!   compute the numerical result,
!   then compare with real observations.

program main

  use ENV_mod
  use Linear_Ode_mod
  implicit none

  integer(kind=KIND_INT) &
    :: NUM_UNIT

  integer(kind=KIND_INT) &
    :: Dim_Ode
  integer(kind=KIND_INT) &
    :: Dim_Time
  real(kind=KIND_REAL) , dimension(:) , allocatable &
    :: Time_Point
  real(kind=KIND_REAL) , dimension(:) , allocatable &
    :: Order1
  real(kind=KIND_REAL) , dimension(:) , allocatable &
    :: Order0
  real(kind=KIND_REAL) , dimension(:) , allocatable &
    :: Observation

  integer(kind=KIND_INT) , dimension(:) , allocatable &
    :: Info
  integer(kind=KIND_INT) , dimension(:) , allocatable &
    :: Max_Fun


  real(kind=KIND_REAL) &
    :: Tol = 1e-5

  real(kind=KIND_REAL) , dimension(:) , allocatable &
    :: Ode_Result

  open ( newunit = NUM_UNIT , file = 'config.data' ,&
    action = 'read' , status = 'old' )
  read ( NUM_UNIT , * ) Dim_Ode , Dim_Time
  allocate ( Time_Point(Dim_Time) )
  read ( NUM_UNIT , * ) Time_Point
  close ( unit = NUM_UNIT )

  allocate ( Order1(Dim_Ode**2) )
  open ( newunit = NUM_UNIT , file = 'coefficient.data' ,&
    action = 'read' , status = 'old' )
  read ( NUM_UNIT , * ) Order1
  close ( unit = NUM_UNIT )

  allocate ( Order0(Dim_Ode) )
  open ( newunit = NUM_UNIT , file = 'intercept.data' ,&
    action = 'read' , status = 'old' )
  read ( NUM_UNIT , * ) Order0
  close ( unit = NUM_UNIT )

  allocate ( Observation(Dim_Ode*Dim_Time) )
  open ( newunit = NUM_UNIT , file = 'observation.data' ,&
    action = 'read' , status = 'old' )
  read ( NUM_UNIT , * ) Observation
  close ( unit = NUM_UNIT )

  call Set_Rkf45_Parameter_sub ( Order1 , Order0 )

  allocate ( Ode_Result(Dim_Ode*Dim_Time) )
  allocate ( Info(Dim_Time-1) )
  allocate ( Max_Fun(Dim_Time-1) )

  Info = 0
  Max_Fun = 0

  call Ode_Result_sub &
    ( &
      Ode_Result , &
      Dim_Ode , Dim_Time , &
      Time_Point , &
      Observation ( 1 : Dim_Ode ) , &
      Tol , Tol , &
      Info , &
      Max_Fun &
    )
  write(*,*) maxval ( abs ( Ode_Result - Observation ) )

  write(*,*) Info
  write(*,*) Max_Fun

  deallocate(Time_Point)
  deallocate(Order1)
  deallocate(Order0)
  deallocate(Observation)
  deallocate(Ode_result)
  deallocate(Info)
  deallocate(Max_Fun)

  stop

end program main
