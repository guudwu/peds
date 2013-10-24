! Module Linear_Ode

! Last modified on 2013-08-12 by GuuD.
! derivative_fun:
!   Function to compute derivative of the linear ODE model.
! derivative_rkf45_sub:
!   Wrapper passed to rkf45-solver.
! Set_Rkf45_Parameter_sub:
!   Set parameters of derivative_rkf45_sub.
! Ode_Result_sub:
!   Subroutine to compute numerical result of an ODE system.

module Linear_Ode_mod

  use ENV_mod
  implicit none
  private

  public &
    :: Set_Rkf45_Parameter_sub
  public &
    :: Ode_Result_sub

! Parameters passed to rkf45 solver!{{{
  real(kind=KIND_REAL) , dimension(:) , pointer &
    :: rkf45_order1
  real(kind=KIND_REAL) , dimension(:) , pointer &
    :: rkf45_order0
!}}}

  contains

  pure &
  function derivative_fun &!{{{
  ( N , X , Order1 , Order0 )

    implicit none

    real(kind=KIND_REAL) , dimension(N) &
      :: derivative_fun

    integer &
      , intent(in) &
      :: N
    real(kind=KIND_REAL) , dimension(N) &
      , intent(in) &
      :: X
    real(kind=KIND_REAL) , dimension(N*N) &
      , intent(in) &
      :: Order1
    real(kind=KIND_REAL) , dimension(N) &
      , intent(in) &
      :: Order0

    derivative_fun = &
      matmul &
      ( &
        reshape ( Order1 , (/N,N/) ) &
        , X &
      ) &
      + Order0

    return

  end function derivative_fun!}}}

  pure &
  subroutine derivative_rkf45_sub &!{{{
  ( T , X , Dx )

    implicit none

    real(kind=KIND_REAL) &
      , intent(in) &
      :: T
    real(kind=KIND_REAL) , dimension(:) &
      , intent(in) &
      :: X
    real(kind=KIND_REAL) , dimension(:) &
      , intent(out) &
      :: Dx

    Dx = &
      derivative_fun ( &
        size(X) , X , rkf45_order1 , rkf45_order0 &
      )

    return

  end subroutine derivative_rkf45_sub!}}}

  subroutine Set_Rkf45_Parameter_sub &!{{{
  ( &
    Order1 , Order0 &
  )

    implicit none

    real(kind=KIND_REAL) , dimension(:) &
      , target &
      , intent(in) &
      :: Order1
    real(kind=KIND_REAL) , dimension(:) &
      , target &
      , intent(in) &
      :: Order0

    rkf45_order1 => Order1
    rkf45_order0 => Order0

    return

  end subroutine Set_Rkf45_Parameter_sub!}}}

  subroutine Ode_Result_sub &!{{{
  ( &
    Ode_Result , &
    Dim_Ode , Dim_Time , &
    Time_Point , &
    Initial , &
    Rel_Tol , Abs_Tol , &
    Info , &
    Num_Fun &
  )

    implicit none

    real(kind=KIND_REAL) , dimension(Dim_Ode*Dim_Time) &
      , intent(out) &
      :: Ode_Result

    integer &
      , intent(in) &
      :: Dim_Ode
    integer &
      , intent(in) &
      :: Dim_Time
    real(kind=KIND_REAL) , dimension(Dim_Time) &
      , intent(in) &
      :: Time_Point
    real(kind=KIND_REAL) , dimension(Dim_Ode) &
      , intent(in) &
      :: Initial
    real(kind=KIND_REAL) &
      , intent(in) &
      :: Rel_Tol
    real(kind=KIND_REAL) &
      , intent(in) &
      :: Abs_Tol
    integer , dimension(Dim_Time-1) &
      , intent(out) &
      :: Info
    integer , dimension(Dim_Time-1) &
      , intent(out) &
      :: Num_Fun

! Local Variables
    real(kind=KIND_REAL) &
      :: time
    integer &
      :: ind
    integer &
      :: iflag
    real(kind=KIND_REAL) , dimension(Dim_Ode) &
      :: state
    real(kind=KIND_REAL) , dimension(3+6*Dim_Ode) &
      :: work_r
    integer , dimension(5) &
      :: work_i

    iflag = 1
    time = Time_Point ( 1 )
    state = Initial

    Ode_Result ( 1 : Dim_Ode ) = Initial

    do ind = 1 , (Dim_Time-1)
      call rkf45 ( &
        derivative_rkf45_sub , &
        Dim_Ode , state , &
        time , Time_Point ( ind + 1 ) , &
        Rel_Tol , Abs_Tol , &
        iflag , &
        work_r , &
        work_i &
      )
      Ode_Result ( (ind*Dim_Ode+1) : ((ind+1)*Dim_Ode) ) &
        = state
      Info ( ind ) = iflag
      Num_Fun ( ind ) = work_i ( 1 )
      work_i ( 1 ) = 0
    end do

    return

  end subroutine Ode_Result_sub!}}}

end module Linear_Ode_mod
