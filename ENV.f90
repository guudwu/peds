! Module ENV

! Last modified on 2013-10-22 by GuuD WU.

! The module contains two public KIND parameters.

module ENV_mod

  implicit none
  private

! Public!{{{
  public :: KIND_INT
  public :: KIND_REAL
!}}}

! Variable!{{{
  double precision :: &
    dp
  integer , parameter :: &
    KIND_INT = 4
  integer , parameter :: &
    KIND_REAL = selected_real_kind ( precision(dp) )
!}}}

end module ENV_mod
