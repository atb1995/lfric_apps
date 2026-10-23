!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Converts restricted (masked_value, mask) into a level index.
!>
!> @details target = nint( masked_value / mask ), clamped to
!>          [level_min, level_max], where mask > 0. Where no fine cell was
!>          valid (mask = 0) target = none_value.
module level_masked_mean_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, GH_REAL, &
                               GH_READ, GH_WRITE, CELL_COLUMN,         &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  type, public, extends(kernel_type) :: level_masked_mean_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                     &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: level_masked_mean_code
  end type

  public :: level_masked_mean_code

contains

!> @param[in]     nlayers       Number of layers (1 for a 2D field)
!> @param[in,out] level         Output level index (stored as a real)
!> @param[in]     masked_value  Restricted mask * field
!> @param[in]     mask          Restricted validity mask (fraction valid)
!> @param[in]     level_min     Lower clamp for valid cells
!> @param[in]     level_max     Upper clamp for valid cells
!> @param[in]     none_value    Value where no fine cell was valid
!> @param[in]     ndf           Number of dofs per cell
!> @param[in]     undf          Number of unique dofs
!> @param[in]     map           Dofmap
subroutine level_masked_mean_code( nlayers,                    &
                                   level, masked_value, mask,  &
                                   level_min, level_max,       &
                                   none_value,                 &
                                   ndf, undf, map )

  implicit none

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf, undf
  integer(kind=i_def), dimension(ndf),  intent(in)    :: map
  real(kind=r_def),    dimension(undf), intent(inout) :: level
  real(kind=r_def),    dimension(undf), intent(in)    :: masked_value, mask
  real(kind=r_def),    intent(in) :: level_min, level_max, none_value

  ! Valid fractions are multiples of 1/(number of fine cells), so any
  ! genuinely valid cell is far above this tolerance
  real(kind=r_def), parameter :: tol = 1.0e-6_r_def
  real(kind=r_def) :: mean
  integer(kind=i_def) :: df

  do df = 1, ndf
    if ( mask( map(df) ) > tol ) then
      mean = masked_value( map(df) ) / mask( map(df) )
      level( map(df) ) = real( nint( mean ), r_def )
      level( map(df) ) = max( level_min, min( level_max, level( map(df) ) ) )
    else
      level( map(df) ) = none_value
    end if
  end do

end subroutine level_masked_mean_code

end module level_masked_mean_kernel_mod
