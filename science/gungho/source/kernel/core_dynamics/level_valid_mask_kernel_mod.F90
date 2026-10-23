!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Builds a validity mask and masked values for a 2D level-index field
!>        so that it can be restricted as a mean over valid cells only.
!>
!> @details A cell is valid if valid_field > threshold. The outputs are
!>          mask = 1 (valid) or 0 (not valid), and masked_value = mask * field.
!>          Restricting both and dividing gives the mean over valid cells.
module level_valid_mask_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, GH_REAL, &
                               GH_READ, GH_WRITE, CELL_COLUMN,         &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  type, public, extends(kernel_type) :: level_valid_mask_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                     &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: level_valid_mask_code
  end type

  public :: level_valid_mask_code

contains

!> @param[in]     nlayers       Number of layers (1 for a 2D field)
!> @param[in,out] mask          1.0 where valid, 0.0 elsewhere
!> @param[in,out] masked_value  mask * field
!> @param[in]     field         Field whose values are to be averaged
!> @param[in]     valid_field   Field that decides validity
!> @param[in]     threshold     A cell is valid if valid_field > threshold
!> @param[in]     ndf           Number of dofs per cell
!> @param[in]     undf          Number of unique dofs
!> @param[in]     map           Dofmap
subroutine level_valid_mask_code( nlayers,                    &
                                  mask, masked_value,         &
                                  field, valid_field,         &
                                  threshold,                  &
                                  ndf, undf, map )

  implicit none

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf, undf
  integer(kind=i_def), dimension(ndf),  intent(in)    :: map
  real(kind=r_def),    dimension(undf), intent(inout) :: mask, masked_value
  real(kind=r_def),    dimension(undf), intent(in)    :: field, valid_field
  real(kind=r_def),    intent(in) :: threshold

  integer(kind=i_def) :: df

  do df = 1, ndf
    if ( valid_field( map(df) ) > threshold ) then
      mask( map(df) )         = 1.0_r_def
      masked_value( map(df) ) = field( map(df) )
    else
      mask( map(df) )         = 0.0_r_def
      masked_value( map(df) ) = 0.0_r_def
    end if
  end do

end subroutine level_valid_mask_code

end module level_valid_mask_kernel_mod
