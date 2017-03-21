module control
   use const, only : dp
   implicit none
   character(len=*), parameter :: file_git_in   = "data_in.cgns"
   character(len=*), parameter :: file_sol_out  = "data_out.cgns"
   character(len=*), parameter :: file_bc       = "bc.bin"
   character(len=*), parameter :: file_sol_in   = "data_in.cgns"
   integer :: n_BC_Cell = 1
   integer :: Dimen
   integer :: nFace
   integer :: nCorner
   integer :: iter
   integer :: inner_iter
   integer :: max_iter = 50
   integer :: n_inner_iter = 1
   integer :: iter_sol_out = 10
   integer :: iter_res_out = 1000
   logical :: sol_out = .false.
!   logical :: given_time_step = .false.
   logical :: given_time_step = .true.
   logical :: const_time_step = .true.
   real(kind=dp) :: time_step = 1.0E-5_dp
   real(kind=dp) :: sol_time = 0.0E0_dp
   real(kind=dp) :: CFL = 1.0E0_dp

   logical :: implicit = .true.
!   logical :: implicit = .false.
end module control
