module time_int
   use const, only : dp
   implicit none

contains
   subroutine update_sol(T,res,dt,cell_vol)
      use control, only : n_BC_Cell
      implicit none
      real(kind=dp),intent(inout)  :: T(1-n_BC_Cell:,1-n_BC_Cell:,1-n_BC_Cell:)
      real(kind=dp),intent(in)  :: res(:,:,:)
      real(kind=dp),intent(in)  :: dt(:,:,:)
      real(kind=dp),intent(in)  :: cell_vol(:,:,:)

      integer :: i,j,k
      integer :: ni,nj,nk

      ni = ubound(T,1)-n_BC_Cell
      nj = ubound(T,2)-n_BC_Cell
      nk = ubound(T,3)-n_BC_Cell

      do k = 1,nk
         do j = 1,nj
            do i = 1,ni
               T(i,j,k) = T(i,j,k) + dt(i,j,k)  * res(i,j,k) * cell_vol(i,j,k)
            end do
         end do
      end do
!      write(*,*) T(0,3,1),T(1,3,1),dt(1,3,1),res(1,3,1),T(1,3,0),T(1,3,2),cell_vol(1,3,1)
   end subroutine update_sol

   subroutine calc_timestep(dt,a,dn2_dt)
      use const, only : dp
      use control, only : time_step,CFL,const_time_step,given_time_step
      implicit none
      real(kind=dp),intent(out)  :: dt(:,:,:)
      real(kind=dp),intent(in)   :: a(:,:,:)
      real(kind=dp),intent(in)   :: dn2_dt(:,:,:)

      integer :: i,j,k
      integer :: ni,nj,nk

      ni = ubound(dt,1)
      nj = ubound(dt,2)
      nk = ubound(dt,3)
      if (.not.given_time_step) then
         do k = 1, nk
            do j = 1,nj
               do i = 1,ni
                  dt(i,j,k) = CFL * dn2_dt(i,j,k) / (a(i,j,k))
               end do
            end do
         end do
         if (const_time_step) then
            time_step = min(minval(dt) , time_step)
         end if
      end if
   end subroutine calc_timestep

   subroutine set_timestep(dt)
      use const, only : dp
      use control, only : time_step,const_time_step,given_time_step
      implicit none
      real(kind=dp),intent(out)  :: dt(:,:,:)
      !!! time_step was calculated in calc_timestep or is set by control
      if (const_time_step .or. given_time_step) then
         dt = time_step
      end if
   end subroutine set_timestep
end module time_int
