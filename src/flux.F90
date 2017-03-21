module flux
   use const, only : dp

   implicit none

contains
   subroutine calc_flux(t,a,flux,dn)
      use control, only : n_BC_Cell
      implicit none
      real(kind=dp),intent(in)  :: T(1-n_BC_Cell:,1-n_BC_Cell:,1-n_BC_Cell:)
      real(kind=dp),intent(in)  :: a(1-n_BC_Cell:,1-n_BC_Cell:,1-n_BC_Cell:)
      real(kind=dp),intent(out) :: flux(:,:,:,:)
      real(kind=dp),intent(in)  :: dn(:,:,:,:)

      integer :: i,j,k
      integer :: ni,nj,nk
!      write(*,*) "CALC_RES"

      ni = ubound(T,1)-n_BC_Cell
      nj = ubound(T,2)-n_BC_Cell
      nk = ubound(T,3)-n_BC_Cell


      do k = 1,nk
         do j = 1,nj
            do i = 1,ni+1
               flux(i,j,k,1) = (a(i,j,k)+a(i-1,j,k)) &
                             * (T(i,j,k)-T(i-1,j,k)) &
                             * dn(i,j,k,1)
            end do
         end do
      end do
      do k = 1,nk
         do j = 1,nj+1
            do i = 1,ni
               flux(i,j,k,2) = (a(i,j,k)+a(i,j-1,k)) &
                             * (T(i,j,k)-T(i,j-1,k)) &
                             * dn(i,j,k,2)
            end do
         end do
      end do
      do k = 1,nk+1
         do j = 1,nj
            do i = 1,ni
               flux(i,j,k,3) = (a(i,j,k)+a(i,j,k-1)) &
                             * (T(i,j,k)-T(i,j,k-1)) &
                             * dn(i,j,k,3)
            end do
         end do
      end do
   end subroutine calc_flux

   subroutine calc_res(res,flux,edge_area)
      implicit none
      real(kind=dp),intent(out) :: res(:,:,:)
      real(kind=dp),intent(out)  :: flux(:,:,:,:)
      real(kind=dp),intent(in)  :: edge_area(:,:,:,:)

      integer :: i,j,k
      integer :: ni,nj,nk
!      write(*,*) "CALC_RES"

      ni = ubound(res,1)
      nj = ubound(res,2)
      nk = ubound(res,3)
      do k = 1,nk
         do j = 1,nj
            do i = 1,ni
               res(i,j,k) = 0.5E0_dp * &
                          ( flux(i+1,j  ,k  ,1) * edge_area(i+1,j  ,k  ,1) &
                          - flux(i  ,j  ,k  ,1) * edge_area(i  ,j  ,k  ,1) &
                          + flux(i  ,j+1,k  ,2) * edge_area(i  ,j+1,k  ,2) &
                          - flux(i  ,j  ,k  ,2) * edge_area(i  ,j  ,k  ,2) &
                          + flux(i  ,j  ,k+1,3) * edge_area(i  ,j  ,k+1,3) &
                          - flux(i  ,j  ,k  ,3) * edge_area(i  ,j  ,k  ,3) )
!               if (j == 2) then
!               write(*,'(3(I3,1X),4(ES15.8,1X) i,j,k,res(i,j,k),t(i,j,k)
!               write(*,'(6(ES15.8,1X))') t(i-1))',j,k),t(i+1,j,k) &
!                                        ,t(i,j-1,k),t(i,j+1,k) &
!                                        ,t(i,j,k-1),t(i,j,k+1)
!               write(*,'(6(ES15.8,1X))') dn(i,j,k,1),dn(i+1,j,k,1) &
!                                        ,dn(i,j,k,2),dn(i,j+1,k,2) &
!                                        ,dn(i,j,k,3),dn(i,j,k+1,3)
!
!               write(*,'(6(ES15.8,1X))')f(i+1,j  ,k  ,1) ,- f(i  ,j  ,k  ,1) &
!                          ,f(i  ,j+1,k  ,2) ,- f(i  ,j  ,k  ,2) &
!                          ,f(i  ,j  ,k+1,3) ,- f(i  ,j  ,k  ,3)
!               end if
            end do
         end do
      end do
!      i = 1
!      j = 3
!      k = 1
!      write(*,*) iter, res(1,3,1)
!      write(*,*) f(i+1,j  ,k  ,1) , edge_area(i+1,j  ,k  ,1)
!      write(*,*) -f(i  ,j  ,k  ,1) , edge_area(i  ,j  ,k  ,1)
!      write(*,*) f(i  ,j+1,k  ,2) , edge_area(i  ,j+1,k  ,2)
!      write(*,*) -f(i  ,j  ,k  ,2) , edge_area(i  ,j  ,k  ,2)
!      write(*,*) f(i  ,j  ,k+1,3) , edge_area(i  ,j  ,k+1,3)
!      write(*,*) -f(i  ,j  ,k  ,3) , edge_area(i  ,j  ,k  ,3)
   end subroutine calc_res

   subroutine res_control()
      use control, only : iter_res_out, iter,const_time_step,given_time_step,time_Step,sol_time
      use data, only : Block,nBlock,nCell
      implicit none
      integer :: b
      real(kind=dp) :: res_max,res_avg
      if (const_time_step .or. given_time_step) then
         sol_time = sol_time + time_Step
      end if
      if (mod(iter,iter_res_out) == 0 &
     .or. iter <= 50) then
         res_avg = 0.0E0_dp
         res_max = -1.0E0_dp
         do b = 1,nBlock
            res_avg = res_avg + sum(abs(block(b) % res))
            res_max = max(res_max,maxval(abs(block(b) % res)))
         end do
         res_avg = res_avg / dble(ncell)
         write(*,'(I8,1X,ES10.3,2(1X,ES12.5))') iter, sol_time,res_avg,res_max
      end if
   end subroutine res_control
end module flux
