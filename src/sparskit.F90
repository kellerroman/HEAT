module solver
   use const
   implicit none
!   integer           :: ipar(16)    !< Parameter for solver
!   real(kind=dp)     :: fpar(16)    !< Real-Parameter for solver
!!  real(kind=dp),allocatable :: mat(:)      !< Matrix
!
!   contains
!   subroutine solve_system(n        &          
!                       ,rhs      &
!                       ,sol      &
!                       ,mat      &
!                       ,jmat     &
!                       ,imat     &
!                       ,solver)
!
!   use const
!
!   implicit none
!!!!!!!! PARAMETER   
!   integer           :: n           !< dimension of matrix
!   real(kind=dp)     :: rhs(n)      !< RHS
!   real(kind=dp)     :: sol(n)      !< Solution
!   real(kind=dp)     :: mat(:)      !< Matrix
!   integer           :: jmat(:)     !< Splate der jeweiligen Einträge
!   integer           :: imat(n+1)   !< Erstes element der N-ten Spalte
!   integer           :: solver
!   
!!!!!! LOCAL VATIABLES
!   logical           :: first_run = .true.
!
!   integer,allocatable  :: ju(:),jau(:)
!   real(kind=dp) ,allocatable     :: wk(:)
!   real(kind=dp) ,allocatable     :: au(:)
!   integer                        :: dim_mat 
!   integer                        :: dim_wk 
!   integer ,allocatable           :: iw(:)
!   
!   integer :: its,iou,i,ierr
!   real(kind=dp) :: time,res
!   real(kind=dp) :: tol = 1.0D-7
!   integer :: lfil = 6
!   integer :: nwk
!   dim_mat = ubound(mat,1)
!   dim_wk = dim_mat * 40
!   allocate (wk  (dim_wk))
!   allocate (au  (dim_mat*3))
!   allocate (ju  (dim_mat*3))
!   allocate (jau (dim_mat*3))
!   allocate (iw  (n*3))
!   ipar(2) = 2
!   ipar(3) = 1
!   ipar(4) = dim_wk
!   ipar(5) = 16
!   ipar(6) = 60
!   fpar(1) = 1.0D-8
!   fpar(2) = 1.0D-10
!   nwk = dim_mat
!
!!  if (first_run) then
!      first_run = .false. 
!      call ilut(n,mat,jmat,imat,lfil,tol    &
!               ,au,jau,ju,nwk,wk,iw,ierr)
!!      write(*,*) rhs
!!      call amux(n, sol, rhs, mat, jmat, imat)
!!      write(*,*) rhs
!!  end if
!   res = 0.0D0
!   iou = 6
!   ipar(2) = 2
! 
!   its = 0
!   ipar(1) = 0
!  ! time = dtime(dt)
!!  write(*,*) "===== SPARSKIT SOLVER ====="
!10   call bcg(n,rhs,sol,ipar,fpar,wk)
!!
!!    output the residuals
!!
!!     write(*,*) "return value", ipar(1)
!      if (ipar(7).ne.its) then
!!         write (iou, *) its, real(res)
!         its = ipar(7)
!      endif
!      res = fpar(5)
!!
!      if (ipar(1).eq.1) then
!      !   write(*,*) mat,jmat,imat
!         !write(*,*)  ipar(8),ipar(9),wk
!         call amux(n, wk(ipar(8)), wk(ipar(9)), mat, jmat, imat)
!         goto 10
!      else if (ipar(1).eq.2) then
!         call atmux(n, wk(ipar(8)), wk(ipar(9)), mat, jmat, imat)
!         goto 10
!      else if (ipar(1).eq.3 .or. ipar(1).eq.5) then
!         call lusol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
!         goto 10
!      else if (ipar(1).eq.4 .or. ipar(1).eq.6) then
!         call lutsol(n,wk(ipar(8)),wk(ipar(9)),au,jau,ju)
!         goto 10
!      else if (ipar(1).le.0) then
!         if (ipar(1).eq.0) then
!!           print *, 'Iterative sovler has satisfied convergence test.'
!         else if (ipar(1).eq.-1) then
!            print *, 'Iterative solver has iterated too many times.'
!         else if (ipar(1).eq.-2) then
!            print *, 'Iterative solver was not given enough work space.'
!            print *, 'The work space should at least have ', ipar(4),&
!                     ' elements.'
!                  stop 1
!         else if (ipar(1).eq.-3) then
!            print *, 'Iterative sovler is facing a break-down.'
!            stop 1
!         else
!            print *, 'Iterative solver terminated. code =', ipar(1)
!            stop 1
!         endif
!      endif
!     ! time2 = dtime(dt)
!!      write (iou, *) ipar(7), real(fpar(6))
!!      write (iou, *) '# retrun code =', ipar(1),&
!!          '	convergence rate =', fpar(7)
!!     write (iou, *) '# total execution time (sec)', time2-time
!!
!!     check the error
!!
!      call amux(n,sol,wk,mat,jmat,imat)
!      do i = 1, n
!         wk(n+i) = sol(i) -1.0D0
!         wk(i) = wk(i) - rhs(i)
!      enddo
!!      write (iou, *) '# the actual residual norm is', dnrm2(n,wk,1)
!!      write (iou, *) '# the error norm is', dnrm2(n,wk(1+n),1)
!!
!      if (iou.ne.6) close(iou)
!      deallocate (wk)  
!      deallocate (au) 
!      deallocate (ju)
!      deallocate (jau)
!      deallocate (iw)
!      return
!!-----end-of-runrc
!!
!end subroutine solve_system
!-----------------------------------------------------------------------
end module solver
!      function distdot(n,x,ix,y,iy)
!         use const
!      integer n, ix, iy
!      real(kind=dp) ::  distdot, x(*), y(*), ddot
!      external ddot
!      distdot = ddot(n,x,ix,y,iy)
!      return
!      end function distdot
!-----end-of-distdot
!
