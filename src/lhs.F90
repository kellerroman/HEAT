module lhs
   use const
   use data
   implicit none
   contains

   subroutine update_implicit ( mat, rhs, sol, block)
      implicit none
      real(kind=dp), intent(out) :: mat (:)
      real(kind=dp), intent(out) :: rhs(:)
      real(kind=dp), intent(out) :: sol(:)
      type(tblock), intent(in) :: block(:)


      integer :: i,j,k,b,n,pos,pos1
      pos = 0
      n = 0
      do b = 1,ubound(block,1)
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nCell(2)
               do i = 1, block(b) % nCell(1)
                  pos = pos + 1
                  n = n + 1
                  pos1 = pos
                  mat(pos) = block(b) %  mat_geo(0,i,j,k)
                  if (k > 1) then
                     pos = pos + 1
                     mat(pos) = block(b) %  mat_geo(-3,i,j,k) * block(b) % dt (i,j,k)
                  else
                     if (block(b) % face(5) % BC(1) % BC_Type == BC_SYMMETRY) then
                        mat(pos1) = mat(pos1) +  block(b) %  mat_geo(-3,i,j,k)
                     else
                        write(*,*)  "BC Back LHS unkown"
                        stop 1
                     end if
                  end if
                  if (j > 1) then
                     pos = pos + 1
                     mat(pos) = block(b) %  mat_geo(-2,i,j,k) * block(b) % dt (i,j,k)
                  else
                     if (block(b) % face(3) % BC(1) % BC_Type == BC_SYMMETRY) then
                        mat(pos1) = mat(pos1) +  block(b) %  mat_geo(-2,i,j,k)
                     else
                        write(*,*)  "BC south LHS unkown"
                        stop 1
                     end if
                  end if
                  if (i > 1) then
                     pos = pos + 1
                     mat(pos) = block(b) %  mat_geo(-1,i,j,k) * block(b) % dt (i,j,k)
                  else
                     if (block(b) % face(1) % BC(1) % BC_Type == BC_SYMMETRY) then
                        mat(pos1) = mat(pos1) +  block(b) %  mat_geo(-1,i,j,k)
                     else
                        write(*,*)  "BC west LHS unkown"
                        stop 1
                     end if
                  end if
                  if (i < block(b) % nCell(1)) then
                     pos = pos + 1
                     mat(pos) = block(b) %  mat_geo(1,i,j,k) * block(b) % dt (i,j,k)
                  else
                     if (block(b) % face(2) % BC(1) % BC_Type == BC_SYMMETRY) then
                        mat(pos1) = mat(pos1) +  block(b) %  mat_geo(1,i,j,k)
                     else
                        write(*,*)  "BC east LHS unkown"
                        stop 1
                     end if
                  end if
                  if (j < block(b) % nCell(2)) then
                     pos = pos + 1
                     mat(pos) = block(b) %  mat_geo(2,i,j,k) * block(b) % dt (i,j,k)
                  else
                     if (block(b) % face(4) % BC(1) % BC_Type == BC_SYMMETRY) then
                        mat(pos1) = mat(pos1) +  block(b) %  mat_geo(2,i,j,k)
                     else
                        write(*,*)  "BC north LHS unkown"
                        stop 1
                     end if
                  end if
                  if (k < block(b) % nCell(3)) then
                     pos = pos + 1
                     mat(pos) = block(b) %  mat_geo(3,i,j,k) * block(b) % dt (i,j,k)
                  else
                     if (block(b) % face(6) % BC(1) % BC_Type == BC_SYMMETRY) then
                        mat(pos1) = mat(pos1) +  block(b) %  mat_geo(3,i,j,k)
                     else
                        write(*,*)  "BC Front LHS unkown"
                        stop 1
                     end if
                  end if
                  mat(pos1) = mat(pos1) * block(b) % dt (i,j,k) + 1.0E0_dp
                  rhs(n) = block(b) % T(i,j,k,1)
                  sol(n) = block(b) % T(i,j,k,1)
               end do
            end do
         end do

      end do
!      k = 0
!      do n=1,n_nonzero_elements
!         if (n == imat(k+1) ) k = k + 1
!         write(*,*) n,mat(n),jmat(n),k
!      end do
!    call dump (1,nCell,.true.,mat,jmat,imat,6)
!      write(*,*) imat(k+1),nCell,k
!      do n = 1 , nCell
!        write(*,*) sol(n),rhs(n)
!     end do
!     stop
      
!      do n=1,n_nonzero_elements
!         if (row_pos(n) == col_pos(n).and. col_pos(n) /= 0) then
!            write(*,*) row_pos(n), col_pos(n),mat(n), rhs(col_pos(n))
!         else
!            write(*,*) row_pos(n), col_pos(n),mat(n)
!         end if
!     Cell1

   end subroutine update_implicit

   subroutine update_sol_imp(rhs,block)
      implicit none
      real(kind=dp), intent(in) :: rhs(:)
      type(tblock), intent(inout) :: block(:)
      integer :: i,j,k,b,n,pos,pos1
      pos = 0
      n = 0
      do b = 1,ubound(block,1)
         do k = 1, block(b) % nCell(3)
            do j = 1, block(b) % nCell(2)
               do i = 1, block(b) % nCell(1)
                  n = n + 1
                  block(b) % T(i,j,k,1) = rhs(n)
               end do
            end do
         end do
      end do

   end subroutine update_sol_imp


end module lhs
