module boundary
   use const, only : dp
   implicit none
!!!!!!!! BOUNDARY CONDITION FILE
!bc.bin
! VERSION - INTEGER
! NUMBER OF BLOCKS - INTEGER
! <BLOCK>
!   <FACE>
!     NUMBER OF BC PER FACE - INTEGER
!     <nBCpF>
!        BC_Type - INTEGER
!        dim1_start,dim1_end,dim2_start,dim2_end - INTEGER
!        {BC_DATA - REAL(:,:)}

!
contains
   subroutine update_boundary(block,nblock,inner_iter)
      use control, only : n_BC_Cell,iter
      use types
      use const
      implicit none
      integer :: nBlock
      type(tBlock) :: block(nBlock)
      integer :: inner_iter
      !< Speicherung in n-tem Variablenvektor (Für Verfahren mit mehreren Auswerteschritten (Runge-Kutta


      integer :: i,j,k,f,bc,b,nb,b1,pos
      integer :: i1,j1,k1
      integer :: ni,nj,nk
      integer :: a(4,3)


      do b = 1 , nBlock
         ni = block(b) % nCell(1)
         nj = block(b) % nCell(2)
         nk = block(b) % nCell(3)

         f = WEST_FACE
         pos = 1
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do j = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                  do nb = 1, n_BC_Cell
                     i = pos-nb
                     i1 = pos-1+nb
                     if (block(b) % face(f) % BC(bc) % BC_Type == BC_SYMMETRY) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i1,j,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_ISOTHERMAL) then
                        block(b) % T(i,j,k,inner_iter) = 2.0D0 * &
                                                   block(b) % face(f) % BC(bc) % data( &
                                                            j-block(b) % face(f) % BC(bc) % dist(1,2)+1   &
                                                           ,k-block(b) % face(f) % BC(bc) % dist(1,1)+1 ) &
                                                 - block(b) % T(i1,j,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i1,j,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type > 0 ) then
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % T(i,j,k,inner_iter) = block(b1) % T(i1,j1,k1,inner_iter)
                     else
                        !stop "BC_WEST"
                     end if
                  end do
               end do
            end do
         end do
         f = EAST_FACE
         pos = ni
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do j = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                  do nb = 1, n_BC_Cell
                     i = pos + nb
                     i1 = pos+1-nb
                     if (block(b) % face(f) % BC(bc) % BC_Type == BC_SYMMETRY) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i1,j,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_ISOTHERMAL) then
                        block(b) % T(i,j,k,inner_iter) = 2.0D0 * &
                                                   block(b) % face(f) % BC(bc) % data( &
                                                            j-block(b) % face(f) % BC(bc) % dist(1,2)+1   &
                                                           ,k-block(b) % face(f) % BC(bc) % dist(1,1)+1 ) &
                                                 - block(b) % T(i1,j,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i1,j,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type > 0 ) then
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % T(i,j,k,inner_iter) = block(b1) % T(i1,j1,k1,inner_iter)
                     else
                        stop "BC_EAST"
                     end if
                  end do
               end do
            end do
         end do
         !!!!!! NORTH & SOUTH !!!!!
         f = SOUTH_FACE
         pos = 1
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do nb = 1, n_BC_Cell
                  j = pos - nb
                  j1 = pos - 1 +nb
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                     if (block(b) % face(f) % BC(bc) % BC_Type == BC_SYMMETRY) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j1,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_ISOTHERMAL) then
                        block(b) % T(i,j,k,inner_iter) = 2.0D0 * &
                                                   block(b) % face(f) % BC(bc) % data( &
                                                            i-block(b) % face(f) % BC(bc) % dist(1,2)+1   &
                                                           ,k-block(b) % face(f) % BC(bc) % dist(1,1)+1 ) &
                                                 - block(b) % T(i,j1,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j1,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type > 0 ) then
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % T(i,j,k,inner_iter) = block(b1) % T(i1,j1,k1,inner_iter)
                     else
                        stop "BC_SOUTH"
                     end if
                  end do
               end do
            end do
         end do
         f = NORTH_FACE
         pos = nj
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do nb = 1, n_BC_Cell
                  j = pos + nb
                  j1 = pos + 1 -nb
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                     if (block(b) % face(f) % BC(bc) % BC_Type == BC_SYMMETRY) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j1,k,inner_iter)
!                        WRITE(*,*) I,J,K,block(b) % T(i,j,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_ISOTHERMAL) then
                        block(b) % T(i,j,k,inner_iter) = 2.0D0 * &
                                                   block(b) % face(f) % BC(bc) % data( &
                                                   i-block(b) % face(f) % BC(bc) % dist(1,2)+1 &
                                                  ,k-block(b) % face(f) % BC(bc) % dist(1,1)+1 ) &
                                                 - block(b) % T(i,j1,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j1,k,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type > 0 ) then
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % T(i,j,k,inner_iter) = block(b1) % T(i1,j1,k1,inner_iter)
                     else
                        stop "NORTH FACE BC UNKNOWN"
                     end if
                  end do
               end do
            end do
         end do

         f = BACK_FACE
         pos = 1
         do bc = 1,block(b) % face(f) % nBC
            do nb = 1, n_BC_Cell
               k = pos - nb
               k1 = pos - 1 + nb
               do j = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)

                     if (block(b) % face(f) % BC(bc) % BC_Type == BC_SYMMETRY) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j,k1,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_ISOTHERMAL) then
                        block(b) % T(i,j,k,inner_iter) = 2.0D0 * &
                                                   block(b) % face(f) % BC(bc) % data( &
                                                            i-block(b) % face(f) % BC(bc) % dist(1,2)+1   &
                                                           ,j-block(b) % face(f) % BC(bc) % dist(1,1)+1 ) &
                                                 - block(b) % T(i,j,k1,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j,k1,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type > 0 ) then
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % T(i,j,k,inner_iter) = block(b1) % T(i1,j1,k1,inner_iter)
                     else
                        stop "BC BACK"
                     end if
                  end do
               end do
            end do
         end do
         f = FRONT_FACE
         pos = nk
         do bc = 1,block(b) % face(f) % nBC
            do nb = 1, n_BC_Cell
               k = pos + nb
               k1 = pos + 1 - nb
               do j = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                     if (block(b) % face(f) % BC(bc) % BC_Type == BC_SYMMETRY) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j,k1,inner_iter)
                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_ISOTHERMAL) then
                        block(b) % T(i,j,k,inner_iter) = 2.0D0 * &
                                                   block(b) % face(f) % BC(bc) % data( &
                                                            i-block(b) % face(f) % BC(bc) % dist(1,2)+1   &
                                                           ,j-block(b) % face(f) % BC(bc) % dist(1,1)+1 ) &
                                                 - block(b) % T(i,j,k1,inner_iter)

                     else if (block(b) % face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
                        block(b) % T(i,j,k,inner_iter) = block(b) % T(i,j,k1,inner_iter)

                     else if (block(b) % face(f) % BC(bc) % BC_Type > 0 ) then
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % T(i,j,k,inner_iter) = block(b1) % T(i1,j1,k1,inner_iter)
                     else
                        stop "BC FRONT"
                     end if
                  end do
               end do
            end do
         end do
      end do
!      b = 1
!      DO I = 1,BLOCK(1) % nCell(1)
!!         if ( block(b) % T(i,1,1,1) /= block(b) % T(i,1,2,1) &
!!         .or. block(b) % T(i,1,1,1) /= block(b) % T(i,2,1,1) &
!!         .or. block(b) % T(i,1,1,1) /= block(b) % T(i,2,2,1) ) then
!            write(*,'(2(I0,1X),8ES10.3)') iter,i,block(b) % T(i,1,1,1) &
!                             ,block(b) % T(i,1,2,1) &
!                             ,block(b) % T(i,2,1,1) &
!                             ,block(b) % T(i,2,2,1)
!
!!         end if
!         if ( block(b) % T(i,0,1,1) /= block(b) % T(i,0,2,1) &
!         .or. block(b) % T(i,0,1,1) /= block(b) % T(i,1,0,1) &
!         .or. block(b) % T(i,0,1,1) /= block(b) % T(i,2,0,1) &
!         .or. block(b) % T(i,0,1,1) /= block(b) % T(i,3,1,1) &
!         .or. block(b) % T(i,0,1,1) /= block(b) % T(i,3,2,1) &
!         .or. block(b) % T(i,0,1,1) /= block(b) % T(i,1,3,1) &
!         .or. block(b) % T(i,0,1,1) /= block(b) % T(i,2,3,1) ) then
!            write(*,'(2(I0,1X),8ES10.3)') iter,i,block(b) % T(i,0,1,1) &
!                             ,block(b) % T(i,0,2,1) &
!                             ,block(b) % T(i,1,0,1) &
!                             ,block(b) % T(i,2,0,1) &
!                             ,block(b) % T(i,3,1,1) &
!                             ,block(b) % T(i,3,2,1) &
!                             ,block(b) % T(i,1,3,1) &
!                             ,block(b) % T(i,2,3,1)
!!            stop 1
!         end if
!      end do

   end subroutine  update_boundary

   subroutine boundary_fluxes(flux,face)
      use types
      use const
      implicit none
      real(kind=dp),intent(inout)  :: flux(:,:,:,:)
      type(tFace),intent(in) :: face(6)
      integer :: i,j,k
      integer :: ni,nj,nk
      integer :: f,bc

      ni = ubound(flux,1)
      nj = ubound(flux,2)
      nk = ubound(flux,3)

      f = WEST_FACE
      do bc = 1,face(f) % nBC
         if (face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
            i = 1
            do k = face(f) % BC(bc) % dist(1,1), face(f) % BC(bc) % dist(2,1)
               do j = face(f) % BC(bc) % dist(1,2) ,face(f) % BC(bc) % dist(2,2)
                     flux(i,j,k,1) = face(f) % BC(bc) % data(j-face(f) % BC(bc) % dist(1,2)+1 &
                                                            ,k-face(f) % BC(bc) % dist(1,1)+1 )
               end do
            end do
         end if
      end do
      f = EAST_FACE
      do bc = 1,face(f) % nBC
         if (face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
            i = ni
            do k = face(f) % BC(bc) % dist(1,1), face(f) % BC(bc) % dist(2,1)
               do j = face(f) % BC(bc) % dist(1,2) ,face(f) % BC(bc) % dist(2,2)
                     flux(i,j,k,1) = face(f) % BC(bc) % data(j-face(f) % BC(bc) % dist(1,2)+1 &
                                                            ,k-face(f) % BC(bc) % dist(1,1)+1 )
               end do
            end do
         end if
      end do
      f = SOUTH_FACE
      do bc = 1,face(f) % nBC
         if (face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
            j = 1
            do k = face(f) % BC(bc) % dist(1,1), face(f) % BC(bc) % dist(2,1)
               do i = face(f) % BC(bc) % dist(1,2) ,face(f) % BC(bc) % dist(2,2)
                     flux(i,j,k,2) = face(f) % BC(bc) % data(i-face(f) % BC(bc) % dist(1,2)+1 &
                                                            ,k-face(f) % BC(bc) % dist(1,1)+1 )
               end do
            end do
         end if
      end do
      f = NORTH_FACE
      do bc = 1,face(f) % nBC
         if (face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
            j = nj
            do k = face(f) % BC(bc) % dist(1,1), face(f) % BC(bc) % dist(2,1)
               do i = face(f) % BC(bc) % dist(1,2) ,face(f) % BC(bc) % dist(2,2)
                     flux(i,j,k,2) = face(f) % BC(bc) % data(i-face(f) % BC(bc) % dist(1,2)+1 &
                                                            ,k-face(f) % BC(bc) % dist(1,1)+1 )
               end do
            end do
         end if
      end do
      f = BACK_FACE
      do bc = 1,face(f) % nBC
         if (face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
            k = 1
            do j = face(f) % BC(bc) % dist(1,1), face(f) % BC(bc) % dist(2,1)
               do i = face(f) % BC(bc) % dist(1,2) ,face(f) % BC(bc) % dist(2,2)
                     flux(i,j,k,3) = face(f) % BC(bc) % data(i-face(f) % BC(bc) % dist(1,2)+1 &
                                                            ,j-face(f) % BC(bc) % dist(1,1)+1 )
               end do
            end do
         end if
      end do
      f = FRONT_FACE
      do bc = 1,face(f) % nBC
         if (face(f) % BC(bc) % BC_Type == BC_HEATFLUX) then
            k = nk
            do j = face(f) % BC(bc) % dist(1,1), face(f) % BC(bc) % dist(2,1)
               do i = face(f) % BC(bc) % dist(1,2) ,face(f) % BC(bc) % dist(2,2)
                     flux(i,j,k,3) = face(f) % BC(bc) % data(i-face(f) % BC(bc) % dist(1,2)+1 &
                                                            ,j-face(f) % BC(bc) % dist(1,1)+1 )
               end do
            end do
         end if
      end do
   end subroutine  boundary_fluxes
end module boundary
