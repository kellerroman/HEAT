module data
   use types
   use const
   use sparse_mat
   implicit none
   integer :: nBlock
   type(tblock), allocatable     :: block(:)

   real(kind = dp), allocatable  :: mat(:)
   real(kind = dp), allocatable  :: mat2(:)
   real(kind = dp), allocatable  :: rhs(:)
   real(kind = dp), allocatable  :: sol(:)
   integer        , allocatable  :: row_pos_save(:)
   integer        , allocatable  :: col_pos_save(:)
   integer        , allocatable  :: imat(:)
   !< Gibt jeweils das erste Element einer Zeile an
   integer        , allocatable  :: jmat(:)
   !< Gibt an in welcher Matrix Spalte der Eintrag steht
   integer                       :: n_nonzero_elements
   integer                       :: nCell
contains
   subroutine init()
      use control, only: n_inner_iter,Dimen,n_BC_Cell,time_step,implicit
      implicit none

      integer :: b,i,j,k,n
      integer :: b1,pos,f,bc,col
      integer :: i1,j1,k1
      integer :: i2,j2,k2
      integer :: n_mat_element,dci,dcj,dck
      real(kind=dp) :: d1(3),d2(3),v1(3)
      integer :: a(4,3)

      write(*,*) "INIT"
      n_mat_element = 0
      do b = 1,nBlock


         allocate (block(b) % a &
                  ( 1-n_BC_Cell:block(b) % nCell(1)+n_BC_Cell &
                  , 1-n_BC_Cell:block(b) % nCell(2)+n_BC_Cell &
                  , 1-n_BC_Cell:block(b) % nCell(3)+n_BC_Cell &
                  , n_inner_iter ) )

         allocate ( block(b) % res &
                  (1 : block(b)%nCell (1)            &
                  ,1 : block(b)%nCell (2)            &
                  ,1 : block(b)%nCell (3)            &
                  ,n_inner_iter ) )

         allocate (block(b) % cell_vol &
                  ( block(b) % nCell(1) &
                  , block(b) % nCell(2) &
                  , block(b) % nCell(3) ) )

         allocate (block(b) % schwerpunkt &
                  ( 1-n_BC_Cell:block(b) % nCell(1)+n_BC_Cell &
                  , 1-n_BC_Cell:block(b) % nCell(2)+n_BC_Cell &
                  , 1-n_BC_Cell:block(b) % nCell(3)+n_BC_Cell &
                  , Dimen ) )

         allocate ( block(b) % Edge_area &
                  (1 : block(b)%nPkt (1)            &
                  ,1 : block(b)%nPkt (2)            &
                  ,1 : block(b)%nPkt (3)            &
                  ,Dimen ) )

         allocate ( block(b) % flux &
                  (1 : block(b)%nPkt (1)            &
                  ,1 : block(b)%nPkt (2)            &
                  ,1 : block(b)%nPkt (3)            &
                  ,Dimen ) )
         allocate ( block(b) % dn &
                  (1            : block(b)%nPkt (1)            &
                  ,1            : block(b)%nPkt (2)            &
                  ,1            : block(b)%nPkt (3)            &
                  ,Dimen ) )

         allocate ( block(b) % dt &
                  (1 : block(b)%nCell (1)            &
                  ,1 : block(b)%nCell (2)            &
                  ,1 : block(b)%nCell (3) ))

         allocate ( block(b) % dn2_dt &
                  (1 : block(b)%nCell (1)            &
                  ,1 : block(b)%nCell (2)            &
                  ,1 : block(b)%nCell (3) ))

         if (implicit) then

            allocate (block(b) % mat_geo(-dimen:dimen   &
                     ,1 : block(b)%nCell (1)            &
                     ,1 : block(b)%nCell (2)            &
                     ,1 : block(b)%nCell (3) ))
         end if
         block(b) % a = 385.0E0_dp
         block(b) % T = 200
         block(b) % dt = time_step
      end do

      do b = 1,nBlock
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Zellseiten in I-Richtung !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do k = 1,block(b) % nCell(3)
            do j = 1,block(b) % nCell(2)
               do i = 1,block(b) % nPkt(1)
                  d1 = block(b) % xyz(i  ,j+1,k+1,:) - block(b) % xyz(i  ,j  ,k  ,:)
                  d2 = block(b) % xyz(i  ,j  ,k+1,:) - block(b) % xyz(i  ,j+1,k  ,:)
                  block(b) % Edge_area(i,j,k,1) = 0.5E0_dp*len_cross(d1,d2)
               end do
            end do
         end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Zellseiten in J-Richtung !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do k = 1,block(b) % nCell(3)
            do j = 1,block(b) % nPkt(2)
               do i = 1,block(b) % nCell(1)
                  d1 = block(b) % xyz(i+1,j  ,k+1,:) - block(b) % xyz(i  ,j  ,k  ,:)
                  d2 = block(b) % xyz(i  ,j  ,k+1,:) - block(b) % xyz(i+1,j  ,k  ,:)
                  block(b) % Edge_area(i,j,k,2) = 0.5E0_dp*len_cross(d1,d2)
               end do
            end do
         end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Zellseiten in K-Richtung !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do k = 1,block(b) % nPkt(3)
            do j = 1,block(b) % nCell(2)
               do i = 1,block(b) % nCell(1)
                  d1 = block(b) % xyz(i+1,j+1,k  ,:) - block(b) % xyz(i  ,j  ,k  ,:)
                  d2 = block(b) % xyz(i+1,j  ,k  ,:) - block(b) % xyz(i  ,j+1,k  ,:)
                  block(b) % Edge_area(i,j,k,3) = 0.5E0_dp*len_cross(d1,d2)
               end do
            end do
         end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!! Schwerpunkt und Volumen der Zellen !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do k = 1,block(b) % nCell(3)
            do j = 1,block(b) % nCell(2)
               do i = 1,block(b) % nCell(1)
                  block(b) % schwerpunkt(i,j,k,:) =  ( block(b) % xyz(i  ,j  ,k  ,:) &
                                                     + block(b) % xyz(i+1,j  ,k  ,:) &
                                                     + block(b) % xyz(i  ,j+1,k  ,:) &
                                                     + block(b) % xyz(i+1,j+1,k  ,:) &
                                                     + block(b) % xyz(i  ,j  ,k+1,:) &
                                                     + block(b) % xyz(i+1,j  ,k+1,:) &
                                                     + block(b) % xyz(i  ,j+1,k+1,:) &
                                                     + block(b) % xyz(i+1,j+1,k+1,:) &
                                                     ) * 0.125E0_dp
                  d1 = block(b) % xyz(i  ,j+1,k+1,:) - block(b) % xyz(i  ,j  ,k  ,:)
                  d2 = block(b) % xyz(i  ,j  ,k+1,:) - block(b) % xyz(i  ,j+1,k  ,:)
                  v1 = abs_cross(d1,d2)
                  d1 = block(b) % xyz(i+1,j  ,k+1,:) - block(b) % xyz(i  ,j  ,k  ,:)
                  d2 = block(b) % xyz(i  ,j  ,k+1,:) - block(b) % xyz(i+1,j  ,k  ,:)
                  v1 = v1 + abs_cross(d1,d2)
                  d1 = block(b) % xyz(i+1,j+1,k  ,:) - block(b) % xyz(i  ,j  ,k  ,:)
                  d2 = block(b) % xyz(i+1,j  ,k  ,:) - block(b) % xyz(i  ,j+1,k  ,:)
                  v1 = v1 + abs_cross(d1,d2)
                  d1 = block(b) % xyz(i+1,j+1,k+1,:) - block(b) % xyz(i  ,j  ,k  ,:)

                  block(b) % cell_vol(i,j,k) = 6.0E0_dp / abs(dot_product(v1,d1))
!                  write(*,'(4(I3,1X),4(ES12.5,1X))') b,i,j,k,block(b) % cell_vol(i,j,k),block(b) % Edge_area(i,j,k,:)
               end do
            end do
         end do
      end do
      do b = 1, nBlock
         f = WEST_FACE
         pos = 1
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do j = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                  do n = 1, n_BC_Cell
                     i = pos - n
                     if (block(b) % face(f) % BC(bc) % BC_Type < 0 ) then
                     i1 = pos + 1 - n
                     i2 = pos + 2 - n
                     block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
                                                       * block(b) % schwerpunkt(i1,j,k,:) &
                                                       - block(b) % schwerpunkt(i2,j,k,:)
                     else
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % schwerpunkt(i,j,k,:) = block(b1) % schwerpunkt(i1,j1,k1,:)
                     end if
                  end do
               end do
            end do
         end do
         f = EAST_FACE
         pos = block(b) % nCell(1)
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do j = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                  do n=1,n_BC_Cell
                     i = pos + n
                     if (block(b) % face(f) % BC(bc) % BC_Type < 0 ) then
                        i1 = pos - 1 + n
                        i2 = pos - 2 + n
                        block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
                                                          * block(b) % schwerpunkt(i1,j,k,:) &
                                                          - block(b) % schwerpunkt(i2,j,k,:)
                     else
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % schwerpunkt(i,j,k,:) = block(b1) % schwerpunkt(i1,j1,k1,:)
                     end if
                  end do
               end do
            end do
         end do
         f = SOUTH_FACE
         pos = 1
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do n=1,n_BC_Cell
                  j = pos - n
                  j1 = pos + 1 - n
                  j2 = pos + 2 - n
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                     if (block(b) % face(f) % BC(bc) % BC_Type < 0 ) then
                        block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
                                                          * block(b) % schwerpunkt(i,j1,k,:) &
                                                          - block(b) % schwerpunkt(i,j2,k,:)
                     else
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % schwerpunkt(i,j,k,:) = block(b1) % schwerpunkt(i1,j1,k1,:)
                     end if
                  end do
               end do
            end do
         end do

         f = NORTH_FACE
         pos = block(b) % nCell(2)
         do bc = 1,block(b) % face(f) % nBC
            do k = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
               do n=1,n_BC_Cell
                  j = pos + n
                  j1 = pos - 1 + n
                  j2 = pos - 2 + n
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                     if (block(b) % face(f) % BC(bc) % BC_Type < 0 ) then
                        block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
                                                          * block(b) % schwerpunkt(i,j1,k,:) &
                                                          - block(b) % schwerpunkt(i,j2,k,:)
                     else
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % schwerpunkt(i,j,k,:) = block(b1) % schwerpunkt(i1,j1,k1,:)
                     end if
                  end do
               end do
            end do
         end do

         f = BACK_FACE
         pos = 1
         do bc = 1,block(b) % face(f) % nBC
            do n = 1, n_BC_Cell
               k = pos - n
               k1 = pos + 1 - n
               k2 = pos + 2 - n
               do j = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                     if (block(b) % face(f) % BC(bc) % BC_Type < 0 ) then
                        if (block(b) % nCell(3) > 1) then
                           block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
                                                             * block(b) % schwerpunkt(i,j,k1,:) &
                                                             - block(b) % schwerpunkt(i,j,k2,:)
                        else
                           block(b) % schwerpunkt(i,j,k,:) = block(b) % schwerpunkt(i,j,1,:)
                           block(b) % schwerpunkt(i,j,k,3) = block(b) % schwerpunkt(i,j,1,3)-dble(n)
                        end if
                     else
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % schwerpunkt(i,j,k,:) = block(b1) % schwerpunkt(i1,j1,k1,:)
                     end if
                  end do
               end do
            end do
         end do

         f = FRONT_FACE
         pos = block(b) % nCell(3)
         do bc = 1,block(b) % face(f) % nBC
            do n = 1, n_BC_Cell
               k = pos + n
               k1 = pos - 1 + n
               k2 = pos - 2 + n
               do j = block(b) % face(f) % BC(bc) % dist(1,1),block(b) % face(f) % BC(bc) % dist(2,1)
                  do i = block(b) % face(f) % BC(bc) % dist(1,2),block(b) % face(f) % BC(bc) % dist(2,2)
                     if (block(b) % face(f) % BC(bc) % BC_Type < 0 ) then
                        if (block(b) % nCell(3) > 1) then
                           block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
                                                             * block(b) % schwerpunkt(i,j,k1,:) &
                                                             - block(b) % schwerpunkt(i,j,k2,:)
                        else
                           block(b) % schwerpunkt(i,j,k,:) = block(b) % schwerpunkt(i,j,1,:)
                           block(b) % schwerpunkt(i,j,k,3) = block(b) % schwerpunkt(i,j,1,3)+dble(n)
                        end if
                     else
                        b1 = block(b) % face(f) % BC(bc) % BC_Type
                        a = block(b) % face(f) % BC(bc) % a
                        i1 = a(1,1) + a(2,1) * i + a(3,1) * j+a(4,1) * k
                        j1 = a(1,2) + a(2,2) * i + a(3,2) * j+a(4,2) * k
                        k1 = a(1,3) + a(2,3) * i + a(3,3) * j+a(4,3) * k
                        block(b) % schwerpunkt(i,j,k,:) = block(b1) % schwerpunkt(i1,j1,k1,:)
                     end if
                  end do
               end do
            end do
         end do


!         if (block(b) % nCell(3) > 1) then
!            pos = 1
!            do n=1,n_BC_Cell
!               k = pos - n
!               k1 = pos + 1 - n
!               k2 = pos + 2 - n
!               do j = 1,block(b) % nCell(2)
!                  do i = 1,block(b) % nCell(1)
!                     block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
!                                                       * block(b) % schwerpunkt(i,j,k1,:) &
!                                                       - block(b) % schwerpunkt(i,j,k2,:)
!
!!                     write(*,*) j,k,block(b) % schwerpunkt(i,j,:,3)
!                  end do
!               end do
!            end do
!            pos = block(b) % nCell(3)
!            do n=1,n_BC_Cell
!               k = pos + n
!               k1 = pos - 1 + n
!               k2 = pos - 2 + n
!               do j = 1,block(b) % nCell(2)
!                  do i = 1,block(b) % nCell(1)
!                     block(b) % schwerpunkt(i,j,k,:) = 2.0E0_dp &
!                                                       * block(b) % schwerpunkt(i,j,k1,:) &
!                                                       - block(b) % schwerpunkt(i,j,k2,:)
!                  end do
!               end do
!            end do
!
!
!            do j = 1,block(b) % nCell(2)
!               do i = 1,block(b) % nCell(1)
!                  block(b) % schwerpunkt(i,j,0,:) = block(b) % schwerpunkt(i,j,1,:)
!                  block(b) % schwerpunkt(i,j,2,:) = block(b) % schwerpunkt(i,j,1,:)
!                  block(b) % schwerpunkt(i,j,0,3) = 0.0E0_dp
!                  block(b) % schwerpunkt(i,j,2,3) = 2.0E0_dp * block(b) % schwerpunkt(i,j,1,3)
!                  do n=2,n_BC_Cell
!                     k = 1
!                     block(b) % schwerpunkt(i,j,k-n,:) = 2.0E0_dp &
!                                                       * block(b) % schwerpunkt(i,j,k+1-n,:) &
!                                                       - block(b) % schwerpunkt(i,j,k+2-n,:)
!                     k = block(b) % nCell(3)
!                     block(b) % schwerpunkt(i,j,k+n,:) = 2.0E0_dp &
!                                                       * block(b) % schwerpunkt(i,j,k-1+n,:) &
!                                                       - block(b) % schwerpunkt(i,j,k-2+n,:)
!
!                  end do
!               end do
!         end do
!
!
!         end if
         do k = 1,block(b) % nCell(3)
            do j = 1,block(b) % nCell(2)
               do i = 1,block(b) % nPkt(1)
                  block(b) % dn(i,j,k,1) = 1.0E0_dp / vec_len( block(b) % schwerpunkt(i  ,j  ,k  ,:) &
                                                             - block(b) % schwerpunkt(i-1,j  ,k  ,:) )
               end do
            end do
         end do
         do k = 1,block(b) % nCell(3)
            do j = 1,block(b) % nPkt(2)
               do i = 1,block(b) % nCell(1)
                  block(b) % dn(i,j,k,2) = 1.0E0_dp / vec_len( block(b) % schwerpunkt(i  ,j  ,k  ,:) &
                                                             - block(b) % schwerpunkt(i  ,j-1,k  ,:) )
               end do
            end do
         end do
         do k = 1,block(b) % nPkt(3)
            do j = 1,block(b) % nCell(2)
               do i = 1,block(b) % nCell(1)
                  block(b) % dn(i,j,k,3) = 1.0E0_dp / vec_len( block(b) % schwerpunkt(i  ,j  ,k  ,:) &
                                                             - block(b) % schwerpunkt(i  ,j  ,k-1,:) )
               end do
            end do
         end do
         block(b) % dn2_dt = 0.0D0
         do n = 1, Dimen
            do k = 1,block(b) % nCell(3)
               do j = 1,block(b) % nCell(2)
                  do i = 1,block(b) % nCell(1)
                     d1(1) =     max( block(b) % xyz(i  ,j  ,k  ,n) &
                                    , block(b) % xyz(i+1,j  ,k  ,n) &
                                    , block(b) % xyz(i  ,j+1,k  ,n) &
                                    , block(b) % xyz(i+1,j+1,k  ,n) &
                                    , block(b) % xyz(i  ,j  ,k+1,n) &
                                    , block(b) % xyz(i+1,j  ,k+1,n) &
                                    , block(b) % xyz(i  ,j+1,k+1,n) &
                                    , block(b) % xyz(i+1,j+1,k+1,n) )

                     d1(2) =     min( block(b) % xyz(i  ,j  ,k  ,n) &
                                    , block(b) % xyz(i+1,j  ,k  ,n) &
                                    , block(b) % xyz(i  ,j+1,k  ,n) &
                                    , block(b) % xyz(i+1,j+1,k  ,n) &
                                    , block(b) % xyz(i  ,j  ,k+1,n) &
                                    , block(b) % xyz(i+1,j  ,k+1,n) &
                                    , block(b) % xyz(i  ,j+1,k+1,n) &
                                    , block(b) % xyz(i+1,j+1,k+1,n) )

                     d1(3) = (d1(1) - d1(2) )

                     block(b) % dn2_dt(i,j,k) = block(b) % dn2_dt(i,j,k) &
                                              + 1.0E0_dp / (d1(3)*d1(3))
                  end do
               end do
            end do
         end do
         block(b) % dn2_dt = 0.5E0_dp / block(b) % dn2_dt
         if (implicit) then
            do k = 1, block(b) % nCell(3)
               do j = 1,block(b) % nCell(2)
                  do i = 1, block(b) % nCell(1)
                     block(b) % mat_geo(-3,i,j,k) = - block(b) % dn       (i  ,j  ,k  ,3 ) &
                                                    * block(b) % edge_area(i  ,j  ,k  ,3 ) &
                                                    * block(b) % cell_vol (i  ,j  ,k     ) &
                                                    * block(b) % a        (i  ,j  ,k  ,1 )

                     block(b) % mat_geo(-2,i,j,k) = - block(b) % dn       (i  ,j  ,k  ,2 ) &
                                                    * block(b) % edge_area(i  ,j  ,k  ,2 ) &
                                                    * block(b) % cell_vol (i  ,j  ,k     ) &
                                                    * block(b) % a        (i  ,j  ,k  ,1 )

                     block(b) % mat_geo(-1,i,j,k) = - block(b) % dn       (i  ,j  ,k  ,1 ) &
                                                    * block(b) % edge_area(i  ,j  ,k  ,1 ) &
                                                    * block(b) % cell_vol (i  ,j  ,k     ) &
                                                    * block(b) % a        (i  ,j  ,k  ,1 )

                     !!! CENTRAL ELEMENT: PLUS 1 is done when multiplying with delta_t

                     block(b) % mat_geo( 0,i,j,k) = ( &
                                                    block(b) % dn       (i  ,j  ,k  ,1 ) &
                                                  * block(b) % edge_area(i  ,j  ,k  ,1 ) &
                                                  + block(b) % dn       (i  ,j  ,k  ,2 ) &
                                                  * block(b) % edge_area(i  ,j  ,k  ,2 ) &
                                                  + block(b) % dn       (i  ,j  ,k  ,3 ) &
                                                  * block(b) % edge_area(i  ,j  ,k  ,3 ) &
                                                  + block(b) % dn       (i+1,j  ,k  ,1 ) &
                                                  * block(b) % edge_area(i+1,j  ,k  ,1 ) &
                                                  + block(b) % dn       (i  ,j+1,k  ,2 ) &
                                                  * block(b) % edge_area(i  ,j+1,k  ,2 ) &
                                                  + block(b) % dn       (i  ,j  ,k+1,3 ) &
                                                  * block(b) % edge_area(i  ,j  ,k+1,3 ) )&
                                                  * block(b) % cell_vol (i  ,j  ,k     ) &
                                                  * block(b) % a        (i  ,j  ,k  ,1 )


                     block(b) % mat_geo( 1,i,j,k) = - block(b) % dn       (i+1,j  ,k  ,1 ) &
                                                    * block(b) % edge_area(i+1,j  ,k  ,1 ) &
                                                    * block(b) % cell_vol (i  ,j  ,k     ) &
                                                    * block(b) % a        (i  ,j  ,k  ,1 )

                     block(b) % mat_geo( 2,i,j,k) = - block(b) % dn       (i  ,j+1,k  ,2 ) &
                                                    * block(b) % edge_area(i  ,j+1,k  ,2 ) &
                                                    * block(b) % cell_vol (i  ,j  ,k     ) &
                                                    * block(b) % a        (i  ,j  ,k  ,1 )

                     block(b) % mat_geo( 3,i,j,k) = - block(b) % dn       (i  ,j  ,k+1,3 ) &
                                                    * block(b) % edge_area(i  ,j  ,k+1,3 ) &
                                                    * block(b) % cell_vol (i  ,j  ,k     ) &
                                                    * block(b) % a        (i  ,j  ,k  ,1 )


                     if (k > 1) then
                        n_mat_element = n_mat_element + 1
                     end if

                     if (j > 1) then
                        n_mat_element = n_mat_element + 1
                     end if

                     if (i > 1) then
                        n_mat_element = n_mat_element + 1
                     end if

                     if (i < block(b) % nCell(1)) then
                        n_mat_element = n_mat_element + 1
                     end if

                     if (j < block(b) % nCell(2)) then
                        n_mat_element = n_mat_element + 1
                     end if

                     if (k < block(b) % nCell(3)) then
                        n_mat_element = n_mat_element + 1
                     end if

                     n_mat_element = n_mat_element + 1

                  end do
               end do
            end do
         end if
      end do
      if (implicit) then
         n_nonzero_elements = n_mat_element

         vec_size1 = n_mat_element 

         vec_size2 = n_mat_element * 2

         vec_size3 = nCell

         write(*,*) ncell,n_mat_element,dble(n_mat_element)/dble(ncell), vec_size1,vec_size2,vec_size3

         allocate ( jmat      (vec_size1)          )
         allocate ( imat      (nCell+1)            )
         allocate ( rhs       (nCell)              )
         allocate ( sol       (nCell)              )
         allocate ( pivot     (nCell)              )
         allocate ( ha        (vec_size3  ,dim2ha) )

         allocate (mat(vec_size1))
         allocate (mat2(vec_size1))
         
         n = 0
         pos = 0
         col = 0
         do b = 1, nBlock
            dci = 1
            dcj = block(b) % nCell(1)
            dck = block(b) % nCell(2) * block(b) % nCell(1)
            do k = 1, block(b) % nCell(3)
               do j = 1,block(b) % nCell(2)
                  do i = 1, block(b) % nCell(1)
                     col = col + 1
                     pos = pos + 1
                     n = n + 1
                     jmat(n) = pos
                     imat(col) = n
                     if (k > 1) then
                        n  = n + 1
                        jmat(n) = pos - dck
                     end if
                     if (j > 1) then
                        n = n + 1
                        jmat(n) = pos - dcj
                     end if
                     if (i > 1) then
                        n = n + 1
                        jmat(n) = pos - dci
                     end if
                     if (i < block(b) % nCell(1)) then
                        n = n + 1
                        jmat(n) = pos + dci
                     end if
                     if (j < block(b) % nCell(2)) then
                        n = n + 1
                        jmat(n) = pos + dcj
                     end if
                     if (k < block(b) % nCell(3)) then
                        n  = n + 1
                        jmat(n) = pos + dck
                     end if

                  end do
               end do
            end do
         end do
         imat(col+1) = n+1
!         row_pos_save = row_pos
!         col_pos_save = col_pos
!         col = 0
!         do n=1,n_nonzero_elements
!            if (n == imat(col+1) ) col = col + 1
!            write(*,*) mat(n),jmat(n),col
!         end do
!         stop
          
      end if


   end subroutine init

   function cross(a, b)
      implicit none
      real(kind=dp), dimension(3) :: cross
      real(kind=dp), dimension(3), INTENT(IN) :: a, b
      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
   end function cross

   function abs_cross(a, b)
      implicit none
      real(kind=dp), dimension(3) :: abs_cross
      real(kind=dp), dimension(3), INTENT(IN) :: a, b
      abs_cross(1) = abs(a(2) * b(3) - a(3) * b(2))
      abs_cross(2) = abs(a(3) * b(1) - a(1) * b(3))
      abs_cross(3) = abs(a(1) * b(2) - a(2) * b(1))
   end function abs_cross

   function len_cross(a, b)
      implicit none
      real(kind=dp) :: len_cross
      real(kind=dp), dimension(3), INTENT(IN) :: a, b
      len_cross = (a(2) * b(3) - a(3) * b(2)) * (a(2) * b(3) - a(3) * b(2))
      len_cross = len_cross + (a(3) * b(1) - a(1) * b(3)) * (a(3) * b(1) - a(1) * b(3))
      len_cross = len_cross + (a(1) * b(2) - a(2) * b(1)) * (a(1) * b(2) - a(2) * b(1))
      len_cross = sqrt(len_cross)
   end function len_cross

   function vec_len(vec)
      implicit none
      real(kind=dp) :: vec_len
      real(kind=dp), dimension(3), INTENT(IN) :: vec

      vec_len = sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))

   end function vec_len
end module data
