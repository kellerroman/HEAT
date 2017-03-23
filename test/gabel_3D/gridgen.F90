program gridgen

   use mod_gridgen
   use const
implicit none


integer, parameter :: NBLOCK       = 3
integer, parameter :: imax(NBLOCK) = (/17,4,5/)
integer, parameter :: jmax(NBLOCK) = (/3,5,4/)
integer, parameter :: kmax(NBLOCK) = (/8,8,8/)
integer, parameter :: di_2         = 3
!< Verschiebung der Blockansatzes in i-Richtung fuer block 2 relativ zu block 1
integer, parameter :: di_3         = 10
!< Verschiebung der Blockansatzes in i-Richtung fuer block 3 relativ zu block 1

INTEGER, PARAMETER :: ioout = 10
integer, parameter :: Version = 1001

integer, parameter :: bc(6) = (/BC_SYMMETRY,BC_SYMMETRY,BC_SYMMETRY,BC_SYMMETRY,BC_SYMMETRY,BC_SYMMETRY/)
!integer, parameter :: bc(6) = (/BC_ISOTHERMAL,BC_SYMMETRY,BC_SYMMETRY,BC_SYMMETRY,BC_SYMMETRY,BC_SYMMETRY/)


real(kind=8) :: temp = 300.0D0

integer :: i,j,k,b
integer :: b2c   !Block to Connect

write(*,'(A)') "SIMPLE GRID GEN"

do b = 1, NBLOCK
   call add_block(imax(b),jmax(b),kmax(b))
end do
call allocate_blocks
b = 1
do i = 1,imax(1)+1
   do j = 1,jmax(1)+1
      do k = 1,kmax(1)+1
         blocks(b) % xyzs(i,j,k,1) = 1.0D-1 * dble(i-1)

         blocks(b) % xyzs(i,j,k,2) = 1.0D-1 * dble(j-1)

         blocks(b) % xyzs(i,j,k,3) = 1.0D-1 * dble(k-1)!+1.0D-1 * dble(i-1)
      end do
   end do
end do

b = 2
blocks(b) % xyzs(:,1,:,:) = blocks(1) % xyzs(di_2:di_2+imax(b)+1,jmax(1)+1,:,:)
do i = 1,imax(b)+1
   do j = 2,jmax(b)+1
      do k = 1,kmax(b)+1
         blocks(b) % xyzs(i,j,k,1) = blocks(b) % xyzs(i,1,k,1)

         blocks(b) % xyzs(i,j,k,2) = blocks(b) % xyzs(i,1,k,2) + 1.0D-1 * dble(j-1)

         blocks(b) % xyzs(i,j,k,3) = blocks(b) % xyzs(i,1,k,3)
      end do
   end do
end do

b = 3
blocks(b) % xyzs(:,1,:,:) = blocks(1) % xyzs(di_3:di_3+imax(b)+1,jmax(1)+1,:,:)
do i = 1,imax(b)+1
   do j = 2,jmax(b)+1
      do k = 1,kmax(b)+1
         blocks(b) % xyzs(i,j,k,1) = blocks(b) % xyzs(i,1,k,1)

         blocks(b) % xyzs(i,j,k,2) = blocks(b) % xyzs(i,1,k,2) + 1.0D-1 * dble(j-1)

         blocks(b) % xyzs(i,j,k,3) = blocks(b) % xyzs(i,1,k,3)
      end do
   end do
end do

do b = 1,NBLOCK
   do i = 1,imax(b)
      do j = 1,jmax(b)
         do k = 1,kmax(b)
            if (b == 1) then
               blocks(b) % temps(i,j,k) = 1000.0D0 !dble(i)
            elseif (b == 2) then
               blocks(b) % temps(i,j,k) = 1100.0D0 !dble(j)
            else
               blocks(b) % temps(i,j,k) = 1200.0D0 !dble(k)
            end if
         end do
      end do
   end do
end do

call write_grid()
call write_xdmf()

open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")

i = 1
write (ioout) Version,nBlock
b = 1
!WEST
write(ioout) i,bc(1),i,kmax(b),1,jmax(b)
if (bc(1) == BC_ISOTHERMAL) then
   write(ioout) ((temp,j=1,jmax(b)),k=b,kmax(b))
end if

!OST
write(ioout) i,bc(2),i,kmax(b),i,jmax(b)

!SÜD
write(ioout) i,bc(3),i,kmax(b),i,imax(b)

!NORD
k = 5
j = 2
! Noedseite hat 5 verschieden Randbedingungen
write(ioout) k
! Erste Randbedingung bis zum Anschluss von B2 
write(ioout) bc(4),i,kmax(b),i,di_2-1
! Anschluss von Block 2
! Ueber die gesammt k dimesnion und in i-Richtung entlang der blocklaenge
b2c = 2
write(ioout) b2c,i,kmax(b),di_2,di_2+imax(b2c)
! Block 2 ist auch auf CPU 0
write(ioout) 0
! Index Verschiebungsmatrix
write(ioout) -di_2+1, 1, 0, 0   !I1
write(ioout)  -jmax(1), 0, 1, 0   !i2
write(ioout)  0, 0, 0, 1   !i3
! Dritte Randbedinugng zwischen den Bloecken
write(ioout) bc(4),i,kmax(b),di_2+imax(b2c)+1,di_3-1
! Anschluss von Block 3
! Ueber die gesammt k dimesnion und in i-Richtung entlang der blocklaenge
b2c = 3
write(ioout) b2c,i,kmax(b),di_3,di_3+imax(b2c)
! Block 3 ist auch auf CPU 0
write(ioout) 0
! Index Verschiebungsmatrix
write(ioout) -di_3+1   , 1, 0, 0   !I1
write(ioout)  -jmax(1) , 0, 1, 0   !i2
write(ioout)  0        , 0, 0, 1   !i3
! Fuenfte Randbedinugng zwischen den Bloecken
write(ioout) bc(4),i,kmax(b),di_3+imax(b2c)+1,imax(b)

!BACK
write(ioout) i,bc(5),i,jmax(b),i,imax(b)

!FRONT
write(ioout) i,bc(6),i,jmax(b),i,imax(b)

b = 2
!WEST
write(ioout) i,bc(1),i,kmax(b),i,jmax(b)
if (bc(1) == BC_ISOTHERMAL) then
   write(ioout) ((temp,j=1,jmax(b)),k=1,kmax(b))
end if
!OST
write(ioout) i,bc(2),i,kmax(b),i,jmax(b)
!SUED
!write(ioout) i,bc(3),i,kmax(b),i,imax(b)
write(ioout) 1,1,i,kmax(b),i,imax(b)
write(ioout) 0
write(ioout)  di_2-1  , 1, 0, 0   !I1
write(ioout)  jmax(1) , 0, 1, 0   !i2
write(ioout)  0       , 0, 0, 1   !i3

!NORD
write(ioout) i,bc(4),i,kmax(b),i,imax(b)
!BACK
write(ioout) i,bc(5),i,jmax(b),i,imax(b)
!FRONT
write(ioout) i,bc(6),i,jmax(b),i,imax(b)

b = 3
!WEST
write(ioout) i,bc(1),i,kmax(b),i,jmax(b)
if (bc(1) == BC_ISOTHERMAL) then
   write(ioout) ((temp,j=1,jmax(b)),k=1,kmax(b))
end if
!OST
write(ioout) i,bc(2),i,kmax(b),i,jmax(b)
!SÜD
!write(ioout) i,bc(3),i,kmax(b),i,imax(b)
write(ioout) 1,1,i,kmax(b),i,imax(b)
write(ioout) 0
write(ioout)  di_3-1  , 1, 0, 0   !I1
write(ioout)  jmax(1) , 0, 1, 0   !i2
write(ioout)  0       , 0, 0, 1   !i3
!NORD
write(ioout) i,bc(4),i,kmax(b),i,imax(b)
!BACK
write(ioout) i,bc(5),i,jmax(b),i,imax(b)
!FRONT
write(ioout) i,bc(6),i,jmax(b),i,imax(b)

close (ioout)



write(*,'(A)') "done"

end program
