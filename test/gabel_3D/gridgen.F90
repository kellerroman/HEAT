program gridgen

   use mod_gridgen
implicit none


integer, parameter :: NBLOCK = 3
integer, parameter :: imax(NBLOCK) = (/17,4,5/)
integer, parameter :: jmax(NBLOCK) = (/3,5,4/)
integer, parameter :: kmax(NBLOCK) = (/8,3,2/)

INTEGER, PARAMETER :: ioout = 10
integer, parameter :: Version = 1001

integer, parameter :: bc(6) = (/-1,-3,-3,-3,-3,-3/)


real(kind=8) :: temp = 300.0D0

integer :: i,j,k,b

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

         blocks(b) % xyzs(i,j,k,3) = 1.0D-1 * dble(k-1)
      end do
   end do
end do

b = 2
do i = 1,imax(b)+1
   do j = 1,jmax(b)+1
      do k = 1,kmax(b)+1
         blocks(b) % xyzs(i,j,k,1) = blocks(1) % xyzs(3,jmax(1)+1,k,1) + 1.0D-1 * dble(i-1)

         blocks(b) % xyzs(i,j,k,2) = blocks(1) % xyzs(i+1,jmax(1)+1,k,2) + 1.0D-1 * dble(j-1)

         blocks(b) % xyzs(i,j,k,3) = blocks(1) % xyzs(i,jmax(1)+1,3,3) + 1.0D-1 * dble(k-1)
      end do
   end do
end do

b = 3
do i = 1,imax(b)+1
   do j = 1,jmax(b)+1
      do k = 1,kmax(b)+1
         blocks(b) % xyzs(i,j,k,1) = blocks(1) % xyzs(10,jmax(1)+1,k,1) + 1.0D-1 * dble(i-1)

         blocks(b) % xyzs(i,j,k,2) = blocks(1) % xyzs(i+7,jmax(1)+1,k,2) + 1.0D-1 * dble(j-1)

         blocks(b) % xyzs(i,j,k,3) = blocks(1) % xyzs(i,jmax(1)+1,4,3) + 1.0D-1 * dble(k-1)
      end do
   end do
end do

do b = 1,NBLOCK
   do i = 1,imax(b)
      do j = 1,jmax(b)
         do k = 1,kmax(b)
            blocks(b) % temps(i,j,k) = dble(i)
            blocks(b) % temps(i,j,k) = dble(j)
            blocks(b) % temps(i,j,k) = dble(k)
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
write(ioout) i,-3,i,kmax(b),1,jmax(b)
!write(ioout) ((temp,j=1,jmax(1)),k=1,kmax(1))
!OST
write(ioout) i,bc(2),i,kmax(b),i,jmax(b)

!SÜD

write(ioout) i,bc(3),i,kmax(b),i,imax(b)

!NORD
k = 3
j = 2
write(ioout) k
write(ioout) bc(4),i,kmax(b),i,5

write(ioout) 2,i,kmax(b),6,10
write(ioout) 0
write(ioout) -5, 1, 0, 0   !I1
write(ioout)  -jmax(1), 0, 1, 0   !i2
write(ioout)  0, 0, 0, 1   !i3
write(ioout) bc(4),i,kmax(b),11,15
!BACK
write(ioout) i,bc(5),i,jmax(b),i,imax(b)
!FRONT
write(ioout) i,bc(6),i,jmax(b),i,imax(b)

b = 2
!WEST
write(ioout) i,bc(2),i,kmax(b),1,jmax(b)
!OST
write(ioout) i,bc(2),i,kmax(b),i,jmax(b)
!SÜD
!write(ioout) i,bc(3),i,kmax(b),i,imax(b)
write(ioout) 1,1,i,kmax(b),i,imax(b)
write(ioout) 0
write(ioout)  5, 1, 0, 0   !I1
write(ioout)  jmax(1), 0, 1, 0   !i2
write(ioout)  0, 0, 0, 1   !i3

!NORD
write(ioout) i,bc(4),i,kmax(b),i,imax(b)
!BACK
write(ioout) i,bc(5),i,jmax(b),i,imax(b)
!FRONT
write(ioout) i,bc(6),i,jmax(b),i,imax(b)
close (ioout)



write(*,'(A)') "done"

end program
