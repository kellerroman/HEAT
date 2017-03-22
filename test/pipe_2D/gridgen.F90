program gridgen

   use mod_gridgen
implicit none

integer, parameter :: imax = 20
integer, parameter :: jmax = 5
integer, parameter :: kmax = 2

integer, parameter :: bc(6) = (/-1,-3,-3,-3,-3,-3/)
integer, parameter :: Version = 1001
integer            :: ioout
real(kind=8) :: temp = 500.0D0

integer :: i,j,k


write(*,'(A)') "SIMPLE GRID GEN"

call add_block(imax-1,jmax-1,kmax-1)
call allocate_blocks

do i = 1,imax
   do j = 1,jmax
      do k = 1,kmax
         blocks(1) % xyzs(i,j,k,1) = 2.0D0/dble(imax) * dble(i-1)

         blocks(1) % xyzs(i,j,k,2) = 1.0D-1/dble(jmax) * dble(j-1)

         blocks(1) % xyzs(i,j,k,3) = 1.0D-2/dble(kmax) * dble(k-1)
      end do
   end do
end do
do i = 1,imax-1
   do j = 1,jmax-1
      do k = 1,kmax-1
         blocks(1) % temps(i,j,k) = 300.0D0 + 200.0D0 / dble(imax-1) * dble(imax-i)
      end do
   end do
end do

call write_grid()
call write_xdmf()

open (newunit=ioout,file="bc.bin",form="unformatted",access="stream",status="replace")

i = 1
write (ioout) Version,1
!WEST
write(ioout) i,bc(1),i,kmax,1,jmax
if (bc(1) == -1) write(ioout) ((temp,j=1,jmax),k=1,kmax)
!OST
write(ioout) i,bc(2),i,kmax,i,jmax
!SÜD
write(ioout) i,bc(3),i,kmax,i,imax
!NORD
write(ioout) i,bc(4),i,kmax,i,imax
!BACK
write(ioout) i,bc(5),i,jmax,i,imax
!FRONT
write(ioout) i,bc(6),i,jmax,i,imax

close (ioout)



write(*,'(A)') "done"

end program
