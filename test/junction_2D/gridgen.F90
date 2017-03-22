program gridgen

use cgnslib
use types

USE cgns_types, ONLY: CGSIZE_T
implicit none

integer, parameter :: nBlock = 2

integer, parameter :: imax(nBlock) = (/15,5/)
integer, parameter :: jmax(nBlock) = (/5,10/)
integer, parameter :: kmax(nBlock) = (/5,3/)

INTEGER, PARAMETER :: ioout = 10
integer, parameter :: Version = 1001
integer, parameter :: Dimen = 3

integer, parameter :: bc(6) = (/-1,-3,-3,-3,-3,-3/)

integer(kind = CGSIZE_T) :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol

real(kind=8) :: xyz (imax(1)+1,jmax(1)+1,kmax(1)+1,Dimen)
real(kind=8) :: xyz1 (imax(2)+1,jmax(2)+1,kmax(2)+1,Dimen)
real(kind=8) :: temp_ini (imax(1),jmax(1),kmax(1))
real(kind=8) :: temp_ini1 (imax(2),jmax(2),kmax(2))

real(kind=8) :: temp = 300.0D0

type(tblock) :: block(nBlock)

integer :: i,j,k,b

integer(kind=CGSIZE_T) :: isize (dimen,3)

write(*,'(A)') "SIMPLE GRID GEN"


do i = 1,imax(1)+1
   do j = 1,jmax(1)+1
      do k = 1,kmax(1)+1
         xyz(i,j,k,1) = 1.5D0/dble(imax(1)) * dble(i-1)

         xyz(i,j,k,2) = 5.0D-1/dble(jmax(1)) * dble(j-1)

         xyz(i,j,k,3) = 5.0D-1/dble(kmax(1)) * dble(k-1)
      end do
   end do
end do

do i = 1,imax(2)+1
   do j = 1,jmax(2)+1
      do k = 1,kmax(2)+1
         xyz1(i,j,k,1) = xyz(6,jmax(1)+1,k,1) + 5.0D-1/dble(imax(2)) * dble(i-1)

         xyz1(i,j,k,2) = xyz(i+5,jmax(1)+1,k,2) + 1.0D-0/dble(jmax(2)) * dble(j-1)

         xyz1(i,j,k,3) = xyz(i,jmax(1)+1,2,3) + 3.0D-1/dble(kmax(2)) * dble(k-1)
      end do
   end do
end do

do i = 1,imax(1)
   do j = 1,jmax(1)
      do k = 1,kmax(1)
         temp_ini(i,j,k) = 100.0D0 !+ 200.0D0 / dble(imax(1)-1) * dble(imax(1)-i)
      end do
   end do
end do

temp_ini1 = 200.0D0

call cg_open_f("git.cgns",CG_MODE_WRITE,cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_base_write_f(cgns_file,"Grid",Dimen,Dimen,cgns_base,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

isize(1,1) = imax(1)+1
isize(2,1) = jmax(1)+1
isize(3,1) = kmax(1)+1

isize(1,2) = isize(1,1) - 1
isize(2,2) = isize(2,1) - 1
isize(3,2) = isize(3,1) - 1

isize(1,3) = 0
isize(2,3) = 0
isize(3,3) = 0

call cg_zone_write_f(cgns_file,cgns_base,"BLOCK 1",isize,Structured,cgns_zone,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateX",xyz(:,:,:,1),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateY",xyz(:,:,:,2),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateZ",xyz(:,:,:,3),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

isize(1,1) = imax(2)+1
isize(2,1) = jmax(2)+1
isize(3,1) = kmax(2)+1

isize(1,2) = isize(1,1) - 1
isize(2,2) = isize(2,1) - 1
isize(3,2) = isize(3,1) - 1

isize(1,3) = 0
isize(2,3) = 0
isize(3,3) = 0

call cg_zone_write_f(cgns_file,cgns_base,"BLOCK 2",isize,Structured,cgns_zone,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateX",xyz1(:,:,:,1),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateY",xyz1(:,:,:,2),cgns_coord,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,"CoordinateZ",xyz1(:,:,:,3),cgns_coord,ierror)

if (ierror /= CG_OK) call cg_error_exit_f()
call cg_close_f(cgns_file,ierror)

if (ierror /= CG_OK) call cg_error_exit_f()


call cg_open_f("restart.cgns",CG_MODE_WRITE,cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_base_write_f(cgns_file,"SOLUTION",Dimen,Dimen,cgns_base,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
isize(1,1) = imax(1)+1
isize(2,1) = jmax(1)+1
isize(3,1) = kmax(1)+1

isize(1,2) = isize(1,1) - 1
isize(2,2) = isize(2,1) - 1
isize(3,2) = isize(3,1) - 1

isize(1,3) = 0
isize(2,3) = 0
isize(3,3) = 0
call cg_zone_write_f(cgns_file,cgns_base,"BLOCK 1",isize,Structured,cgns_zone,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

!call cg_goto_f(cgns_file,cgns_base,ierror,"BLOCK 1",0,"end")
!if (ierror /= CG_OK) call cg_error_exit_f()
!call cg_link_write_f("GridCoordinates","git.cgns","Grid/BLOCK 1/GridCoordinates",ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()

call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,"0",CellCenter,cgns_sol,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
                  ,"Temperature",temp_ini,cgns_var,ierror)

isize(1,1) = imax(2)+1
isize(2,1) = jmax(2)+1
isize(3,1) = kmax(2)+1

isize(1,2) = isize(1,1) - 1
isize(2,2) = isize(2,1) - 1
isize(3,2) = isize(3,1) - 1

isize(1,3) = 0
isize(2,3) = 0
isize(3,3) = 0
call cg_zone_write_f(cgns_file,cgns_base,"BLOCK 2",isize,Structured,cgns_zone,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

!call cg_goto_f(cgns_file,cgns_base,ierror,"BLOCK 2",0,"end")
!if (ierror /= CG_OK) call cg_error_exit_f()
!call cg_link_write_f("GridCoordinates","git.cgns","Grid/BLOCK 2/GridCoordinates",ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()

call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,"0",CellCenter,cgns_sol,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
                  ,"Temperature",temp_ini1,cgns_var,ierror)
call cg_close_f(cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()


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
