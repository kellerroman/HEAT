program gridgen

use cgnslib

USE cgns_types, ONLY: CGSIZE_T
implicit none

integer, parameter :: imax = 20
integer, parameter :: jmax = 5
integer, parameter :: kmax = 1

INTEGER, PARAMETER :: ioout = 10
integer, parameter :: Version = 1001
integer, parameter :: Dimen = 3
integer, parameter :: nBlock = 1

integer, parameter :: bc(6) = (/-1,-3,-3,-3,-3,-3/)

integer(kind = CGSIZE_T) :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol

real(kind=8) :: xyz (imax+1,jmax+1,kmax+1,Dimen)
real(kind=8) :: temp_ini (imax,jmax,kmax)

real(kind=8) :: temp = 500.0D0

integer :: i,j,k

integer(kind=CGSIZE_T) :: isize (dimen,3)

write(*,'(A)') "SIMPLE GRID GEN"

do i = 1,imax+1
   do j = 1,jmax+1
      do k = 1,kmax+1
         xyz(i,j,k,1) = 2.0D0/dble(imax) * dble(i-1)

         xyz(i,j,k,2) = 1.0D-1/dble(jmax) * dble(j-1)

         xyz(i,j,k,3) = 1.0D-2/dble(kmax) * dble(k-1)
      end do
   end do
end do
do i = 1,imax
   do j = 1,jmax
      do k = 1,kmax
         temp_ini(i,j,k) = 300.0D0 + 200.0D0 / dble(imax-1) * dble(imax-i)
      end do
   end do
end do
call cg_open_f("git.cgns",CG_MODE_WRITE,cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_base_write_f(cgns_file,"Grid",Dimen,Dimen,cgns_base,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

isize(1,1) = imax+1
isize(2,1) = jmax+1
isize(3,1) = kmax+1

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

call cg_close_f(cgns_file,ierror)

if (ierror /= CG_OK) call cg_error_exit_f()


call cg_open_f("restart.cgns",CG_MODE_WRITE,cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_base_write_f(cgns_file,"SOLUTION",Dimen,Dimen,cgns_base,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_zone_write_f(cgns_file,cgns_base,"BLOCK 1",isize,Structured,cgns_zone,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,"0",CellCenter,cgns_sol,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
                  ,"Temperature",temp_ini,cgns_var,ierror)
call cg_close_f(cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()


open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")

i = 1
write (ioout) Version,nBlock
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
