program gridgen

use cgnslib
use types
use const

USE cgns_types, ONLY: CGSIZE_T
implicit none

integer, parameter :: nBlock = 1

integer, parameter :: Dimen = 3

integer :: imax = 100

integer, parameter :: jmax = 2

integer, parameter :: kmax = 2

real(kind=dp), parameter :: x0  =  0.0E0_dp
real(kind=dp), parameter :: x1  =  1.0E0_dp

real(kind=dp), parameter :: y0  = -1.0E-1_dp
real(kind=dp), parameter :: y1  =  1.0E-1_dp

real(kind=dp), parameter :: z0  = -1.0E-1_dp
real(kind=dp), parameter :: z1  =  1.0E-1_dp

INTEGER, PARAMETER :: ioout = 10
integer, parameter :: Version = 1001

real(kind=dp), parameter :: pi = 2.E0_dp*DASIN(1.D0)

integer(kind = CGSIZE_T) :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol

real(kind=dp) :: x, temp_wall = 300.0E0_dp, hf_wall = 385.0E2_dp

real(kind=dp) :: pos_start(Dimen,nBlock),pos_end(Dimen,nBlock)

real(kind=dp) :: dn

type(tblock) :: block(nBlock)

integer :: i,j,k,b,d,ii,arg_count

integer :: n2ii(Dimen,Dimen)

integer(kind=CGSIZE_T) :: isize (Dimen,3)

character(len=100) :: arg



character(len=7) :: zonename
character(len=11),parameter :: coordnames(3) = (/"CoordinateX","CoordinateY","CoordinateZ"/)


arg_count=command_argument_count()
if (arg_count >= 1) then
   CALL get_command_argument(1, arg)
   read(arg,*) imax
   write(*,*) "IMAX = ",imax
end if

n2ii= 0
do d = 1,Dimen
   n2ii(d,d) = 1
end do

write(*,'(A)') "SIMPLE GRID GEN"

block(1) % nCell = (/ imax,jmax,kmax/)


pos_start(:,1) = (/x0, y0, z0/)

pos_end(:,1) = (/x1, y1, z1/)

do b = 1, nBlock
   block(b) % nPkt = block(b) % nCell+1
   allocate (block(b) % xyz(block(b) % nPkt(1),block(b) % nPkt(2),block(b) % nPkt(3),3))
   allocate (block(b) % T(block(b) % nCell(1),block(b) % nCell(2),block(b) % nCell(3),1))
   do d = 1,Dimen
      dn = (pos_end(d,b) - pos_start(d,b)) / block(b) % nCell(d)
      do k = 1, block(b) % nPkt(3)
         do j = 1, block(b) % nPkt(2)
            do i = 1, block(b) % nPkt(1)
               ii = n2ii(1,d) * i + n2ii(2,d) * j + n2ii(3,d) * k -1
               block(b) % xyz(i,j,k,d) = pos_start(d,b) + dn * dble(ii)
            end do
         end do
      end do
   end do
   dn = (pos_end(1,b) - pos_start(1,b)) / block(b) % nCell(1)
   do i = 1,block(b) % nCell(1)
      x = pos_start(1,b) + 0.5E0_dp * dn + dn * dble(i-1)
      block(b) % T(i,:,:,:) = 300.0D0 !dble(i)!sin((pi * x) / (pos_end(1,b) - pos_start(1,b)) )
   end do
end do

call cg_open_f("data_in.cgns",CG_MODE_WRITE,cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()
call cg_base_write_f(cgns_file,"Grid",Dimen,Dimen,cgns_base,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

do b = 1,nBlock
   write(zonename,'("Zone",I3)') b

   do d = 1, Dimen
      isize(d,1) = block(b) % nPkt(d)
      isize(d,2) = block(b) % nCell(d)
      isize(d,3) = 0
   end do
   call cg_zone_write_f(cgns_file,cgns_base,zonename,isize,Structured,cgns_zone,ierror)
   if (ierror /= CG_OK) call cg_error_exit_f()
   do d = 1, Dimen
      call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,coordnames(d) &
                           ,block(b) % xyz(:,:,:,d),cgns_coord,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
   end do
   call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,"0",CellCenter,cgns_sol,ierror)
   if (ierror /= CG_OK) call cg_error_exit_f()
   call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
                     ,"Temperature",block(b) % T,cgns_var,ierror)
end do

call cg_close_f(cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()

open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")
b = 1
write (ioout) Version,nBlock
!!!!! BLOCK 1
!!! WEST
!write(ioout) 1,BC_SYMMETRY,1,block(b) % nCell(3),1,block(b) % nCell(2)
write(ioout) 1,BC_HEATFLUX,1,block(b) % nCell(3),1,block(b) % nCell(2)
write(ioout) ((hf_wall,j = 1,block(b) % nCell(2)),k=1,block(b) % nCell(3))
!!! OST
write(ioout) 1,BC_ISOTHERMAL,1,block(b) % nCell(3),1,block(b) % nCell(2)
write(ioout) ((temp_wall,j = 1,block(b) % nCell(2)),k=1,block(b) % nCell(3))
!!! SÜD
write(ioout) 1,BC_SYMMETRY,1,block(b) % nCell(3),1,block(b) % nCell(1)
!!! NORD
write(ioout) 1,BC_SYMMETRY,1,block(b) % nCell(3),1,block(b) % nCell(1)
!!! BACK
write(ioout) 1,BC_SYMMETRY,1,block(b) % nCell(2),1,block(b) % nCell(1)
!!! FRONT
write(ioout) 1,BC_SYMMETRY,1,block(b) % nCell(2),1,block(b) % nCell(1)

!close (ioout)
!b = 1
!!WEST
!write(ioout) i,-3,i,kmax(b),1,jmax(b)
!!write(ioout) ((temp,j=1,jmax(1)),k=1,kmax(1))
!!OST
!write(ioout) i,bc(2),i,kmax(b),i,jmax(b)
!
!!SÜD
!
!write(ioout) i,bc(3),i,kmax(b),i,imax(b)
!
!!NORD
!k = 3
!j = 2
!write(ioout) k
!write(ioout) bc(4),i,kmax(b),i,5
!
!write(ioout) 2,i,kmax(b),6,10
!write(ioout) 0
!write(ioout) -5, 1, 0, 0   !I1
!write(ioout)  -jmax(1), 0, 1, 0   !i2
!write(ioout)  0, 0, 0, 1   !i3
!write(ioout) bc(4),i,kmax(b),11,15
!!BACK
!write(ioout) i,bc(5),i,jmax(b),i,imax(b)
!!FRONT
!write(ioout) i,bc(6),i,jmax(b),i,imax(b)
!
!b = 2
!!WEST
!write(ioout) i,bc(2),i,kmax(b),1,jmax(b)
!!OST
!write(ioout) i,bc(2),i,kmax(b),i,jmax(b)
!!SÜD
!!write(ioout) i,bc(3),i,kmax(b),i,imax(b)
!write(ioout) 1,1,i,kmax(b),i,imax(b)
!write(ioout) 0
!write(ioout)  5, 1, 0, 0   !I1
!write(ioout)  jmax(1), 0, 1, 0   !i2
!write(ioout)  0, 0, 0, 1   !i3
!
!!NORD
!write(ioout) i,bc(4),i,kmax(b),i,imax(b)
!!BACK
!write(ioout) i,bc(5),i,jmax(b),i,imax(b)
!!FRONT
!write(ioout) i,bc(6),i,jmax(b),i,imax(b)


write(*,'(A)') "done"

end program
