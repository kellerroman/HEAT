program gridgen

use cgnslib
use types

USE cgns_types, ONLY: CGSIZE_T
implicit none

integer, parameter :: nBlock = 4

integer, parameter :: Dimen = 3

integer, parameter :: imax = 290!580

integer, parameter :: j1 = 19!38
integer, parameter :: j2 = 12!24
integer, parameter :: j3 = 19!38
integer, parameter :: jmax = j1+j2+j3

integer, parameter :: k1 = 73
integer, parameter :: k2 = 24
integer, parameter :: k3 = 73
integer, parameter :: kmax = k1+k2+k3
! _____________________________ y3
!|                             |  j3
!|             4               |
!|_____________________________|y2
!|    |                  |     |
!|    |                  |     |
!| 2  |                  |  3  |   j2
!|    |                  |     |
!|____|__________________|_____|y1
!|                             |
!|             1               |   j1
!|_____________________________|y0
!z0   z1                 z2   z3
!  k1       k2             k3
!

real(kind=dp), parameter :: x0  =   0.0E-3_dp
real(kind=dp), parameter :: x1  = 290.0E-3_dp

real(kind=dp), parameter :: y0  = -25.0E-3_dp
real(kind=dp), parameter :: y1  =  -6.0E-3_dp
real(kind=dp), parameter :: y2  =   6.0E-3_dp
real(kind=dp), parameter :: y3  =  25.0E-3_dp

real(kind=dp), parameter :: z0  = -42.5E-3_dp
real(kind=dp), parameter :: z1  =  -6.0E-3_dp
real(kind=dp), parameter :: z2  =   6.0E-3_dp
real(kind=dp), parameter :: z3  =  42.5E-3_dp

INTEGER, PARAMETER :: ioout = 10
integer, parameter :: Version = 1001

integer(kind = CGSIZE_T) :: ierror,cgns_file,cgns_base,cgns_zone,cgns_coord,cgns_var,cgns_sol

real(kind=dp) :: temp_ini = 300.0D0

real(kind=dp) :: pos_start(Dimen,nBlock),pos_end(Dimen,nBlock)

real(kind=dp) :: dn

type(tblock) :: block(nBlock)

integer :: i,j,k,b,d,ii

integer :: n2ii(Dimen,Dimen)

integer(kind=CGSIZE_T) :: isize (Dimen,3)

character(len=7) :: zonename
character(len=11),parameter :: coordnames(3) = (/"CoordinateX","CoordinateY","CoordinateZ"/)

n2ii= 0
do d = 1,Dimen
   n2ii(d,d) = 1
end do

write(*,'(A)') "SIMPLE GRID GEN"

block(1) % nCell = (/ imax,j1,kmax/)
block(2) % nCell = (/ imax,j2,k1/)
block(3) % nCell = (/ imax,j2,k3/)
block(4) % nCell = (/ imax,j3,kmax/)


pos_start(:,1) = (/x0, y0, z0/)
pos_start(:,2) = (/x0, y1, z0/)
pos_start(:,3) = (/x0, y1, z2/)
pos_start(:,4) = (/x0, y2, z0/)

pos_end(:,1) = (/x1, y1, z3/)
pos_end(:,2) = (/x1, y2, z1/)
pos_end(:,3) = (/x1, y2, z3/)
pos_end(:,4) = (/x1, y3, z3/)

do b = 1, nBlock
   block(b) % nPkt = block(b) % nCell+1
   allocate (block(b) % xyz(block(b) % nPkt(1),block(b) % nPkt(2),block(b) % nPkt(3),3))
   do d = 1,3
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
end do

call cg_open_f("git.cgns",CG_MODE_WRITE,cgns_file,ierror)
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
end do

call cg_close_f(cgns_file,ierror)
if (ierror /= CG_OK) call cg_error_exit_f()


!call cg_open_f("restart.cgns",CG_MODE_WRITE,cgns_file,ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()
!
!call cg_base_write_f(cgns_file,"SOLUTION",Dimen,Dimen,cgns_base,ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()
!isize(1,1) = imax(1)+1
!isize(2,1) = jmax(1)+1
!isize(3,1) = kmax(1)+1
!
!isize(1,2) = isize(1,1) - 1
!isize(2,2) = isize(2,1) - 1
!isize(3,2) = isize(3,1) - 1
!
!isize(1,3) = 0
!isize(2,3) = 0
!isize(3,3) = 0
!call cg_zone_write_f(cgns_file,cgns_base,"BLOCK 1",isize,Structured,cgns_zone,ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()
!
!!call cg_goto_f(cgns_file,cgns_base,ierror,"BLOCK 1",0,"end")
!!if (ierror /= CG_OK) call cg_error_exit_f()
!!call cg_link_write_f("GridCoordinates","git.cgns","Grid/BLOCK 1/GridCoordinates",ierror)
!!if (ierror /= CG_OK) call cg_error_exit_f()
!
!call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,"0",CellCenter,cgns_sol,ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()
!call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
!                  ,"Temperature",temp_ini,cgns_var,ierror)
!
!isize(1,1) = imax(2)+1
!isize(2,1) = jmax(2)+1
!isize(3,1) = kmax(2)+1
!
!isize(1,2) = isize(1,1) - 1
!isize(2,2) = isize(2,1) - 1
!isize(3,2) = isize(3,1) - 1
!
!isize(1,3) = 0
!isize(2,3) = 0
!isize(3,3) = 0
!call cg_zone_write_f(cgns_file,cgns_base,"BLOCK 2",isize,Structured,cgns_zone,ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()
!
!!call cg_goto_f(cgns_file,cgns_base,ierror,"BLOCK 2",0,"end")
!!if (ierror /= CG_OK) call cg_error_exit_f()
!!call cg_link_write_f("GridCoordinates","git.cgns","Grid/BLOCK 2/GridCoordinates",ierror)
!!if (ierror /= CG_OK) call cg_error_exit_f()
!
!call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,"0",CellCenter,cgns_sol,ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()
!call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
!                  ,"Temperature",temp_ini1,cgns_var,ierror)
!call cg_close_f(cgns_file,ierror)
!if (ierror /= CG_OK) call cg_error_exit_f()
!
!
open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")
b = 1
write (ioout) Version,nBlock
!!!!! BLOCK 1
!!! WEST
write(ioout) 1,-3,1,block(b) % nCell(3),1,block(b) % nCell(2)
!!! OST
write(ioout) 1,-3,1,block(b) % nCell(3),1,block(b) % nCell(2)

!!! SÜD
write(ioout) 1,-3,1,block(b) % nCell(3),1,block(b) % nCell(1)
!!! NORD
write(ioout) 3
write(ioout)  1,1    ,k1                ,1,block(b) % nCell(1)
write(ioout)  0,1,0,0   !I1
write(ioout)
write(ioout) -3,k1+1,k2                ,1,block(b) % nCell(1)
write(ioout)  1,k2+1,block(b) % nCell(3),1,block(b) % nCell(1)
!!! BACK
write(ioout) 1,-3,1,block(b) % nCell(2),1,block(b) % nCell(1)
!!! FRONT
write(ioout) 1,-3,1,block(b) % nCell(2),1,block(b) % nCell(1)

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
