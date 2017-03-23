program gridgen

   use mod_gridgen
   use const
implicit none

integer, parameter :: nBlock = 2

integer, parameter :: Dimen = 3

integer, parameter :: imax = 290!580

integer, parameter :: j1 = 38*2
integer, parameter :: j2 = 24
integer, parameter :: jmax = j1+j2

integer, parameter :: k1 = 73*2
integer, parameter :: k2 = 24
integer, parameter :: kmax = k1+k2
!  UNTERE LINKE HÄLFTE DER BRENNKAMMER
!
!
!______                         y2
!|    |
!|    |
!| 2  |                            j2
!|    |
!|____|_________________________y1
!|                             |
!|             1               |   j1
!|_____________________________|y0
!z0   z1                      z2
!  k1              k2
!

real(kind=dp), parameter :: x0  =   0.0E-3_dp
real(kind=dp), parameter :: x1  = 290.0E-3_dp

real(kind=dp), parameter :: y0  = -25.0E-3_dp
real(kind=dp), parameter :: y1  =  -6.0E-3_dp
real(kind=dp), parameter :: y2  =   0.0E-3_dp

real(kind=dp), parameter :: z0  = -42.5E-3_dp
real(kind=dp), parameter :: z1  =  -6.0E-3_dp
real(kind=dp), parameter :: z2  =   0.0E-3_dp

INTEGER, PARAMETER :: ioout = 10
integer, parameter :: Version = 1001


real(kind=dp) :: temp_ini = 300.0D0

real(kind=dp) :: pos_start(Dimen,nBlock),pos_end(Dimen,nBlock)

real(kind=dp) :: hf(imax,j2)

real(kind=dp) :: dn

integer :: i,j,k,b,d,ii

integer :: n2ii(Dimen,Dimen)


n2ii= 0
do d = 1,Dimen
   n2ii(d,d) = 1
end do

write(*,'(A)') "SIMPLE GRID GEN"

call add_block(imax,j1,kmax)
call add_block(imax,j2,k1)

call allocate_blocks

pos_start(:,1) = (/x0, y0, z0/)
pos_start(:,2) = (/x0, y1, z0/)

pos_end(:,1) = (/x1, y1, z2/)
pos_end(:,2) = (/x1, y2, z1/)

do i = 1,imax
   do j = 1,j2
      hf(i,j) = 1+4/dble(imax-1)*dble(i-1)
      hf(i,j) = hf(i,j) * (0.5D0+0.5D0/ dble(j2) * dble(j))
   end do
!   write(*,*) hf(i,:)
end do
hf = hf * -1E6_dp
do b = 1, nBlock
   blocks(b) % temps = 300.0D0 !+ 100.0D0 * dble(b)
   do d = 1,3
      dn = (pos_end(d,b) - pos_start(d,b)) / blocks(b) % nCells(d)
      do k = 1, blocks(b) % nPkts(3)
         do j = 1, blocks(b) % nPkts(2)
            do i = 1, blocks(b) % nPkts(1)
               ii = n2ii(1,d) * i + n2ii(2,d) * j + n2ii(3,d) * k -1
               blocks(b) % xyzs(i,j,k,d) = pos_start(d,b) + dn * dble(ii)
            end do
         end do
      end do
   end do
end do


call write_grid()
call write_xdmf()
open (ioout,file="bc.bin",form="unformatted",access="stream",status="replace")
b = 1
write (ioout) Version,nBlock


do  b = 1,nBlock
   !!! WEST
   write(ioout) 1,BC_SYMMETRY,1,blocks(b) % nCells(3),1,blocks(b) % nCells(2)
   !!! OST
   write(ioout) 1,BC_SYMMETRY,1,blocks(b) % nCells(3),1,blocks(b) % nCells(2)
   !!! SÜD
   if (b == 1) then
      write(ioout) 1,BC_SYMMETRY,1,blocks(b) % nCells(3),1,blocks(b) % nCells(1)
   else if (b == 2) then
      !!! CONNECTION TO BLOCK 1
      write(ioout) 1 !NUMBER OF BOUNDARY CONDITIONS
      write(ioout) 1 !BLOCK NUMBER
      write(ioout) 1,blocks(b) % nCells(3),1,blocks(b) % nCells(1)
      write(ioout) 1 ! CPU_ID
      write(ioout) 0,1,0,0
      write(ioout) blocks(1) % nCells(2),0,1,0
      write(ioout) 0,0,0,1
   else
      write(*,*) "BLOCK SÜD RB NOT DEFINDED",b
      stop 1
   end if
   !!! NORD
   if (b == 1) then
      write(ioout) 2 !NUMBER OF BOUNDARY CONDITIONS

      !!! CONNECTION TO BLOCK 2
      write(ioout)  2 !BLOCK NUMBER
      write(ioout)  1,blocks(2) % nCells(3),1,blocks(b) % nCells(1)
      write(ioout) 1 ! CPU_ID
      write(ioout) 0,1,0,0
      write(ioout) -blocks(1) % nCells(2),0,1,0
      write(ioout) 0,0,0,1

      !!!!!! HEATFLUX
      write(ioout) BC_HEATFLUX,blocks(2) % nCells(3)+1,blocks(1) % nCells(3),1,blocks(b) % nCells(1)
      write(ioout) ((hf(i,j),i=1,imax),j=1,j2)

      !!!!!! SYMMETRIE
      !write(ioout) BC_SYMMETRY,block(2) % nCell(3)+1,block(1) % nCell(3),1,block(b) % nCell(1)
   else if (b == 2) then
      write(ioout) 1,BC_SYMMETRY,1,blocks(b) % nCells(3),1,blocks(b) % nCells(1)
   else
      write(*,*) "BLOCK NORD RB NOT DEFINDED",b
      stop 1
   end if
   !!! BACK
   write(ioout) 1,BC_SYMMETRY,1,blocks(b) % nCells(2),1,blocks(b) % nCells(1)
   !!! FRONT
   if (b == 1) then
      write(ioout) 1,BC_SYMMETRY,1,blocks(b) % nCells(2),1,blocks(b) % nCells(1)
   else if (b == 2) then
      !!!!!! HEATFLUX
      write(ioout) 1,BC_HEATFLUX,1,blocks(b) % nCells(2),1,blocks(b) % nCells(1)
      write(ioout) ((hf(i,j),i=1,imax),j=1,j2)
   else
      write(*,*) "BLOCK FRONT RB NOT DEFINDED",b
      stop 1
   end if
end do
close (ioout)

write(*,'(A)') "done"

end program
