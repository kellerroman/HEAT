module io
   use cgnslib
   use cgns_types, ONLY: CGSIZE_T
   use const
   use control, only: Dimen,nFace,nCorner,n_BC_Cell
   implicit none
   character(len=32),allocatable :: cgns_git_zonename(:)
   character(len=32) :: cgns_git_basename
   character(len=32) :: cgns_basename
   logical :: write_sol_header = .true.

contains
   subroutine error_out (text,file,line)
      implicit none
      character ( len = *) ,intent(in) :: text
      character ( len = *) ,intent(in) :: file
      integer              ,intent(in) :: line

      write(*,'(2(/100("!")))')
      write(*,'(10("!"),2X,A,1X,A,1X,A,1X,I0)') "ERROR IN",file(8:),"@",line
      write(*,'(10("!"),2X,A)') text
      write(*,'(100("!")/100("!"))')
      stop 1
   end subroutine error_out

   subroutine read_git()
      use data, only : nCell, nBlock, block
      use control, only: file_git_in
      implicit none

      integer(kind=CGSIZE_T) :: b ,d
      logical :: fexists
      integer(kind=CGSIZE_T) :: cgns_file,ierror,cgns_base,PhysDim,cgns_zone,zonetype
      integer(kind=CGSIZE_T),allocatable :: isize(:,:),istart(:)
      character(len=32),parameter :: coord_name(3) = (/ "CoordinateX","CoordinateY","CoordinateZ" /)
      real(kind=dp), allocatable :: temp_coord(:,:,:)

      nCell = 0

      inquire(file=trim(file_git_in),exist=fexists)

      if(.not. fexists) then
        call error_out("GITTER INPUT DATEI konnte nicht gefunden werden: "//TRIM(file_git_in),__FILE__,__LINE__)
      end if


      call cg_open_f(trim(file_git_in),CG_MODE_READ,cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      call cg_nbases_f(cgns_file,cgns_base,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      if (cgns_base /= 1) then
         call error_out("Input Grid File has more than one base",__FILE__,__LINE__)
      end if

      call cg_base_read_f(cgns_file,cgns_base,cgns_git_basename,Dimen,PhysDim,ierror)

      allocate(isize(Dimen,3))
      allocate(istart(Dimen))

      nFace = Dimen * 2
      nCorner = 2**Dimen
      istart = 1

      call cg_nzones_f(cgns_file,cgns_base,nBlock,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      allocate(block(nBlock))
      allocate(cgns_git_zonename(nBlock))

      do b = 1,nBlock
         cgns_zone = b
         call cg_zone_read_f(cgns_file,cgns_base,cgns_zone,cgns_git_zonename(b),isize,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()

         call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
         if (zonetype /= Structured) then
            call error_out("Only Structured Grid supported."//TRIM(file_git_in),__FILE__,__LINE__)
         end if
         block(b) % nPkt = 1
         block(b) % nCell = 1
         block(b) % nPkt(1:Dimen) = int(isize(1:Dimen,1))
         block(b) % nCell(1:Dimen) = int(isize(1:Dimen,2))
         nCell = nCell + product(block(b) % nCell)
         allocate (block(b) % xyz(1-n_BC_Cell : block(b)%nPkt(1)+n_BC_Cell &
                                 ,1-n_BC_Cell : block(b)%nPkt(2)+n_BC_Cell &
                                 ,1-n_BC_Cell : block(b)%nPkt(3)+n_BC_Cell &
                                 ,Dimen))
         allocate (temp_coord(block(b)%nPkt(1),block(b)%nPkt(2),block(b)%nPkt(3)))
         do d = 1,Dimen
            call cg_coord_read_f(cgns_file,cgns_base,cgns_zone,coord_name(d),RealDouble &
                                ,istart,isize(:,1),temp_coord,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()
            block(b) % xyz(1:block(b)%nPkt(1) &
                          ,1:block(b)%nPkt(2) &
                          ,1:block(b)%nPkt(3),d) = temp_coord
         end do
         deallocate (temp_coord)

      end do

      call cg_close_f(cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
      write(*,*) "GIT IN"
      write(*,*) dimen,nblock,ncell
      do b = 1,nBlock
         write(*,*) b,block(b) % nCell
      end do

   end subroutine read_git

   subroutine write_sol()
      use data, only : nBlock, block
      use control, only : iter,sol_out,file_sol_out,iter_sol_out,max_iter,file_git_in
      implicit none
      integer(kind=CGSIZE_T) :: cgns_file,ierror,cgns_base,cgns_zone,cgns_sol,cgns_var,cgns_coord
      integer(kind=CGSIZE_T),allocatable :: isize(:,:)
      integer(kind=CGSIZE_T) :: iter_in_file
      character(len=32)  :: solname
      character(len=100) :: linkname
      character(len=11),parameter :: coordnames(3) = (/"CoordinateX","CoordinateY","CoordinateZ"/)
      integer(kind=CGSIZE_T) :: b,d
      if (mod(iter,iter_sol_out) == 0 .or. iter == max_iter) then
         sol_out = .true.
      end if
      if (sol_out) then
         sol_out = .false.
         write(*,*) "Writing Solution to File"

         write(solname,'(I0)') iter
         if (write_sol_header) then
!            write(*,*) "Writing New Solution to File "//trim(file_sol_out)
            write_sol_header = .false.

            call cg_open_f(trim(file_sol_out),CG_MODE_WRITE,cgns_file,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()

            call cg_base_write_f(cgns_file,"SOLUTION",Dimen,Dimen,cgns_base,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()

            allocate(isize(Dimen,3))
            isize = 0

            do b = 1, nBlock

               isize(:,1) = block(b) % nPkt(1:Dimen)
               isize(:,2) = block(b) % nCell(1:Dimen)

               call cg_zone_write_f(cgns_file,cgns_base,cgns_git_zonename(b),isize,Structured,cgns_zone,ierror)

               if (ierror /= CG_OK) call cg_error_exit_f()
               do d = 1, Dimen
                  call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,coordnames(d) &
                                       ,block(b) % xyz(1:block(b) % nPkt(1)            &
                                                      ,1:block(b) % nPkt(2)            &
                                                      ,1:block(b) % nPkt(3),d),cgns_coord,ierror)
                  if (ierror /= CG_OK) call cg_error_exit_f()
               end do

!               write(*,*) "WRITING GRID LINK"
!               write(linkname,'(A,"/",A,"/GridCoordinates")') trim(cgns_git_basename),trim(cgns_git_zonename(b))
!               write(*,*) trim(linkname)
!               call cg_goto_f(cgns_file,cgns_base,ierror,cgns_git_zonename(b),0,"end")
!               if (ierror /= CG_OK) call cg_error_exit_f()
!               call cg_link_write_f("GridCoordinates",trim(file_git_in),trim(linkname),ierror)
!               if (ierror /= CG_OK) call cg_error_exit_f()

      !         call cg_ziter_write_f(cgns_file,cgns_base,cgns_zone,'ZoneIterativeData',ierror)
      !         call cg_goto_f(cgns_file,cgns_base,ierror,'Zone_t',cgns_zone,'ZoneIterativeData_t',1,'end')
      !         call cg_array_write_f('FlowSolutionPointers',Character,1,32,solname,ierror)

            end do
            call cg_biter_write_f(cgns_file,cgns_base,"TimeIterValues",1,ierror)


      !      call cg_goto_f(cgns_file,cgns_base,ierror,'BaseIterativeData_t',1,'end')
      !      allocate(iter_time(2))
      !      iter_time = dble(iteration)
      !      call cg_array_write_f('TimeValues',RealDouble,1,2,iter_time,ierror)

         else
!            write(*,*) "Appending Solution to existing File: "//trim(file_sol_out)
            call cg_open_f(trim(file_sol_out),CG_MODE_MODIFY,cgns_file,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()
!            write(*,*) "reading nbases"
            call cg_nbases_f(cgns_file,cgns_base,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()

      !      call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
      !      if (ierror /= CG_OK) call cg_error_exit_f()
      !      if (zonetype /= Structured) then
      !         call error_out("Only Structured Grid supported."//TRIM(file_git_in),__FILE__,__LINE__)
      !      end if

            call cg_biter_read_f(cgns_file,cgns_base,linkname,iter_in_file,ierror)
            iter_in_file = iter_in_file + 1
!            write(*,*) "Iterations on File:",iter_in_file
            call cg_biter_write_f(cgns_file,cgns_base,linkname,iter_in_file,ierror)

      !      if (iter_in_file == 2) then
      !
      !      else
      !         call cg_goto_f(cgns_file,cgns_base,ierror,'BaseIterativeData_t',1,'end')
      !         allocate(iter_time(iter_in_file))
      !
      !         iter_time = dble(iteration)
      !         call cg_array_write_f('TimeValues',RealDouble,1,2,iter_time,ierror)
      !
      !      end if
      !      if (iter_in_file == 2) then
      !         do b = 1,nBlock
      !            cgns_zone = b
      !            call cg_ziter_write_f(cgns_file,cgns_base,cgns_zone,'ZoneIterativeData',ierror)
      !            call cg_goto_f(cgns_file,cgns_base,ierror,'Zone_t',cgns_zone,'ZoneIterativeData_t',1,'end')
      !            call cg_array_write_f('FlowSolutionPointers',Character,1,32,solname,ierror)
      !         end do
      !      else
      !      end if
         end if

         do b = 1,nBlock
            cgns_zone = b

            call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,solname,CellCenter,cgns_sol,ierror)
            if (ierror /= CG_OK) call cg_error_exit_f()
            call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
                                 ,"Temperature",block(b) % T(1:block(b) % nCell(1)            &
                                                            ,1:block(b) % nCell(2)            &
                                                            ,1:block(b) % nCell(3),1),cgns_var,ierror)
!            call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
!                                 ,"Heat Transfer Koeff",block(b) % a(1:block(b) % nCell(1)            &
!                                                            ,1:block(b) % nCell(2)            &
!                                                            ,1:block(b) % nCell(3),1),cgns_var,ierror)
!            call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
!                                 ,"Temperature Residual",block(b) % res(1:block(b) % nCell(1)            &
!                                                            ,1:block(b) % nCell(2)            &
!                                                            ,1:block(b) % nCell(3),1),cgns_var,ierror)
         end do
         call cg_close_f(cgns_file,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
      end if
   end subroutine
   subroutine read_bc()
      use control, only : file_bc
      use data, only : nBlock, block
      implicit none
      integer, parameter :: fu = 99
      integer, parameter :: bc_file_version_number = 1001
      real(kind=dp), parameter :: vorfaktor = -2.0E0_dp
      !< vorfaktor für den per Randbedingung gegebenen Heatflux, da ein Faktor 0.5 bei der RES-Berechnung verwendet wird
      !< und die Gradienten in negativer Richtung berechnet werden.

      integer :: b, f, bc
      integer :: r_bc_file_version_number,r_nBlock,nBC,BC_TYPE

      logical :: fexists

      inquire(file=trim(file_bc),exist=fexists)

      if(.not. fexists) then
        call error_out("GITTER INPUT DATEI konnte nicht gefunden werden: "//TRIM(file_bc),__FILE__,__LINE__)
      end if

      open(unit=fu,file=trim(file_bc),form="UNFORMATTED",access="STREAM",status="OLD")
         read(fu) r_bc_file_version_number,r_nBlock
         if (r_bc_file_version_number /= bc_file_version_number) then
            call error_out("BC_File Version is not correct!",__FILE__,__LINE__)
         end if
         if (r_nBlock /= nBlock) then
            call error_out("Number of Blocks in BC File is wrong!",__FILE__,__LINE__)
         end if
         do b = 1,nBlock
            do f = 1,6
               read(fu) nBC
               block(b) % face(f) % nBC = nBC
               allocate(block(b) % face(f) % BC(nBC))
               do bc = 1, nBC
                  read(fu) BC_Type
                  block(b) % face(f) % BC(bc) % BC_Type = BC_Type
                  read(fu) block(b) % face(f) % BC(bc) % dist
                  block(b) % face(f) % BC(bc) % ncell = block(b) % face(f) % BC(bc) % dist(2,:) &
                                                      - block(b) % face(f) % BC(bc) % dist(1,:) + 1

                  if (BC_Type == BC_ISOTHERMAL .or. BC_Type == BC_HEATFLUX) then
                     allocate(block(b) % face(f) % BC(bc) % data(block(b) % face(f) % BC(bc) % ncell(2) &
                                                                ,block(b) % face(f) % BC(bc) % ncell(1)))
                     read(fu) block(b) % face(f) % BC(bc) % data
                     if (BC_Type == BC_HEATFLUX) then
                        block(b) % face(f) % BC(bc) % data = block(b) % face(f) % BC(bc) % data * vorfaktor
                     end if
                  else if (BC_Type > 0) then
                     read(fu) block(b) % face(f) % BC(bc) % CPU_id
                     read(fu) block(b) % face(f) % BC(bc) % a

                  end if

                  write(*,*) b,f,bc,BC_Type,block(b) % face(f) % BC(bc) % dist, block(b) % face(f) % BC(bc) % ncell

               end do
            end do
         end do
      close(fu)
   end subroutine read_bc

   subroutine read_sol()
      use data, only : nBlock, block
      use control, only : file_sol_in
      implicit none
      logical :: fexists, t_found
      integer(kind=CGSIZE_T) :: cgns_file,ierror,cgns_base,PhysDim,cgns_zone
      integer(kind=CGSIZE_T),allocatable :: isize(:,:),istart(:)
      integer(kind=CGSIZE_T) :: data_location,zonetype,nSol,nVar_in,datatype,var
      integer(kind=CGSIZE_T) :: rDimen,rnBlock
      integer(kind=CGSIZE_T) :: b,n
      character(len=32)  :: solname,varname_in

      inquire(file=trim(file_sol_in),exist=fexists)

      if(.not. fexists) then
        call error_out("RESTART DATEI konnte nicht gefunden werden: "//TRIM(file_sol_in),__FILE__,__LINE__)
      end if


      call cg_open_f(trim(file_sol_in),CG_MODE_READ,cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      call cg_nbases_f(cgns_file,cgns_base,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      if (cgns_base /= 1) then
         call error_out("Restart File has more than one base",__FILE__,__LINE__)
      end if
      call cg_base_read_f(cgns_file,cgns_base,cgns_basename,rDimen,PhysDim,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

      if (rDimen /= Dimen) then
         call error_out("Restart File has differet dimensions than the grid.",__FILE__,__LINE__)
      end if

      allocate(isize(Dimen,3))
      allocate(istart(Dimen))

      istart = 1

      call cg_nzones_f(cgns_file,cgns_base,rnBlock,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()
      if (nBlock /= rnBlock) then
         call error_out("Sol #Blocks is different than Grid #Blocks",__FILE__,__LINE__)
      end if
      do b = 1,nBlock
         cgns_zone = b
         call cg_zone_read_f(cgns_file,cgns_base,cgns_zone,cgns_git_zonename(b),isize,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()

         call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
         if (zonetype /= Structured) then
            call error_out("Only Structured Grid supported."//TRIM(file_sol_in),__FILE__,__LINE__)
         end if
         do n = 1,Dimen
            if (block(b) % nPkt(n) /= isize(n,1)) then
               call error_out("nPkt do not match",__FILE__,__LINE__)
            end if
            if (block(b) % nCell(n) /= isize(n,2)) then
               call error_out("nCell do not match.",__FILE__,__LINE__)
            end if
         end do
         !!!! CHECKING if more than one solution.
         call cg_nsols_f(cgns_file,cgns_base,cgns_zone,nSol,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
         if (nSol /= 1) then
            call error_out("More than One Solution in Restart-File",__FILE__,__LINE__)
         end if
         !!!! Checking if Cell-Centered
         call cg_sol_info_f(cgns_file,cgns_base,cgns_zone,1,solname,data_location,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
         if (data_location .ne. CellCenter) then
            call error_out("Not Cell-Centered Data",__FILE__,__LINE__)
         end if
         call cg_nfields_f(cgns_file,cgns_base,cgns_zone,nSol,nVar_in,ierror)
         if (ierror /= CG_OK) call cg_error_exit_f()
         t_found = .false.
         do var = 1, nVar_in
            call cg_field_info_f(cgns_file,cgns_base,cgns_zone,nSol,var,datatype,varname_in,ierror)
            if (varname_in == "Temperature") then
            call cg_field_read_f(cgns_file,cgns_base,cgns_zone,nSol,varname_in,datatype       &
                                   ,istart,isize(:,2),block(b) % t(1:block(b) % nCell(1)            &
                                                                  ,1:block(b) % nCell(2)            &
                                                                  ,1:block(b) % nCell(3),1),ierror)
            t_found = .true.
            end if
         end do
         if (.not.t_found) then
            call error_out("Temperature not found in Restart File",__FILE__,__LINE__)

         end if
      end do

      call cg_close_f(cgns_file,ierror)
      if (ierror /= CG_OK) call cg_error_exit_f()

   end subroutine
end module io
