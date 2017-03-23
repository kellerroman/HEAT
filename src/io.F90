module io
!   use cgnslib
!   use cgns_types, ONLY: CGSIZE_T
   use hdf5
   use const
   use control, only: Dimen,nFace,nCorner,n_BC_Cell
   implicit none
   logical :: write_sol_header = .true.

   character(len=*), parameter :: GROUP_GRID            = "grid"
   character(len=*), parameter :: GROUP_DATA            = "data"
   character(len=*), parameter :: GROUP_BLOCK           = "block"
   character(len=*) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]
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

   subroutine read_sol()
      use data, only : nCell, nBlock, block
      use control, only: file_git_in,n_inner_iter
      implicit none

      integer(hid_t) :: file_id       ! file identifier
      integer(hid_t) :: dset_id       ! dataset identifier
      integer(hid_t) :: group_id      ! dataset identifier
      integer(hid_t) :: group_id1     ! dataset identifier
      integer(hid_t) :: group_id2     ! dataset identifier
      integer(hid_t) :: dspace_id     ! dataspace identifier
      integer        :: solution_type ! dataspace identifier
      integer        :: var_type ! dataspace identifier
      integer        ::   hdf5_nSol, nVar_in
      character(len=len(GROUP_BLOCK)+2) :: block_group
      integer(HSIZE_T) :: dims(3)
      integer(HSIZE_T) :: maxdims(3)
      integer(kind=8) :: b ,d
      logical :: fexists
      integer     ::   error ! Error flag
      integer(kind=8),allocatable :: isize(:,:),istart(:)
      real(kind=dp), allocatable :: data_in(:,:,:)
      character(len=20) :: varname_in
      character(len=10) :: solution_name

      nCell = 0

      inquire(file=trim(file_git_in),exist=fexists)

      if(.not. fexists) then
        call error_out("GITTER INPUT DATEI konnte nicht gefunden werden: " &
                       //TRIM(file_git_in),__FILE__,__LINE__)
      end if


      ! Initialize FORTRAN interface.
      call h5open_f(error)

      ! Open an existing file.
      call h5fopen_f (trim(file_git_in), h5f_acc_rdwr_f, file_id, error)

      call h5gopen_f(file_id,GROUP_GRID,group_id,error)

      call h5gn_members_f(file_id, GROUP_GRID, nBlock, error)
      !!!!!!!!! ONLY ONE BLOCK SUPPORTED AT THE MOMENT
      Dimen = 3

      allocate(isize(Dimen,3))
      allocate(istart(Dimen))

      nFace = Dimen * 2
      nCorner = 2**Dimen
      istart = 1


      allocate(block(nBlock))

      do b = 1,nBlock
         write(block_group,'(A,I0)') GROUP_BLOCK, b
         call h5gopen_f(group_id,block_group,group_id2,error)
         block(b) % nPkt = 1
         block(b) % nCell = 1
         call h5dopen_f(group_id2, COORD_NAME(1), dset_id, error)
         call h5dget_space_f(dset_id,dspace_id,error)
         
         call h5sget_simple_extent_ndims_f(dspace_id,dimen,error)
         call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)

         call h5dclose_f(dset_id, error)
         block(b) % nPkt = INT(dims,ip)
         block(b) % nCell(1:Dimen) = max(0,block(b) % nPkt - 1)
         nCell = nCell + product(block(b) % nCell)
         allocate (block(b) % xyz(1-n_BC_Cell : block(b)%nPkt(1)+n_BC_Cell &
                                 ,1-n_BC_Cell : block(b)%nPkt(2)+n_BC_Cell &
                                 ,1-n_BC_Cell : block(b)%nPkt(3)+n_BC_Cell &
                                 ,Dimen))
         allocate (block(b) % T  ( 1-n_BC_Cell:block(b) % nCell(1)+n_BC_Cell &
                                 , 1-n_BC_Cell:block(b) % nCell(2)+n_BC_Cell &
                                 , 1-n_BC_Cell:block(b) % nCell(3)+n_BC_Cell &
                                 , n_inner_iter ) )
         allocate (data_in(block(b)%nPkt(1),block(b)%nPkt(2),block(b)%nPkt(3)))

         do d = 1,Dimen
            call h5dopen_f(group_id2, COORD_NAME(d), dset_id, error)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_in, dims, error)
            block(b) % xyz(1:block(b)%nPkt(1) &
                          ,1:block(b)%nPkt(2) &
                          ,1:block(b)%nPkt(3),d) = data_in
            call h5dclose_f(dset_id, error)
         end do
         deallocate (data_in)
         call h5gclose_f(group_id2, error)

      end do
      call h5gclose_f(group_id, error) ! CLOSE GRID GROUP
      write(*,*) 
      do b = 1, nBlock
         write(*,*) b, block(b) % nCell
      end do
      !!!!
      !!!! READING GRID DONE
      !!!!

      call h5gopen_f(file_id,GROUP_DATA,group_id,error)
      call h5gn_members_f(file_id, GROUP_DATA, hdf5_nSol, error)
      !!!!!!!!! ONLY ONE BLOCK SUPPORTED AT THE MOMENT
      if (hdf5_nSol > 1) then
         call error_out("More than One Solution in Restart-File",__FILE__,__LINE__)
      end if
      call h5gget_obj_info_idx_f(file_id, GROUP_DATA, 0,solution_name, solution_type, error)
      !write(*,*) solution_name
      call h5gopen_f(group_id,solution_name,group_id1,error) ! OPEN TIMESTEP GROUP
      do b = 1, nBlock
         write(block_group,'(A,I0)') GROUP_BLOCK,b
         !write(*,*) "opening", block_group
         call h5gopen_f(group_id1,block_group,group_id2,error)
         call h5gn_members_f(group_id1, block_group, nVar_in, error)
         allocate (data_in(block(b)%nCell(1),block(b)%nCell(2),block(b)%nCell(3)))
         call h5gget_obj_info_idx_f(group_id1, block_group, 0,varName_in, var_type, error)
         call h5dopen_f(group_id2, varName_in, dset_id, error)
         call h5dget_space_f(dset_id,dspace_id,error)
         
         call h5sget_simple_extent_dims_f(dspace_id,dims,maxdims,error)

         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_in, dims, error)
         call h5dclose_f(dset_id, error)
         block(b) % T (1:block(b) % nCell(1)            &
                      ,1:block(b) % nCell(2)            &
                      ,1:block(b) % nCell(3),1) = data_in
         deallocate(data_in)
         call h5gclose_f(group_id2, error)
      end do
      
      call h5gclose_f(group_id1, error) !CLOSE TIMESTEP GROUP
      call h5gclose_f(group_id, error)  !CLOSE DATA GROUP
      ! close the file.
      call h5fclose_f(file_id, error)
      
      ! close fortran interface.
      call h5close_f(error)
      write(*,*) "Datin DONE"
   end subroutine read_sol

   subroutine write_sol()
      use data, only : nBlock, block
      use control, only : iter,sol_out,file_sol_out,iter_sol_out,max_iter,file_git_in
      implicit none
      integer(kind=8) :: cgns_file,ierror,cgns_base,cgns_zone,cgns_sol,cgns_var,cgns_coord
      integer(kind=8),allocatable :: isize(:,:)
      integer(kind=8) :: iter_in_file
      character(len=32)  :: solname
      character(len=100) :: linkname
      character(len=11),parameter :: coordnames(3) = (/"CoordinateX","CoordinateY","CoordinateZ"/)
      integer(kind=8) :: b,d
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

!            call cg_open_f(trim(file_sol_out),CG_MODE_WRITE,cgns_file,ierror)
!            if (ierror /= CG_OK) call cg_error_exit_f()

!            call cg_base_write_f(cgns_file,"SOLUTION",Dimen,Dimen,cgns_base,ierror)
!            if (ierror /= CG_OK) call cg_error_exit_f()

            allocate(isize(Dimen,3))
            isize = 0

            do b = 1, nBlock

               isize(:,1) = block(b) % nPkt(1:Dimen)
               isize(:,2) = block(b) % nCell(1:Dimen)

!               call cg_zone_write_f(cgns_file,cgns_base,cgns_git_zonename(b),isize,Structured,cgns_zone,ierror)

!               if (ierror /= CG_OK) call cg_error_exit_f()
               do d = 1, Dimen
!                  call cg_coord_write_f(cgns_file,cgns_base,cgns_zone,RealDouble,coordnames(d) &
!                                       ,block(b) % xyz(1:block(b) % nPkt(1)            &
!                                                      ,1:block(b) % nPkt(2)            &
!                                                      ,1:block(b) % nPkt(3),d),cgns_coord,ierror)
!                  if (ierror /= CG_OK) call cg_error_exit_f()
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
!            call cg_biter_write_f(cgns_file,cgns_base,"TimeIterValues",1,ierror)


      !      call cg_goto_f(cgns_file,cgns_base,ierror,'BaseIterativeData_t',1,'end')
      !      allocate(iter_time(2))
      !      iter_time = dble(iteration)
      !      call cg_array_write_f('TimeValues',RealDouble,1,2,iter_time,ierror)

         else
!            write(*,*) "Appending Solution to existing File: "//trim(file_sol_out)
!            call cg_open_f(trim(file_sol_out),CG_MODE_MODIFY,cgns_file,ierror)
!            if (ierror /= CG_OK) call cg_error_exit_f()
!            write(*,*) "reading nbases"
!            call cg_nbases_f(cgns_file,cgns_base,ierror)
!            if (ierror /= CG_OK) call cg_error_exit_f()

      !      call cg_zone_type_f(cgns_file,cgns_base,cgns_zone,zonetype,ierror)
      !      if (ierror /= CG_OK) call cg_error_exit_f()
      !      if (zonetype /= Structured) then
      !         call error_out("Only Structured Grid supported."//TRIM(file_git_in),__FILE__,__LINE__)
      !      end if

!            call cg_biter_read_f(cgns_file,cgns_base,linkname,iter_in_file,ierror)
            iter_in_file = iter_in_file + 1
!            write(*,*) "Iterations on File:",iter_in_file
!            call cg_biter_write_f(cgns_file,cgns_base,linkname,iter_in_file,ierror)

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

!            call cg_sol_write_f(cgns_file,cgns_base,cgns_zone,solname,CellCenter,cgns_sol,ierror)
!            if (ierror /= CG_OK) call cg_error_exit_f()
!            call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
!                                 ,"Temperature",block(b) % T(1:block(b) % nCell(1)            &
!                                                            ,1:block(b) % nCell(2)            &
!                                                            ,1:block(b) % nCell(3),1),cgns_var,ierror)
!            call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
!                                 ,"Heat Transfer Koeff",block(b) % a(1:block(b) % nCell(1)            &
!                                                            ,1:block(b) % nCell(2)            &
!                                                            ,1:block(b) % nCell(3),1),cgns_var,ierror)
!            call cg_field_write_f(cgns_file,cgns_base,cgns_zone,cgns_sol,RealDouble         &
!                                 ,"Temperature Residual",block(b) % res(1:block(b) % nCell(1)            &
!                                                            ,1:block(b) % nCell(2)            &
!                                                            ,1:block(b) % nCell(3),1),cgns_var,ierror)
         end do
!         call cg_close_f(cgns_file,ierror)
!         if (ierror /= CG_OK) call cg_error_exit_f()
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
        call error_out("GITTER INPUT DATEI konnte nicht gefunden werden: " &
                       //TRIM(file_bc),__FILE__,__LINE__)
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

end module io
