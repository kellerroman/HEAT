module io
   use hdf5
   use const
   use control, only: Dimen,nFace,nCorner,n_BC_Cell
   implicit none
   logical :: write_sol_header = .true.

   integer         , parameter :: VARNAME_LENGTH        = 20

   character(len=*), parameter :: GROUP_GRID            = "grid"
   character(len=*), parameter :: GROUP_DATA            = "data"
   character(len=*), parameter :: GROUP_BLOCK           = "block"
   character(len=*) , parameter   :: COORD_NAME(3)      = [ "CoordinateX","CoordinateY","CoordinateZ" ]
   character(len=VARNAME_LENGTH), parameter :: VARNAME   = "Temperatur"     ! dataset name

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
      character(len=32)  :: solname
      
      integer(hsize_t) :: dims(3)
      real(kind=dp), allocatable :: data_out(:,:,:)

      character(len=7)           :: block_group
      integer(hid_t)             :: file_id                          ! file identifier
      integer(hid_t)             :: dset_id                          ! dataset identifier
      integer(hid_t)             :: group_id                         ! dataset identifier
      integer(hid_t)             :: group_id1                        ! dataset identifier
      integer(hid_t)             :: group_id2                        ! dataset identifier
      integer(hid_t)             :: dspace_id                        ! dataspace identifier
      integer                    :: error                            ! error flag
      integer                    :: b,d
      if (mod(iter,iter_sol_out) == 0 .or. iter == max_iter) then
         sol_out = .true.
      end if
      if (sol_out) then
         CALL h5open_f(error)
         sol_out = .false.
         write(*,*) "Writing Solution to File"

         if (write_sol_header) then
            write_sol_header = .false.

            call h5fcreate_f(file_sol_out, h5f_acc_trunc_f, file_id, error)

            call h5gcreate_f(file_id,  GROUP_GRID, group_id,  error)

            do b = 1, nBlock
               dims(:) = block(b) % nPkt  (1:Dimen)
               call h5screate_simple_f(Dimen, dims, dspace_id, error)
               write(block_group,'(A,I0)') GROUP_BLOCK,b
               call h5gcreate_f(group_id, block_group, group_id2, error)
               allocate (data_out ( block(b) % nPkt(1), block(b) % nPkt(2), block(b) % nPkt(3) ))
               do d = 1, Dimen
                  data_out = block(b) % xyz(1:block(b) % nPkt(1)            &
                                           ,1:block(b) % nPkt(2)            &
                                           ,1:block(b) % nPkt(3),d)
                  call h5dcreate_f(group_id2, COORD_NAME(d), h5t_native_double, dspace_id, &
                                  dset_id, error)
                  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dims, error)
                  call h5dclose_f(dset_id, error)
               end do
            ! Close the group.
            call h5gclose_f(group_id2, error)
            deallocate(data_out)
            end do
            call h5gclose_f(group_id, error)
            call h5gcreate_f(file_id,  GROUP_DATA, group_id,  error)
         else
            !write(*,*) "Open excisting file"
            CALL h5fopen_f (file_sol_out, H5F_ACC_RDWR_F, file_id, error)
            call h5gopen_f(file_id,GROUP_DATA,group_id,error)
         end if

         write(solname,'(I10.10)') iter
         call h5gcreate_f(group_id,  trim(solname), group_id1,  error)
         do b = 1,nBlock

            write(block_group,'(A,I0)') GROUP_BLOCK,b
            dims = block(b) % nCell 
            call h5screate_simple_f(Dimen, dims, dspace_id, error)

            ! Create a group named for block1 in the file.
            call h5gcreate_f(group_id1, block_group, group_id2, error)
            allocate (data_out(block(b) % nCell(1),block(b) % nCell(2),block(b) % nCell(3))) 
            data_out = block(b) % T(1:block(b) % nCell(1)  &
                                   ,1:block(b) % nCell(2)  &
                                   ,1:block(b) % nCell(3),1)
            call h5dcreate_f(group_id2, VARNAME, h5t_native_double, dspace_id, &
                 dset_id, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dims, error)
            call h5dclose_f(dset_id, error)
            ! terminate access to the data space.
            call h5sclose_f(dspace_id, error)
            ! Close the group.
            call h5gclose_f(group_id2, error)
            deallocate(data_out)
         end do
         call h5gclose_f(group_id1, error)
         ! Close the group.
         call h5gclose_f(group_id, error)
         ! close the file.
         call h5fclose_f(file_id, error)
         
         ! close fortran interface.
         call h5close_f(error)
      end if
   end subroutine write_sol
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
