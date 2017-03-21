module mpnt
   use const
   use data
   use control, only : iter
   implicit none
   type tmpnt
      integer                    :: block
      integer                    :: ijk(8,3)
      real(kind=dp)              :: w(8)
   end type tmpnt

   type(tmpnt) , allocatable     :: mpnts(:)
   integer :: nmpnt
   character(len=100) :: file_mpnt_out = "mpnt.dat"
   character(len=100) :: file_mpnt_in  = "mpnt.cfg"
   integer, parameter :: fu = 666
   contains
      subroutine write_mpnt()
         implicit none
         integer :: p
         open(fu,file=trim(file_mpnt_out), status = "old", position = "append" )

         write(fu,'(I0)',ADVANCE="NO") iter
         do p = 1, nmpnt
            write(fu,'(1X,ES12.5)',ADVANCE="NO") block(mpnts(p) % block) % t( &
                                                       mpnts(p) % ijk(1,1)    &
                                                      ,mpnts(p) % ijk(1,2)    &
                                                      ,mpnts(p) % ijk(1,3),1  )
         end do

         write(fu,*)

         close(fu)
      end subroutine write_mpnt

      subroutine init_mpnt()
         implicit none
         integer :: p

         nmpnt = 4
         allocate( mpnts(nmpnt) )

         p = 1
         mpnts(p) % ijk(1,:) = (/1                  ,block(1) % nCell(2),1                  /)
         mpnts(p) % block = 1
         p = 2
         mpnts(p) % ijk(1,:) = (/block(1) % nCell(1),block(1) % nCell(2),1                  /)
         mpnts(p) % block = 1
         p = 3
         mpnts(p) % ijk(1,:) = (/1                  ,block(1) % nCell(2),block(1) % nCell(3)/)
         mpnts(p) % block = 1
         p = 4
         mpnts(p) % ijk(1,:) = (/block(1) % nCell(1),block(1) % nCell(2),block(1) % nCell(3)/)
         mpnts(p) % block = 1
         open(fu,file=trim(file_mpnt_out))

!         do p = 1, nmpnt
!            write(fu,'(1X,4(I3.3))',ADVANCE="NO")       mpnts(p) % block      &
!                                                      ,mpnts(p) % ijk(1,1)    &
!                                                      ,mpnts(p) % ijk(1,2)    &
!                                                      ,mpnts(p) % ijk(1,3)
!         end do

         write(fu,*)

         close(fu)
      end subroutine
end module mpnt
