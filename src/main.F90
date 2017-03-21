program HEAT
   use io, only : read_git,write_sol,read_bc,read_sol
   use control
   !use sparse_mat
   use solver
   use data!, only: init,block,nBlock
   use flux
   use time_int, only: update_sol, calc_timestep, set_timestep
   use boundary, only: update_boundary, boundary_fluxes
   use mpnt, only : init_mpnt, write_mpnt
   use lhs
   implicit none
   integer :: b,n

   write(*,'(A)') "HEAT SOLVER"
   call read_git()
   call read_bc()
   call init()
   call read_sol()
   call init_mpnt()
   do iter = 1, max_iter
      do inner_iter = 1, n_inner_iter
         if (.not. given_time_step) time_step = 1E10_dp
         do b = 1, nBlock
            if (.not.implicit) then
               !!! UPDATE BOUNDARY CELLS ACCORDING TO BOUNDARY CODITIONS
               call update_boundary(block       = block                             &
                                   ,nblock      = nBlock                            &
                                   ,inner_iter  = inner_iter)

               !!! CALCULATE FLUXES
               call calc_flux ( T               = block(b) % T   (:,:,:,inner_iter)   &
                              , a               = block(b) % a   (:,:,:,inner_iter)   &
                              , flux            = block(b) % flux                     &
                              , dn              = block(b) % dn                       )

               !!! CALCULATE BOUDARY FLUXES (HEATFLUX BC)
               call boundary_fluxes(flux        = block(b) % flux                     &
                                   ,face        = block(b) % face                     )

               !!! CALCULATE CELL RESIDIUMS
               call calc_res  ( res             = block(b) % res (:,:,:,inner_iter)   &
                              , flux            = block(b) % flux                     &
                              , edge_area       = block(b) % edge_area              )
            end if
            call calc_timestep ( dt          = block(b) % dt                     &
                               , a           = block(b) % a (:,:,:,inner_iter)   &
                               , dn2_dt      = block(b) % dn2_dt                 )
         end do
         do b = 1, nBlock
            call set_timestep (  dt          = block(b) % dt                     )
            if (.not.implicit) then
               !!!! UPDATE SOLUTION
               call update_sol( T               = block(b) % T  (:,:,:,inner_iter)  &
                              , res             = block(b) % res(:,:,:,inner_iter)  &
                              , dt              = block(b) % dt                     &
                              , cell_vol        = block(b) % cell_vol               )
            end if
         end do
         if (implicit) then
!            row_pos = row_pos_save
!            col_pos = col_pos_save
            call update_implicit ( mat = mat &
                                 , rhs = rhs &
                                 , sol = sol &
                                 , block = block)
!            mat2 = mat
!            write(*,*) rhs


            call solve_system(nCell,rhs,sol,mat,jmat,imat,1)
!           call solve_sys (nCell,n_nonzero_elements,vec_size1,vec_size2,vec_size3 &
!                           ,mat,row_pos,col_pos,rhs &
!                           ,pivot,ha,real_flag,int_flag,error_flag)

!            do n = 1,vec_size1
!               write(*,*) mat(n),mat2(n)
!            end do
!            stop
!            write(*,*) rhs
!            write(*,*) block(1) % T(1,1,1,1)
            call update_sol_imp( rhs = sol &
                           , block = block)
!            write(*,*) block(1) % T(1,1,1,1)
         end if
      end do
      call res_control()
      call write_sol()
      call write_mpnt()
   end do
   write(*,'(A)') "DONE!"
end program HEAT
