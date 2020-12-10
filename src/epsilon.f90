!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE epsilon
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the dielectric function
  !
  USE kinds,                  ONLY : dp
  USE constants,              ONLY : pi
  USE cell_base,              ONLY : omega
  USE io_global,              ONLY : stdout, ionode, ionode_id
  USE io_files,               ONLY : nwordwfc, iunwfc
  USE wavefunctions,          ONLY : evc
  USE noncollin_module,       ONLY : npol, noncolin
  USE klist,                  ONLY : nks, nkstot, ngk, wk
  USE wvfct,                  ONLY : nbnd, npwx, et
  USE lsda_mod,               ONLY : nspin
  USE orbm_module,            ONLY : ry2ha, ci, vel_evc, eta, wmin, wmax, &
                                     nw, nbnd_occ, sigma, ry2ev
  USE buffers,                ONLY : get_buffer
  USE mp_pools,               ONLY : my_pool_id, me_pool, root_pool,  &
                                     inter_pool_comm, intra_pool_comm, npool
  USE mp,                     ONLY : mp_sum, mp_bcast
  USE mp_world,               ONLY : world_comm
  USE symme,                  ONLY : symmatrix

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! the following three quantities are for norm-conserving PPs
  complex(dp), allocatable, dimension(:,:,:) :: vmat            ! <n|v|m> 
  
  real(dp) :: eps2(nw, 3, 3)
  
  integer :: ik, ios, iunout
  integer :: i, j, ibnd, jbnd, n, ie
  real(dp) :: det, wgrid(nw), pref
  real(dp), external :: get_clock, delta
  integer, external :: find_free_unit
  integer :: npw

 
  call start_clock('epsilon')
  !-----------------------------------------------------------------------
  ! allocate memory
  !-----------------------------------------------------------------------
  allocate ( vmat(nbnd, nbnd, 3), vel_evc(npwx*npol, nbnd, 3))

  ! print memory estimate
  call orbm_memory_report

  write(stdout, '(5X,''Computing the dielectric function:'',$)')
  write(stdout, *)
 
  do i = 1, nw
     wgrid(i) = wmin + (wmax-wmin)/(nw-1)*(i-1)
  enddo
  eps2 = 0.d0 
  !====================================================================
  ! loop over k-points on the pool
  !====================================================================
  do ik = 1, nks

    npw = ngk(ik)

    ! read wfcs from file and compute becp
    call get_buffer (evc, nwordwfc, iunwfc, ik)
    
    
    ! calculate du/dk    
    vel_evc(:,:,:) = (0.d0,0.d0)
    vmat(:,:,:)= (0.d0,0.d0)
    do i = 1,3
       call apply_vel(evc, vel_evc(1,1,i), ik, i)
       
       ! calculate the velocity matrix ps(nbnd,nbnd)
       if (noncolin) then
     
          CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), evc(1,1), &
                    npwx*npol, vel_evc(1,1,i), npwx*npol, (0.d0,0.d0), vmat(1,1,i), nbnd)
       else
       
          CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), evc(1,1), &
                    npwx, vel_evc(1,1,i), npwx, (0.d0,0.d0), vmat(1,1,i), nbnd)
       endif
       
    enddo
   
    
#ifdef __MPI
    call mp_sum(vmat, intra_pool_comm)
#endif 

    do ie = 1, nw
       pref = 4*pi**2/omega/wgrid(ie)**2*wk(ik)
       
       do i = 1, 3
          do j = 1, 3
          
             do ibnd = 1, nbnd_occ(ik)
                do jbnd = 1, nbnd_occ(ik)
                
                   det = (et(jbnd,ik) - et(ibnd,ik))*ry2ha
                   if (det > wmin - 8*sigma .and. det < wmax + 8*sigma &
                       .and. det-wgrid(ie) < 8*sigma) then
                       
                       eps2(ie, j, i) = eps2(ie, j, i) + pref* &
                               CONJG(vmat(jbnd, ibnd, j))*vmat(jbnd, ibnd, i)* &
                               delta(det-wgrid(ie))
                   endif
                   
                enddo
             enddo
             
          enddo
       enddo
       
    enddo
                    

  enddo ! ik
  
#ifdef __MPI
    call mp_sum(eps2, inter_pool_comm)
#endif 
  
  do ie = 1, nw
     call symmatrix( eps(ie, :, :) )
  enddo
  
  ios = 0
  if ( ionode ) then
     iunout = find_free_unit( )
     open (unit = iunout, file = 'epsilon_imag.dat', status = 'unknown', iostat = ios)
     rewind (iunout)
  endif

  call mp_bcast (ios, ionode_id, world_comm)
  if ( ios/=0 ) call errore ('epsilon', 'Opening file epsilon_imag.dat', abs (ios) )

  if (ionode) then
     write(iunout, '(7A10)') '#   ENERGY','XX','YY','ZZ','XY','YZ','ZX' 
     do ie = 1, net
        write(iunout, '(7F10.4)') wgrid(ie)*ry2ev/ry2ha, eps2(ie,1,1),eps2(ie,2,2), &
                       eps2(ie,3,3), eps2(ie,1,2),eps2(ie,2,3),eps2(ie,3,1)
     enddo
     close(iunout)
  endif
  
  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,*)
  write(stdout,'(5X,''End of dielectric function calculation'')')
  write(stdout,*)


  ! free memory as soon as possible
  deallocate( vel_evc, vmat)

  
  !call restart_cleanup ( )
  call stop_clock('epsilon')

END SUBROUTINE epsilon
