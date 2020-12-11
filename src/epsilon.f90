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
  ! eps2 = 4*pi^2*e^2/V/omega^2 sum_k |vk|^2 delta(omega_cv - omega)
  ! eps1 can be obtained by eps2 through Kramers-Kronig transformation
  ! eps1 = 1 + 8*pi*e^2/V/omega_cv sum_k |vk|^2 /(omega_cv^2 - omega^2) 
  ! delta function replaced by lorentzian
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
  USE symm_base,              ONLY : nsym
  USE symme,                  ONLY : symmatrix

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! the following three quantities are for norm-conserving PPs
  complex(dp), allocatable, dimension(:,:,:) :: vmat            ! <n|v|m> 
  
  real(dp) :: eps2(nw, 3, 3)
  real(dp) :: eps1(nw, 3, 3)
  
  integer :: ik, ios, iunout
  integer :: i, j, ibnd, jbnd, n, ie
  real(dp) :: det, wgrid(nw), pref
  real(dp), external :: get_clock
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
  eps1 = 0.d0 
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
       pref = 4*pi**2/omega*wk(ik)
       
       do i = 1, 3
          do j = 1, 3
          
             do ibnd = 1, nbnd_occ(ik)
                do jbnd = nbnd_occ(ik)+1, nbnd
                
                   det = (et(jbnd,ik) - et(ibnd,ik))*ry2ha
                   !if (det > wmin - 8*sigma .and. det < wmax + 8*sigma &
                   !    .and. ABS(det-wgrid(ie)) < 8*sigma) then
                       
                       eps2(ie,j,i) = eps2(ie,j,i) + pref*CONJG(vmat(jbnd,ibnd,j))* &
                                      vmat(jbnd,ibnd,i)*sigma/pi/((det-wgrid(ie))**2 &
                                      + sigma**2)/det**2
                       eps1(ie,j,i) = eps1(ie,j,i) + pref*CONJG(vmat(jbnd,ibnd,j))* &
                                      vmat(jbnd,ibnd,i)*2/pi*(det-wgrid(ie))/((det-wgrid(ie))**2 &
                                      + sigma**2)/det/(det+wgrid(ie))
                   !endif
                   
                enddo
             enddo
             
          enddo
       enddo
       
    enddo
                    

  enddo ! ik
  
#ifdef __MPI
    call mp_sum(eps1, inter_pool_comm)
    call mp_sum(eps2, inter_pool_comm)
#endif 

  do i = 1,3
     eps1(:,i,i) = 1.d0 + eps1(:,i,i)
  enddo
  if (nsym > 1) then
     do ie = 1, nw
        call symmatrix( eps2(ie, :, :) )
        call symmatrix( eps1(ie, :, :) )
     enddo
  endif
  

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
     do ie = 1, nw
        write(iunout, '(7F10.4)') wgrid(ie)*ry2ev/ry2ha, eps2(ie,1,1),eps2(ie,2,2), &
                       eps2(ie,3,3), eps2(ie,1,2),eps2(ie,2,3),eps2(ie,3,1)
     enddo
     close(iunout)
  endif

   ios = 0
  if ( ionode ) then
     iunout = find_free_unit( )
     open (unit = iunout, file = 'epsilon_real.dat', status = 'unknown', iostat = ios)
     rewind (iunout)
  endif

  call mp_bcast (ios, ionode_id, world_comm)
  if ( ios/=0 ) call errore ('epsilon', 'Opening file epsilon_real.dat', abs (ios) )

  if (ionode) then
     write(iunout, '(7A10)') '#   ENERGY','XX','YY','ZZ','XY','YZ','ZX' 
     do ie = 1, nw
        write(iunout, '(7F10.4)') wgrid(ie)*ry2ev/ry2ha, eps1(ie,1,1),eps1(ie,2,2), &
                       eps1(ie,3,3), eps1(ie,1,2),eps1(ie,2,3),eps1(ie,3,1)
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


