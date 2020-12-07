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
  USE orbm_module,            ONLY : ry2ha, ci, vel_evc, eta, etmin, etmax, &
                                     net, nbnd_occ, sigma, ry2ev
  USE buffers,                ONLY : get_buffer
  USE mp_pools,               ONLY : my_pool_id, me_pool, root_pool,  &
                                     inter_pool_comm, intra_pool_comm, npool
  USE mp,                     ONLY : mp_sum, mp_bcast
  USE mp_world,               ONLY : world_comm

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! the following three quantities are for norm-conserving PPs
  complex(dp), allocatable, dimension(:,:,:) :: vmat            ! <n|v|m> 
  
  real(dp) :: eps_imag(net, 6)
  
  integer :: ik, ios, iunout
  integer :: i, j, ibnd, jbnd, n, ie
  real(dp) :: det, egrid(net), pref
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
 
  do i = 1, net
     egrid(i) = etmin + (etmax-etmin)/(net-1)*(i-1)
  enddo
  eps_imag = 0.d0 
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
    
    do ibnd = 1, nbnd_occ(ik)
       do jbnd = nbnd_occ(ik)+1, nbnd
          det = (et(jbnd,ik) - et(ibnd,ik))*ry2ha   
          if (det > etmin - 8*sigma .and. det < etmax + 8*sigma) then
          ! xx xy xz yy yz zz
             n = 0
             do i = 1,3
                do j = i,3
                   n = n + 1
                   do ie = 1, net
                      !if (abs(det-egrid(ie)) < 8*sigma) then
                         pref = 4*pi**2/omega/egrid(ie)**2*wk(ik)
                         !pref = 4*pi**2/omega/det**2*wk(ik)
                         eps_imag(ie,n) = eps_imag(ie,n) + pref*CONJG(vmat(jbnd,ibnd,j))* &
                                          vmat(jbnd,ibnd,i)*delta(det-egrid(ie))
                      !endif
                   enddo
                enddo
             enddo
          endif        
       enddo
    enddo
    

  enddo ! ik
  
#ifdef __MPI
    call mp_sum(eps_imag, inter_pool_comm)
#endif 
  
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
        write(iunout, '(7F10.4)') egrid(ie)*ry2ev/ry2ha, eps_imag(ie,1),eps_imag(ie,4), &
                       eps_imag(ie,6), eps_imag(ie,2),eps_imag(ie,5),eps_imag(ie,3)
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
