!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE rotation
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the cd spectra
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
                                     nw, nbnd_occ, sigma, ry2ev, evc1, vel2_evc, &
                                     invmass
  USE buffers,                ONLY : get_buffer
  USE mp_pools,               ONLY : my_pool_id, me_pool, root_pool,  &
                                     inter_pool_comm, intra_pool_comm, npool
  USE mp,                     ONLY : mp_sum, mp_bcast
  USE mp_world,               ONLY : world_comm
  !USE symme,                  ONLY : symmatrix

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! the following three quantities are for norm-conserving PPs
  complex(dp), allocatable, dimension(:,:,:) :: vmat            ! <n|v|m> 
  complex(dp), allocatable, dimension(:,:,:,:) :: v2mat
  complex(dp), allocatable, dimension(:,:,:, :) :: ps 
  
  real(dp) :: rot(nw, 3, 3)
  
  integer :: ik, ios, iunout,id
  integer :: i, j, k, ibnd, jbnd, n, ie
  real(dp) :: det, wgrid(nw), pref, mel
  real(dp), external :: get_clock
  integer, external :: find_free_unit
  integer :: npw, ind(2,3)

 
  call start_clock('rotation')
  !-----------------------------------------------------------------------
  ! allocate memory
  !-----------------------------------------------------------------------
  allocate ( vmat(nbnd, nbnd, 3), v2mat(nbnd, nbnd, 3, 3), ps(nbnd, nbnd, 3, 3))
  allocate ( vel_evc(npwx*npol, nbnd, 3), evc1(npwx*npol, nbnd,3) )
  allocate ( vel2_evc(npwx*npol, nbnd, 3, 3) )
  allocate ( invmass(nbnd,nbnd,3,3) )

  ! print memory estimate
  call orbm_memory_report

  write(stdout, '(5X,''Computing the cd spectra:'',$)')
  write(stdout, *)
 
  do i = 1, nw
     wgrid(i) = wmin + (wmax-wmin)/(nw-1)*(i-1)
  enddo
  rot = 0.d0 

  ind(:,1) = (/2,3/)
  ind(:,2) = (/3,1/)
  ind(:,3) = (/1,2/)
  !====================================================================
  ! loop over k-points on the pool
  !====================================================================
  do ik = 1, nks

    npw = ngk(ik)

    ! read wfcs from file and compute becp
    call get_buffer (evc, nwordwfc, iunwfc, ik)
    
    
    ! calculate du/dk    
    vel_evc(:,:,:) = (0.d0,0.d0)
    vel2_evc(:,:,:) = (0.d0,0.d0)
    vmat(:,:,:)= (0.d0,0.d0)
    v2mat(:,:,:)= (0.d0,0.d0)
    ps(:,:,:,:) = (0.d0,0.d0)
    
    call apply_vel2(evc, vel_evc, vel2_evc, ik)
    
    do i = 1, 3
       if (noncolin) then
          CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), evc(1,1), &
                    npwx*npol, vel_evc(1,1,i), npwx*npol, (0.d0,0.d0), vmat(1,1,i), nbnd)
       else
          CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), evc(1,1), &
                    npwx, vel_evc(1,1,i), npwx, (0.d0,0.d0), vmat(1,1,i), nbnd)
       endif
       do j = 1, 3
          if (noncolin) then
             CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), evc(1,1), &
                    npwx*npol, vel2_evc(1,1,j,i), npwx*npol, (0.d0,0.d0), v2mat(1,1,j,i), nbnd)
          else
             CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), evc(1,1), &
                    npwx, vel_evc(1,1,j,i), npwx, (0.d0,0.d0), v2mat(1,1,j,i), nbnd)
          endif
       enddo
    enddo
   
    
#ifdef __MPI
    call mp_sum(vmat, intra_pool_comm)
    call mp_sum(v2mat, intra_pool_comm)
#endif

    ! calculate evc1 and invmass
    call covariant_der(ik)
    
    do i = 1,3
       do j = 1,3
       
          if (noncolin) then
     
             CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), evc1(1,1,j), &
                    npwx*npol, vel_evc(1,1,i), npwx*npol, (0.d0,0.d0), ps(1,1,j,i), nbnd)
          else
       
             CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), evc1(1,1,j), &
                    npwx, vel_evc(1,1,i), npwx, (0.d0,0.d0), ps(1,1,j,i), nbnd)
             
          endif
          
       enddo 
    enddo
   
    
#ifdef __MPI
    call mp_sum(ps, intra_pool_comm)
#endif    
    
    do ibnd = 1, nbnd_occ(ik)
       do jbnd = nbnd_occ(ik)+1, nbnd
          do i = 1, 3
             do j = 1, 3
       
                ps(jbnd,ibnd,j,i) = invmass(jbnd,ibnd,j,i) - v2mat(jbnd,ibnd,j,i) - CONJG(ps(ibnd,jbnd,j,i))
                
             enddo
          enddo
       enddo
    enddo

    do ie = 1, nw
       pref = 4*pi**2/omega/137*wk(ik)
       
       do id = 1, 3
          i = ind(1,id); j = ind(2,id)
          !do j = 1, 3
             do k = 1, 3
          
                do ibnd = 1, nbnd_occ(ik)
                   do jbnd = nbnd_occ(ik)+1, nbnd
                
                      det = (et(jbnd,ik) - et(ibnd,ik))*ry2ha
                      mel = REAL( vmat(ibnd,jbnd,j)*( ps(jbnd,ibnd,k,i) - CONJG(ps(ibnd,jbnd,k,i)) )  - &
                            vmat(ibnd,jbnd,i)*( ps(jbnd,ibnd,k,j) - CONJG(ps(ibnd,jbnd,k,j)) ), DP)
                      !print *, mel
                      !if (det > wmin - 8*sigma .and. det < wmax + 8*sigma &
                      !    .and. ABS(det-wgrid(ie)) < 8*sigma) then
                       
                      rot(ie,id,k) = rot(ie,id,k) + pref*mel*sigma/pi/((det-wgrid(ie))**2 &
                                      + sigma**2)/det
                     ! rot(ie,i,j,k) = rot(ie,i,j,k) + pref*mel*sigma/pi/((det-wgrid(ie))**2 &
                     !                 + sigma**2)/det
                      !endif
                   enddo
                enddo

             enddo
          !enddo
       enddo
       
    enddo
                    

  enddo ! ik
  
#ifdef __MPI
  call mp_sum(rot, inter_pool_comm)
#endif 


  !do ie = 1, nw
  !   call sym3tensor( rot(ie, :, :, :) )
  !enddo
  do ie = 1, nw
     call symaxialtensor2( rot(ie, :, :)  )
  enddo

  ios = 0
  if ( ionode ) then
     iunout = find_free_unit( )
     open (unit = iunout, file = 'rotation.dat', status = 'unknown', iostat = ios)
     rewind (iunout)
  endif

  call mp_bcast (ios, ionode_id, world_comm)
  if ( ios/=0 ) call errore ('rotation', 'Opening file rotatory.dat', abs (ios) )

  if (ionode) then
     write(iunout, '(10A10)') '#   ENERGY','YZX','YZY','YZZ','ZXX','ZXY','ZXZ','XYX', 'XYY', 'XYZ' 
     do ie = 1, nw
        !write(iunout, '(F10.4,9E12.4)') wgrid(ie)*ry2ev/ry2ha, rot(ie,2,3,1),rot(ie,3,1,2),rot(ie,1,2,3), &
        !               rot(ie,2,3,2),rot(ie,3,1,3),rot(ie,1,2,1),rot(ie,2,3,3),rot(ie,3,1,1),rot(ie,1,2,2)
        write(iunout, '(F10.4,9E12.4)') wgrid(ie)*ry2ev/ry2ha, ((rot(ie,i,j),j=1,3),i=1,3)  
        !write(iunout, '(F10.4,3E12.4)') wgrid(ie)*ry2ev/ry2ha, rot(ie,2,3,1),rot(ie,3,2,1),rot(ie,1,1,1)
     enddo
     close(iunout)
  endif


  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,*)
  write(stdout,'(5X,''End of CD calculation'')')
  write(stdout,*)


  ! free memory as soon as possible
  deallocate( vel_evc, evc1, vel2_evc, invmass, vmat, v2mat, ps)

  
  !call restart_cleanup ( )
  call stop_clock('rotation')

END SUBROUTINE rotation

