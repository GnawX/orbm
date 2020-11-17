!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE calc_mag_dipole
  !-----------------------------------------------------------------------
  !
  ! This routine calculates the magnetic transition dipole matrix
  !    i/2 ( <dn|\times v|m> + <n|v \times|dm> ) in AU
  ! 
  USE kinds,                  ONLY : dp
  USE io_global,              ONLY : stdout, ionode, ionode_id
  USE io_files,               ONLY : nwordwfc, iunwfc
  USE cell_base,              ONLY : tpiba2
  USE wavefunctions,          ONLY : evc
  USE noncollin_module,       ONLY : npol, noncolin
  USE uspp,                   ONLY : nkb, vkb
  USE gvect,                  ONLY : ngm, g
  USE gvecw,                  ONLY : gcutw
  USE klist,                  ONLY : nks, nkstot, ngk, xk, igk_k 
  USE wvfct,                  ONLY : nbnd, npwx, et, current_k, g2kin
  USE lsda_mod,               ONLY : nspin, lsda, isk, current_spin
  USE orbm_module,            ONLY : ry2ha, ci, nbnd_occ
  USE buffers,                ONLY : get_buffer
  USE mp_pools,               ONLY : my_pool_id, me_pool, root_pool,  &
                                     inter_pool_comm, intra_pool_comm, npool
  USE mp,                     ONLY : mp_sum, mp_bcast
  USE mp_world,               ONLY : world_comm

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! the following three quantities are for norm-conserving PPs
  complex(dp), allocatable, dimension(:,:,:) :: vel_evc       ! v_{k,k}|evc>
  complex(dp), allocatable, dimension(:,:,:) :: evc1          ! du/dk
  ! temporary working array, same size as evc/evq
  complex(dp), allocatable :: aux(:,:)
  complex(dp), allocatable, dimension(:,:,:,:) :: ps 
  complex(dp), allocatable, dimension(:,:,:) :: ps1
  complex(dp), allocatable, dimension(:,:,:,:) :: ps2
  complex(dp), allocatable, dimension(:,:,:,:) :: mmat

  integer :: ik, ios, iunout
  integer :: i, j, ibnd, jbnd
  real(dp), external :: get_clock
  integer, external :: find_free_unit
  integer :: npw

 
  call start_clock('calc_elec_dipole')
  !-----------------------------------------------------------------------
  ! allocate memory
  !-----------------------------------------------------------------------
  allocate ( vel_evc(npwx*npol,nbnd,3), evc1(npwx*npol,nbnd,3) )
  allocate ( ps(nbnd,nbnd,3,3),  ps1(nbnd,nbnd,3), ps2(nbnd,nbnd,nks,3) )
  allocate ( aux(npwx*npol, nbnd), mmat(nbnd,nbnd,nkstot,3) )

  ! print memory estimate
  call orbm_memory_report

  write(stdout, '(5X,''Computing the magnetic dipole matrix (bohr mag):'',$)')
  write(stdout, *)


  !====================================================================
  ! loop over k-points on the pool
  !====================================================================
  do ik = 1, nks

#ifdef __MPI
    if (me_pool == root_pool) &
    write(*, '(5X,''k-point #'',I5,'' of '',I5,6X,''pool #'',I3,4X,''cpu time:'',F10.1)') &
      ik, nks, my_pool_id+1, get_clock('orbm')
#else
    write(stdout, '(5X,''k-point #'',I5,'' of '',I5,4X,''cpu time:'',F10.1)') &
      ik, nks, get_clock('orbm')
#endif

    ! initialize k, spin, g2kin used in h_psi    
    current_k = ik
    if (lsda) current_spin = isk(ik)
    npw = ngk(ik)
    call gk_sort(xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin)
    g2kin(:) = g2kin(:) * tpiba2
    call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)

    ! read wfcs from file and compute becp
    call get_buffer (evc, nwordwfc, iunwfc, ik)
    
    
    ! calculate du/dk    
    vel_evc(:,:,:) = (0.d0,0.d0)
    ps(:,:,:,:)= (0.d0,0.d0)
    ps1(:,:,:)= (0.d0,0.d0)
    
    do i = 1,3
    
       ! calculate evc1=du/dk
       call apply_vel(evc, vel_evc(1,1,i), ik, i)
       aux(:,:) = vel_evc(:,:,i)
       call greenfunction(ik, aux, evc1(1,1,i))
       
       ! calculate <n|v \times |dm>
       do j = 1,3
       
          if (noncolin) then
     
             CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), vel_evc(1,1,j), &
                    npwx*npol, evc1(1,1,i), npwx*npol, (0.d0,0.d0), ps(1,1,j,i), nbnd)
          else
       
             CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), vel_evc(1,1,j), &
                    npwx, evc1(1,1,i), npwx, (0.d0,0.d0), ps(1,1,j,i), nbnd)
          endif
          
       enddo 
    enddo
    
#ifdef __MPI
    call mp_sum(ps, intra_pool_comm)
#endif 
    
    ps1(:,:,1) = ps(:,:,2,3) - ps(:,:,3,2)
    ps1(:,:,2) = ps(:,:,3,1) - ps(:,:,1,3)
    ps1(:,:,3) = ps(:,:,1,2) - ps(:,:,2,1)
    
    do ibnd = 1, nbnd
       do jbnd = 1, nbnd
          ps2(jbnd,ibnd,ik,:) = (ps1(jbnd,ibnd,:) - CONJG(ps1(ibnd,jbnd,:)))*ci
       endo
    endo
    

  if ( npool == 1 ) then
     mmat = ps2
  else
     call poolcollect_z( nbnd, nks, ps2, nkstot, mmat)
  endif
  
  ios = 0
  if ( ionode ) then
     iunout = find_free_unit( )
     open (unit = iunout, file = 'mdipole', status = 'unknown', form = &
          'unformatted', iostat = ios)
     rewind (iunout)
  endif

  call mp_bcast (ios, ionode_id, world_comm)
  if ( ios/=0 ) call errore ('calc_mag_dipole', 'Opening file mdipole', abs (ios) )

  if (ionode) then
     write(iunout) nbnd, nkstot, nspin
     write(iunout) mmat
     close(iunout)
  endif
  
  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,*)
  write(stdout,'(5X,''End of magnetic dipole calculation'')')
  write(stdout,*)
  write(stdout,'(5X,''Matrix elements dumped in mdipole'')')
  write(stdout,*)

  ! free memory as soon as possible
  deallocate( vel_evc, evc1, aux, ps, ps1, ps2, mmat )

  
  !call restart_cleanup ( )
  call stop_clock('calc_mag_dipole')

END SUBROUTINE calc_mag_dipole
