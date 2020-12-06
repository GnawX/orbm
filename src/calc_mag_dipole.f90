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
  !        M_nm = ie/4 (<dn|\times v|m> + <n|v \times|dm> )
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
  USE orbm_module,            ONLY : ry2ha, ci, nbnd_occ, dudk_method, vel_evc, evc1
  USE buffers,                ONLY : get_buffer
  USE mp_pools,               ONLY : my_pool_id, me_pool, root_pool,  &
                                     inter_pool_comm, intra_pool_comm, npool
  USE mp,                     ONLY : mp_sum, mp_bcast
  USE mp_world,               ONLY : world_comm

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE

  ! temporary working array, same size as evc/evq
  complex(dp), allocatable :: aux(:,:)
  complex(dp), allocatable, dimension(:,:,:,:) :: ps
  complex(dp), allocatable, dimension(:,:,:,:) :: ps1
  complex(dp), allocatable, dimension(:,:,:,:) :: ps2
  complex(dp), allocatable, dimension(:,:,:,:) :: ps3
  complex(dp), allocatable, dimension(:,:,:,:) :: mmat

  integer :: ik, ios, iunout
  integer :: i, j, ibnd, jbnd
  real(dp), external :: get_clock
  integer, external :: find_free_unit
  integer :: npw

 
  call start_clock('calc_mag_dipole')
  
  !-----------------------------------------------------------------------
  ! allocate memory
  !-----------------------------------------------------------------------
  allocate ( vel_evc(npwx*npol, nbnd, 3), evc1(npwx*npol, nbnd, 3) )
  allocate ( ps1(nbnd,nbnd,3,3), ps2(nbnd, nbnd, 3,3), ps3(nbnd,nbnd,nks,3) ) !, ps(nbnd,nbnd,nks,3) )
  allocate ( aux(npwx*npol, nbnd), mmat(nbnd,nbnd,nkstot,3) )

  ! print memory estimate
  call orbm_memory_report

  write(stdout, '(5X,''Computing the magnetic dipole matrix (bohr mag):'',$)')
  write(stdout, *)


  !====================================================================
  ! loop over k-points on the pool
  !====================================================================
  ps3 = (0.d0,0.d0)
  mmat = (0.d0,0.d0)

  do ik = 1, nks

#ifdef __MPI
    if (me_pool == root_pool) &
    write(*, '(5X,''k-point #'',I5,'' of '',I5,6X,''pool #'',I3,4X,''cpu time:'',F10.1)') &
      ik, nks, my_pool_id+1, get_clock('orbm')
#else
    write(stdout, '(5X,''k-point #'',I5,'' of '',I5,4X,''cpu time:'',F10.1)') &
      ik, nks, get_clock('orbm')
#endif

    npw = ngk(ik)

    ! read wfcs from file 
    call get_buffer (evc, nwordwfc, iunwfc, ik)
    
    ! calculate velocity operator
    vel_evc(:,:,:) = (0.d0,0.d0)
    do i = 1,3
       call apply_vel(evc, vel_evc(1,1,i), ik, i)
    enddo
    
    ! calculate du/dk
    call compute_dudk(ik)
    
    !vel_evc(:,:,:) = (0.d0,0.d0)
    !do i = 1,3
    !   call apply_vel_NL(evc, vel_evc(1,1,i), ik, i)
    !enddo
       !do i = 1,3 
       !! transition dipole matrix only accurate for valence to conduction transitions
       !   if (noncolin) then
    ! 
     !        CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), evc(1,1), &
      !              npwx*npol, evc1(1,1,i), npwx*npol, (0.d0,0.d0), ps(1,1,ik,i), nbnd)
      !    else
       
      !       CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), evc(1,1), &
      !              npwx, evc1(1,1,i), npwx, (0.d0,0.d0), ps(1,1,ik,i), nbnd)
             
      !    endif
      ! enddo

    
    ps1(:,:,:,:)= (0.d0,0.d0)
    ps2(:,:,:,:)= (0.d0,0.d0)
       ! calculate <dpsi_v | \times v|psi_c>
    do i = 1,3
       do j = 1,3
       
          if (noncolin) then
     
             CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), evc1(1,1,j), &
                    npwx*npol, vel_evc(1,1,i), npwx*npol, (0.d0,0.d0), ps1(1,1,j,i), nbnd)
             !CALL zgemm('C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), vel_evc(1,1,j), &
             !       npwx*npol, evc1(1,1,i), npwx*npol, (0.d0,0.d0), ps2(1,1,j,i), nbnd)
          else
       
             CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), evc1(1,1,j), &
                    npwx, vel_evc(1,1,i), npwx, (0.d0,0.d0), ps1(1,1,j,i), nbnd)
             !CALL zgemm('C', 'N', nbnd, nbnd, npw, (1.d0,0.d0), vel_evc(1,1,j), &
             !       npwx, evc1(1,1,i), npwx, (0.d0,0.d0), ps2(1,1,j,i), nbnd)
             
          endif
          
       enddo 
    enddo
    
#ifdef __MPI
    !call mp_sum(ps, intra_pool_comm)
    call mp_sum(ps1, intra_pool_comm)
    !call mp_sum(ps2, intra_pool_comm)
#endif 
    !ps2 = CONJG(ps1)
    ! in unit of Bohr mag 
    do ibnd = 1, nbnd
       do jbnd = 1, nbnd
          ps3(jbnd,ibnd,ik,1) = (ps1(jbnd,ibnd,2,3)-ps1(jbnd,ibnd,3,2) - &
                                 CONJG(ps1(ibnd,jbnd,2,3)-ps1(ibnd,jbnd,3,2)))/2*ci 
          ps3(jbnd,ibnd,ik,2) = (ps1(jbnd,ibnd,3,1)-ps1(jbnd,ibnd,1,3) - &
                                 CONJG(ps1(ibnd,jbnd,3,1)-ps1(ibnd,jbnd,1,3)))/2*ci 
          ps3(jbnd,ibnd,ik,3) = (ps1(jbnd,ibnd,1,2)-ps1(jbnd,ibnd,2,1) - &
                                 CONJG(ps1(ibnd,jbnd,1,2)-ps1(ibnd,jbnd,2,1)))/2*ci
       enddo
    enddo 
    !ps3(:,:,ik,1) = ( ps1(:,:,2,3) - ps1(:,:,3,2) + ps2(:,:,2,3) - ps2(:,:,3,2) )/2 *ci 
    !ps3(:,:,ik,2) = ( ps1(:,:,3,1) - ps1(:,:,1,3) + ps2(:,:,3,1) - ps2(:,:,1,3) )/2 *ci 
    !ps3(:,:,ik,3) = ( ps1(:,:,1,2) - ps1(:,:,2,1) + ps2(:,:,1,2) - ps2(:,:,2,1) )/2 *ci 

    !ps3(:,:,ik,1) = ( ps1(:,:,2,3)  - ps2(:,:,3,2) ) *ci 
    !ps3(:,:,ik,2) = ( ps1(:,:,3,1)  - ps2(:,:,1,3) ) *ci 
    !ps3(:,:,ik,3) = ( ps1(:,:,1,2)  - ps2(:,:,2,1) ) *ci 
    
    !do ibnd = nbnd_occ(ik)+1, nbnd
    !   do jbnd = 1, nbnd_occ(ik)
    !      ps2(ibnd,jbnd,ik,:) = CONJG(ps2(jbnd,ibnd,ik,:))
    !   enddo
    !enddo
  
    ! these parts are not correct and not needed
    !ps2(1:nbnd_occ(ik),1:nbnd_occ(ik),ik,:) = (0.d0,0.d0) 
    !ps2(nbnd_occ(ik)+1:nbnd,nbnd_occ(ik)+1:nbnd,ik,:) = (0.d0,0.d0) 
    !
  enddo !ik
  
  if ( npool == 1 ) then
     mmat = ps3
  else
     call poolcollect_z( nbnd, nks, ps3, nkstot, mmat)
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
  deallocate( vel_evc, evc1, aux, ps1, ps2, mmat )

  
  !call restart_cleanup ( )
  call stop_clock('calc_mag_dipole')

END SUBROUTINE calc_mag_dipole
