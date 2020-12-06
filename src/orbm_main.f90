!
! Copyright (C) 2001-2009 GIPAW and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM orbm_main
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the orbm program. 
  ! ... It controls the initialization routines.
  ! ... Features: electric transition dipole
  !...            orbital magnetization
  ! ...
  ! ...
  USE kinds,           ONLY : DP
  USE mp,              ONLY : mp_bcast
  USE cell_base,       ONLY : tpiba
  USE orbm_module,     ONLY : job, q_orbm, max_seconds, emin, emax, sigma, &
                              ry2ev
  USE check_stop  ,    ONLY : check_stop_init
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
  USE mp_global,       ONLY : mp_startup, nproc_pool_file
  USE mp_world,        ONLY : world_comm
  USE mp_pools,        ONLY : intra_pool_comm
  USE mp_bands,        ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_bands,        ONLY : nbgrp
  USE mp_pools,        ONLY : nproc_pool
  USE environment,     ONLY : environment_start, environment_end
  USE lsda_mod,        ONLY : nspin
  USE wvfct,           ONLY : nbnd
  USE uspp,            ONLY : okvan
  USE io_global,       ONLY : stdout

  USE command_line_options, ONLY: input_file_, command_line, ndiag_
  ! for pluginization
  USE input_parameters, ONLY : nat_ => nat, ntyp_ => ntyp
  USE input_parameters, ONLY : assume_isolated_ => assume_isolated, &
                               ibrav_ => ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp
  USE cell_base,        ONLY : ibrav
  ! end
  USE iotk_module  
  !------------------------------------------------------------------------
  IMPLICIT NONE

  include 'laxlib.fh'

  CHARACTER (LEN=9)   :: code = 'orbm'
  CHARACTER (LEN=10)  :: dirname = 'dummy'
  LOGICAL, EXTERNAL  :: check_para_diag
  !------------------------------------------------------------------------

  ! begin with the initialization part

  call mp_startup(start_images=.false.)
  
  ! no band parallelization
  
  call laxlib_start(ndiag_, world_comm, intra_pool_comm, do_distr_diag_inside_bgrp_=.true.)

  call set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm)

  call environment_start (code)



  write(stdout,*)
  write(stdout,'(5X,''***** This is orbm git revision '',A,'' *****'')') 'beta version'
  write(stdout,'(5X,''***** you can cite: ------  *****'')')
  write(stdout,'(5X,''***** in publications or presentations arising from this work.            *****'')')
  write(stdout,*)
 

  call orbm_readin()
  call check_stop_init( max_seconds )

  io_level = 1
 
  ! read ground state wavefunctions
  call read_file
#ifdef __MPI
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  call orbm_openfil

  if (gamma_only) call errore ('orbm_main', 'Cannot run orbm with gamma_only == .true. ', 1)
  if (okvan) call errore('orbm_main', 'USPP not supported yet', 1)

  nat_ = nat
  ntyp_ = ntyp
  ibrav_ = ibrav
  assume_isolated_ = 'none'

  call orbm_allocate()
  call orbm_setup()
  call orbm_summary()
  
  ! convert q_orbm into units of tpiba
  q_orbm = q_orbm / tpiba
  
  emin = emin / ry2ev / 2.d0
  emax = emax / ry2ev / 2.d0
  sigma = sigma /ry2ev/ 2.d0

  ! calculation
  select case ( trim(job) )
  case ( 'orbm' )
     call calc_orb_magnetization 
  case ( 'elec_dipole' )
     call calc_elec_dipole
  case ( 'epsilon' )
     call epsilon
  case ( 'mag_dipole' )
     call calc_mag_dipole
  case ( 'elec_quadrupole' )
     !call calc_elec_quadrupole
  case default
     call errore('orbm_main', 'wrong or undefined job in input', 1)
  end select
  
  ! print timings and stop the code
  call orbm_closefil
  call print_clock_orbm
  call environment_end(code)
  call stop_code( .true. )
  
  STOP
  
END PROGRAM orbm_main

