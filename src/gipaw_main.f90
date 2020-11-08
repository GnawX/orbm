!
! Copyright (C) 2001-2009 GIPAW and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM gipaw_main
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the magnetic response program. 
  ! ... It controls the initialization routines.
  ! ... Features: NMR chemical shifts
  !...            EPR g-tensor
  ! ... Ported to Espresso by:
  ! ... D. Ceresoli, A. P. Seitsonen, U. Gerstamnn and  F. Mauri
  ! ...
  ! ... References (NMR):
  ! ... F. Mauri and S. G. Louie Phys. Rev. Lett. 76, 4246 (1996)
  ! ... F. Mauri, B. G. Pfrommer, S. G. Louie, Phys. Rev. Lett. 77, 5300 (1996)
  ! ... T. Gregor, F. Mauri, and R. Car, J. Chem. Phys. 111, 1815 (1999)
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. B 63, 245101 (2001)
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. Lett. 91, 196401 (2003)
  ! ...
  ! ... References (g-tensor):
  ! ... C. J. Pickard and F. Mauri, Phys. Rev. Lett. 88, 086403 (2002)
  ! ...
  USE kinds,           ONLY : DP
  USE io_files,        ONLY : tmp_dir, create_directory
  USE mp,              ONLY : mp_bcast
  USE cell_base,       ONLY : tpiba
  USE gipaw_module,    ONLY : job, q_gipaw, max_seconds
  USE check_stop  ,    ONLY : check_stop_init
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
  USE mp_global,       ONLY : mp_startup, nproc_pool_file
  USE mp_world,        ONLY : world_comm
  USE mp_images,       ONLY : nimage, my_image_id
  USE mp_pools,        ONLY : intra_pool_comm
  USE mp_bands,        ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_diag,         ONLY : mp_start_diag
  USE mp_bands,        ONLY : nbgrp
  USE mp_pools,        ONLY : nproc_pool
  USE environment,     ONLY : environment_start, environment_end
  USE lsda_mod,        ONLY : nspin
  USE wvfct,           ONLY : nbnd
  USE uspp,            ONLY : okvan
  USE io_global,       ONLY : stdout
  USE noncollin_module,ONLY : noncolin
  USE command_line_options, ONLY: input_file_, command_line, ndiag_
  ! for pluginization
  USE input_parameters, ONLY : nat_ => nat, ntyp_ => ntyp
  USE input_parameters, ONLY : assume_isolated_ => assume_isolated, &
                               ibrav_ => ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp
  USE cell_base,        ONLY : ibrav
  ! end
  USE gipaw_version
  USE iotk_module  
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   :: code = 'GIPAW'
  CHARACTER (LEN=10)  :: dirname = 'dummy'
  LOGICAL, EXTERNAL  :: check_para_diag
  !------------------------------------------------------------------------

  ! begin with the initialization part

  call mp_startup(start_images=.false.)
  
  ! no band parallelization
  
  call mp_start_diag(ndiag_, world_comm, intra_pool_comm, do_distr_diag_inside_bgrp_=.true.)

  call set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm)

  call environment_start (code)



  write(stdout,*)
  write(stdout,'(5X,''***** This is GIPAW git revision '',A,'' *****'')') gipaw_git_revision
  write(stdout,'(5X,''***** you can cite: N. Varini et al., Comp. Phys. Comm. 184, 1827 (2013)  *****'')')
  write(stdout,'(5X,''***** in publications or presentations arising from this work.            *****'')')
  write(stdout,*)
 

  call gipaw_readin()
  call check_stop_init( max_seconds )

  io_level = 1
 
  ! read ground state wavefunctions
  call read_file
#ifdef __MPI
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  call gipaw_openfil

  if (gamma_only) call errore ('gipaw_main', 'Cannot run GIPAW with gamma_only == .true. ', 1)
  if (okvan) call errore('gipaw_main', 'USPP not supported yet', 1)

  nat_ = nat
  ntyp_ = ntyp
  ibrav_ = ibrav
  assume_isolated_ = 'none'

  call gipaw_allocate()
  call gipaw_setup()
  call gipaw_summary()
  
  ! convert q_gipaw into units of tpiba
  q_gipaw = q_gipaw / tpiba
  

  ! calculation
  
  call orbm
  
  ! print timings and stop the code
  call gipaw_closefil
  call print_clock_gipaw
  call environment_end(code)
  call stop_code( .true. )
  
  STOP
  
END PROGRAM gipaw_main

