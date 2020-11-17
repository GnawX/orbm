!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE orbm_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for GIPAW calculations
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : a0_to_cm => bohr_radius_cm
  USE parameters, ONLY : npk, ntypx
  
  IMPLICIT NONE
  SAVE


  ! Physical constants:
  ! fine structure constant alpha (c = 1/alpha)
  REAL(DP), PARAMETER :: alpha = 1.0_dp / 137.03599911_dp
  
  ! avogadro number
  REAL(DP), PARAMETER :: avogadro = 6.022142e23_dp
  
  ! g_e and gprime
  REAL(DP), PARAMETER :: g_e = 2.0023192778_DP
  REAL(DP), PARAMETER :: gprime = 2.d0 * (g_e - 1.d0)
 
  ! rydberg to Hartree
  REAL(DP), PARAMETER :: ry2ha = 0.5_DP
  
  ! imaginary unit
  COMPLEX(DP), PARAMETER :: ci = (0.0_DP, 1.0_DP)
 
  ! number of occupied bands at each k-point
  INTEGER :: nbnd_occ(npk)
  
  ! alpha shift of the projector on the valence wfcs
  REAL(DP) :: alpha_pv

  ! convergence threshold for diagonalizationa and greenfunction
  REAL(DP) :: conv_threshold
  
  ! q for the perturbation (in bohrradius^{-1})
  REAL(DP) :: q_orbm
  

  ! restart mode: 'from_scratch' or 'restart'
  CHARACTER(80) :: restart_mode

  ! max CPU time, in s
  REAL(dp) :: max_seconds
  
  ! verbosity
  INTEGER :: iverbosity
 
  ! diagonalization method
  INTEGER :: isolve
 
  ! job: nmr, g_tensor, efg, hyperfine
  CHARACTER(80) :: job

 
  

!-----------------------------------------------------------------------
END MODULE orbm_module
!-----------------------------------------------------------------------
