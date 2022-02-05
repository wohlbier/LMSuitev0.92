!**************************************************************************!
! This file is part of LATTE/MUSE Numerical Suite (lmsuite)
!
! lmsuite is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! lmsuite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with lmsuite; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!**************************************************************************!
subroutine solveTWT
  use parameters
  use working_variables, only : modelID
  use working_variables, only : M
  use working_variables, only : run_one, run_two, run_three
  use working_variables, only : derived_qtys
  !** modules below share routines commonTWT **!
  !** model specific information also in module **!
  use LATTE_module, only : init_latte, kill_latte
  use LATTE_module, only : solvelatte => commonTWT
  use MUSE_module, only : init_muse, kill_muse
  use MUSE_module, only : solvemuse => commonTWT
  use SMUSE_module, only : init_smuse, kill_smuse
  use SMUSE_module, only : solvesmuse => commonTWT

  implicit none

  !** solve LATTE **!
  if (run_parameters % select_code == 'L' &
       .or. run_parameters % select_code == 'LM' &
       .or. run_parameters % select_code == 'LS' &
       .or. run_parameters % select_code == 'LMS') then

     if (run_two) then
        print*, 'Running LATTE'
     end if

     !** set modelID to LATTE **!
     modelID = 'L'
     call init_latte
     !within solvelatte are the scans, unnormalization and plotting
     call solvelatte
     call kill_latte
  end if

  !** solve MUSE **!
  if (run_parameters % select_code == 'M' &
       .or. run_parameters % select_code == 'LM' &
       .or. run_parameters % select_code == 'MS' &
       .or. run_parameters % select_code == 'LMS') then

     if (run_two) then
        print*, 'Running MUSE'
     end if

     !** set modelID to MUSE **!
     modelID = 'M'
     call init_muse
     call solvemuse
     !within solvelatte are the scans, unnormalization and plotting
     call kill_muse
  end if

  !** solve S-MUSE **!
  if (run_parameters % select_code == 'S' &
       .or. run_parameters % select_code == 'MS' &
       .or. run_parameters % select_code == 'LS' &
       .or. run_parameters % select_code == 'LMS') then

     if (run_two) then
        print*, 'Running S-MUSE'
     end if

     !** set modelID to S-MUSE **!
     modelID = 'S'
     call init_smuse
     call solvesmuse
     !within solvelatte are the scans, unnormalization and plotting
     call kill_smuse
  end if
end subroutine solveTWT
