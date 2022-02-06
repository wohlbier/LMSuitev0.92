!***************************************************************************!
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
!***************************************************************************!



!**************************************************************************!
!**  external subroutine to check for existence of input namelist, read  **!
!**  the namelist
!**************************************************************************!
subroutine process_inputs
  !** use modules **!
  use parameters

  !** no implicit variable naming **!
  implicit none

  !** required variables **!
  character junk !used to process <return> from keyboard
  logical :: file_stat !variable to check existence of input file

  !**** get command line arguments ****!
  if (command_argument_count() /= 1) then
     print*, 'Usage: ./bin/lmsuite /path/to/inputs'
     stop
  endif
  call get_command_argument(1, path_to_inputs_directory)

  path_to_lmsuite_nml = trim(path_to_inputs_directory)//"lmsuite.nml"

  !**** check for the file lmsuite.nml and read the namelist ****!
  inquire(file=path_to_lmsuite_nml, exist=file_stat)
  if (file_stat) then
     call read_parameters
  else
     !! if lmsuite.nml doesn't exist, set the default values
     call set_default_parameters !! still need to check that this works!!
  end if

  !****  print details of LMSuite  ****!
  call about_lmsuite
10 format(A)
20 format(/,/,/,/,/,/,/,/)
  if (.not.file_stat) then
     print 20
     print*, './inputs/lmsuite.nml does not exist, running default example'
     print 20
     print*, '(press enter to continue)'
     read 10, junk
  end if

!** if interface_parameters % echo_namelists is true, then echo to screen **!
  if (interface_parameters % echo_namelists) then
     call echo_parameters
  end if



contains
  !****  subroutine to print details of LMSuite  ****!
  subroutine about_lmsuite
    character junk
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, 'LATTE/MUSE Numerical Suite v0.92.'
    print*, 'Copyright 2003 John G. Wohlbier'
    print*, 'This software is protected under the Gnu General Public License (GPL)'
    print*, 'For the text of the GPL, see the file COPYING or go to'
    print*, 'http://www.gnu.org/licenses/licenses.html'
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
30  format(A)
    if (interface_parameters % with_pauses) then
       print*, '(press enter to continue)'
       read 30, junk
    end if
  end subroutine about_lmsuite
  
end subroutine process_inputs
