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
!! routine runs two passes, the first to determine the dc content of
!! the beam, the second with the dc removed
subroutine make_frames(plot_cube, n_eqns)
  use precision95
  use parameters, only : pi
  use parameters, only : smo
  use parameters, only : eps0
  use parameters, only : run_parameters
  use parameters, only : output_parameters
  use parameters, only : circuit_parameters
  use parameters, only : units_structure
  use parameters, only : numerical_parameters
  use parameters, only : beam_parameters
  use parameters, only : frequency_parameters 
  use working_variables, only : derived_qtys
  use working_variables, only : w0
  use working_variables, only : modelID
  use working_variables, only : run_two
  use working_variables, only : M
  use working_variables, only : fl
  use working_variables, only : first_ckt_index
  use working_variables, only : last_ckt_index
  use working_variables, only : num_ckt_freqs
  use working_variables, only : first_spch_index
  use working_variables, only : last_spch_index
  use working_variables, only : num_spch_freqs
  use working_variables, only : run_three
  use working_variables, only : N_disks
  use working_variables, only : K
  use working_variables, only : pC
  use working_variables, only : vph
  use working_variables, only : find_freq
  
  integer, intent(in) :: n_eqns
  complex(dp), dimension(1,output_parameters % num_axial_points,n_eqns), &
       intent(in) :: plot_cube
  

  character(len=30) :: filename
  real(dp), dimension(:,:), allocatable :: mag_phase
  integer :: i,j,kk,l
  real(dp) :: z, dz, t, dt,amplitude
  
  if (run_two) then
     print*, '  making frames for animation'
  end if

  !! check to make sure just running one model
  if (.not. (run_parameters % select_code == 'L' .or. &
       run_parameters % select_code == 'M' .or. &
       run_parameters % select_code == 'S')) then
     print*, '**********************************************************'
     print*, 'WARNING: When running multiple models, frames will only be'
     print*, 'generated for the last model run.'
     print*, '**********************************************************'
  end if

  
  dz = circuit_parameters % circuit_length &
       / (output_parameters % num_axial_points-1)
  
  do i = 1, run_parameters % num_movie_namelists
     !! compute delta t
     dt = run_parameters % movie_data_array(i) % time &
          / (run_parameters % movie_data_array(i) % num_frames - 1)
     t = 0.0
     
     !! allocate the mag_phase matrix
     if (run_parameters % movie_data_array(i) % movie_type < 2) then
        allocate(mag_phase(output_parameters % num_axial_points, &
             2*num_ckt_freqs))
     else
        allocate(mag_phase(output_parameters % num_axial_points, &
             2*num_spch_freqs))
     end if
     
     !! compute magnitudes and phases
     call compute_mag_phase(mag_phase, &
          run_parameters % movie_data_array(i) % movie_type)
     
     if (valid_movie(run_parameters % movie_data_array(i) % movie_type) .and. &
          run_parameters % movie_data_array(i) % movie_type /= 7) then
        !! loop on num_frames
        do j = 1, run_parameters % movie_data_array(i) % num_frames
           z = 0.0
           !! create the file for frame j in movie_i
           call create_file(i,j)

           !! compute the wave value
           do kk = 1, output_parameters % num_axial_points
              amplitude = 0.0
              !! loop on the frequencies to make the f(z,t) value
              if (run_parameters % movie_data_array(i) % movie_type <= 2) then
                 !! circuit voltage or current
                 do l = first_ckt_index, last_ckt_index
                    if (run_parameters % svea) then
                       !! use vph for reference wave
                       amplitude=amplitude+2*mag_phase(kk,2*l-1)*cos(2*pi*&
                            fl(l)*frequency_parameters%base_frequency*&
                            (z/vph(z,l)-t)+mag_phase(kk,2*l))*exp(&
                            -run_parameters%movie_data_array(i)%scale*z)
                    else
                       !! use u0 for reference wave
                       amplitude=amplitude+2*mag_phase(kk,2*l-1)*cos(2*pi*&
                            fl(l)*frequency_parameters%base_frequency*&
                            (z/derived_qtys%u0-t)+mag_phase(kk,2*l))*exp(&
                            -run_parameters%movie_data_array(i)%scale*z)
                    end if
                 end do
              else
                 !! beam quantities
                 do l = first_spch_index, last_spch_index
                    !! use u0 for reference wave
                    amplitude=amplitude+2*mag_phase(kk,2*l-1)*&
                         cos(2*pi*fl(l)*frequency_parameters%base_frequency*&
                         (z/derived_qtys%u0 - t)+mag_phase(kk,2*l))*&
                         exp(-run_parameters%movie_data_array(i)%scale*z)
                 end do
              end if
              
              !! write the value to the file
              !! note: if units_structure % length_cm then write z*100.0
              if (units_structure % length_cm) then
                 write(1, "(en20.8, en20.8)") z*100.0, amplitude
              else
                 write(1, "(en20.8, en20.8)") z, amplitude
              endif
              
              z = z + dz
           end do
           
           close(1)
           t = t + dt
        end do
     else if (valid_movie(run_parameters % movie_data_array(i) % movie_type) &
          .and. run_parameters % movie_data_array(i) % movie_type == 7) then
        !phase space movie
        call phase_space(plot_cube, i)
     end if
     deallocate(mag_phase)
  end do
contains
  
  function valid_movie(movie_type)
    logical :: valid_movie
    integer, intent(in) :: movie_type

    valid_movie = .true.

    !! there are certain movie types that are not valid or are not implemented
    !! if one of these are chosen, do not generate frames
    if (modelID == 'L' .and. movie_type == 4) then
       print*, '*********************************************************'
       print*, '  WARNING: Velocity wave forms not implemented for LATTE.'
       print*, '  Not generating LATTE velocity frames.'
       print*, '*********************************************************'
       valid_movie = .false.
    end if

    if (run_parameters % svea .and. movie_type == 2) then
       print*, '*********************************************************'
       print*, '  WARNING: Circuit current not computed with SVEA = TRUE.'
       print*, '  Not generating circuit current frames.'
       print*, '*********************************************************'
       valid_movie = .false.
    end if

    if ((modelID == 'M' .or. modelID == 'S') .and. movie_type == 7) then
       print*, '***********************************************************'
       print*, '  WARNING: Phase space plots are only available from LATTE.'
       print*, '  Not generating LATTE velocity frames.'
       print*, '***********************************************************'
       valid_movie = .false.
    end if

  end function valid_movie

  subroutine compute_mag_phase(mag_phase, movie_type)
    real(dp), dimension(:,:), intent(inout) :: mag_phase
    integer, intent(in) :: movie_type
    
    complex(dp), dimension(:), allocatable :: tmpY
    complex(dp) :: rhol, ibeaml
    integer :: i,j,kk,l
    
    !! pick out a vector
    allocate(tmpY(n_eqns))
    
    !! loop on number axial points
    do i = 1, output_parameters % num_axial_points
       tmpY = plot_cube(1,i,:)
       
       !! unnormalize it
       if (modelID == 'L') then
          call unnorm_latte_vector(tmpY)
       else
          call unnorm_muse_vector(tmpY)
       end if
       
       if (movie_type == 1) then
          !! circuit voltage
          !! put mags and phases into mag_phase
          do j = first_ckt_index, last_ckt_index
             mag_phase(i,2*j-1) = abs(tmpY(j))
             mag_phase(i,2*j) = atan2(aimag(tmpY(j)),real(tmpY(j)))
          end do
       else if (movie_type == 2) then
          !! circuit current
          !! put mags and phases into mag_phase
          do j = first_ckt_index, last_ckt_index
             if (modelID=='L') then
                mag_phase(i,2*j-1) = abs(tmpY(j+M))
                mag_phase(i,2*j) = atan2(aimag(tmpY(j+M)),real(tmpY(j+M)))
             else
                mag_phase(i,2*j-1) = abs(tmpY(j+M+1))
                mag_phase(i,2*j) = atan2(aimag(tmpY(j+M+1)),real(tmpY(j+M+1)))
             end if
          end do
       else if (movie_type == 3) then
          !! space charge
          !! put mags and phases into mag_phase
          do j = first_spch_index, last_spch_index
             if (modelID=='L') then
                mag_phase(i,2*j-1) = abs(tmpY(j+2*M))
                mag_phase(i,2*j) = atan2(aimag(tmpY(j+2*M)),real(tmpY(j+2*M)))
             else
                mag_phase(i,2*j-1) = abs(tmpY(j+2*(M+1)))
                mag_phase(i,2*j) = &
                     atan2(aimag(tmpY(j+2*(M+1))), real(tmpY(j+2*(M+1))))
             end if
          end do
       else if (movie_type == 4) then
          !! velocity
          if (modelID == 'L') then
             !! need to add this function
          else
             do j = first_spch_index, last_spch_index
                mag_phase(i,2*j-1) = abs(tmpY(j+3*(M+1)))
                mag_phase(i,2*j) = &
                     atan2(aimag(tmpY(j+3*(M+1))), real(tmpY(j+3*(M+1))))
             end do
          end if
       else if (movie_type == 5) then
          !! density
          do j = first_spch_index, last_spch_index
             if (modelID == 'L') then
                !! density
                !! compute rhol
                rhol = (0.0,0.0)
                do l = 3*M+1, 3*M+N_disks
                   !! unnormalized, no w0
                   rhol = &
                        rhol + exp(-smo*fl(j)*tmpY(l+N_disks))/tmpY(l)
                end do
                rhol = derived_qtys % rho0 * derived_qtys % u0 &
                     * 1.0/real(N_disks) * rhol
                !! assign the magnitude and phase
                mag_phase(i,2*j-1) = abs(rhol)
                mag_phase(i,2*j) = atan2(aimag(rhol),real(rhol))
             else
                mag_phase(i,2*j-1) = abs(tmpY(j+4*(M+1)))
                mag_phase(i,2*j) = &
                     atan2(aimag(tmpY(j+4*(M+1))), real(tmpY(j+4*(M+1))))
             end if
          end do
       else if (movie_type == 6) then
          !! beam current
          do j = first_spch_index, last_spch_index
             if (modelID == 'L') then
                !! compute ibeaml
                ibeaml = (0.0,0.0)
                do l = 3*M+1, 3*M+N_disks
                   !! unnormalized, no w0
                   ibeaml = &
                        ibeaml + exp(-smo*fl(j)*tmpY(l+N_disks))
                end do
                ibeaml = derived_qtys % rho0 * derived_qtys % u0 &
                     * derived_qtys % beam_area &
                     * 1.0/real(N_disks) * ibeaml
                mag_phase(i,2*j-1) = abs(ibeaml)
                mag_phase(i,2*j) = atan2(aimag(ibeaml),real(ibeaml))
             else
                !! compute beam current for MUSE, S-MUSE
                ibeaml = (0.0,0.0)
                !! j = frequency index, loop on other frequencies
                do l = -last_spch_index, last_spch_index
                   !!f_m /= 0 and f_n /= 0
                   if (l /= j .and. l /= 0 .and. &
                        abs(l) >= first_spch_index) then
                      kk = find_freq(fl(j)-fl(l))
                      if (abs(kk) >= first_spch_index) then
                         if (l < 0 .and. kk > 0) then
                            !!f_m < 0 and f_n > 0
                            ibeaml = ibeaml &
                                 + conjg(tmpY(abs(l)+4*(M+1))) &
                                 *tmpY(kk+3*(M+1))
                         else if (l > 0 .and. kk < 0) then
                            !!f_m > 0 and f_n < 0
                            ibeaml = ibeaml &
                                 + tmpY(l+4*(M+1)) &
                                 *conjg(tmpY(abs(kk)+3*(M+1)))
                         else if (l > 0 .and. kk > 0) then
                            !!f_m > 0, f_n > 0
                            ibeaml = ibeaml + tmpY(l+4*(M+1)) &
                                 *tmpY(kk+3*(M+1))
                         end if
                      end if
                   else if (l == j .and. &
                        l >= first_spch_index) then
                      !!f_m = f_\ell /= 0 (since i /= M+1), f_n = 0
                      ibeaml = ibeaml + tmpY(j+4*(M+1))*tmpY(4*(M+1))
                   else if (l == 0) then
                      !! f_m = 0, f_n = f_\ell /= 0 (since i /= M+1)
                      ibeaml = ibeaml + tmpY(5*(M+1))*tmpY(j+3*(M+1))
                   end if
                end do
                ibeaml = ibeaml * derived_qtys % beam_area

                mag_phase(i,2*j-1) = abs(ibeaml)
                mag_phase(i,2*j) = atan2(aimag(ibeaml),real(ibeaml))
             end if
          end do
       end if
    end do
    deallocate(tmpY)
  end subroutine compute_mag_phase
  
  subroutine unnorm_latte_vector(y)
    complex(dp), dimension(:), intent(inout) :: y
    integer i,j

    !unnormalize voltage
    do i = 1, M !!loop on frequencies, pick out drives and assign voltage
       !positive frequency
       if (i >= first_ckt_index .and. i <= last_ckt_index) then
	  y(i) = y(i) * K(0.0d0,i)*beam_parameters % current / pC(0.0d0,i)
       end if
    end do
    
    !unnormalize circuit current
    do i = 1, M !!loop on frequencies
       if (i >= first_ckt_index .and. i <= last_ckt_index) then
          y(i+M) = y(i+M) * beam_parameters % current / pC(0.0d0,i)
       end if
    end do

    !unnormalize space charge field
    do i = 1, M !!loop on frequencies
       if (i >= first_spch_index .and. i <= last_spch_index) then
	  y(i+2*M) = y(i+2*M) * circuit_parameters % circuit_length &
               * derived_qtys % rho0 / eps0
       end if
    end do

    !unnormalize disk velocities
    do i = 3*M + 1, 3*M + N_disks
       y(i) = y(i) * derived_qtys % u0
    end do
    
    !unnormalize disk phases
    do i = 3*M + N_disks + 1, 3*M + 2*N_disks
       y(i) = y(i) * w0
    end do
  end subroutine unnorm_latte_vector

  subroutine unnorm_muse_vector(y)
    complex(dp), dimension(:), intent(inout) :: y
    integer i,j

    !unnormalize voltage
    do i = 1, M !!loop on frequencies, pick out drives and assign voltage
       !positive frequency
       if (i >= first_ckt_index &
            .and. i <= first_ckt_index + num_ckt_freqs - 1) then
	  y(i) = y(i) * K(0.0d0,i)*beam_parameters % current / pC(0.0d0,i)
       end if
    end do
    
    !unnormalize circuit current
    do i = 1, M !!loop on frequencies
       if (i >= first_ckt_index &
            .and. i <= first_ckt_index + num_ckt_freqs - 1) then
          y(i+M+1) = y(i+M+1) * beam_parameters % current / pC(0.0d0,i)
       end if
    end do

    !unnormalize space charge field
    do i = 1, M !!loop on frequencies
       if (i >= first_spch_index &
            .and. i <= first_spch_index + num_spch_freqs - 1) then
	  y(i+2*(M+1)) = y(i+2*(M+1)) * circuit_parameters % circuit_length &
               * derived_qtys % rho0 / eps0
       end if
    end do

    !unnormalize velocity
    do i = 1, M !!loop on frequencies
       y(i+3*(M+1)) = y(i+3*(M+1)) * derived_qtys % u0
    end do
    !! dc portion
    y(4*(M+1)) = y(4*(M+1)) * derived_qtys % u0

    !unnormalize density
    do i = 1, M !!loop on frequencies
       y(i+4*(M+1)) = y(i+4*(M+1)) * derived_qtys % rho0
    end do
    !! dc portion
    y(5*(M+1)) = y(5*(M+1)) * derived_qtys % rho0    
  end subroutine unnorm_muse_vector

  subroutine phase_space(plot_cube, movie_number)
    complex(dp), dimension(:,:,:), intent(in) :: plot_cube
    integer, intent(in) :: movie_number

    character(19) :: file_name
    character(70) :: full_path
    logical :: part_plotted
    integer :: i, j, kk
    complex(dp), dimension(:), allocatable :: tmpPsi, tmpV
    real(dp) :: z, dz, omega0, betae, T, dt, dPsi, dv, fn, z_part, v_part


    !! tmpY takes a row of plot_cube rather than a column
    allocate(tmpPsi(output_parameters % num_axial_points))
    allocate(tmpV(output_parameters % num_axial_points))

    omega0 = 2.0*pi*frequency_parameters % base_frequency
    betae = omega0/derived_qtys % u0
    dz = circuit_parameters % circuit_length &
         / (output_parameters % num_axial_points - 1)

    dt = run_parameters % movie_data_array(movie_number) % time &
         / (run_parameters % movie_data_array(movie_number) % num_frames - 1)
    T = 0.0
    
    !! make a ps plot
    do i = 1, run_parameters % movie_data_array(movie_number) % num_frames
       !! create the file
       if (i < 10) then
          write(file_name,fmt="(A12, I1, A4)") 'phase_space.', &
            i, '.dat'
       else if (i < 100) then
          write(file_name,fmt="(A12, I2, A4)") 'phase_space.', &
            i, '.dat'
       else if (i < 1000) then
          write(file_name,fmt="(A12, I3, A4)") 'phase_space.', &
            i, '.dat'
       end if

       !! write the full file path
       if (i < 10) then
          write(full_path,fmt="(A16,I1,A1,A17)") &
               './outputs/movie_',movie_number,'/',file_name
       else if (i < 100) then
          write(full_path,fmt="(A16,I1,A1,A18)") &
               './outputs/movie_',movie_number,'/',file_name
       else if (i < 1000) then
          write(full_path,fmt="(A16,I1,A1,A19)") &
               './outputs/movie_',movie_number,'/',file_name
       endif
       
       !! open the file
       open(1,file=full_path,action='write')       

       !! loop on particles
       do j = 1, N_disks
          !! take the next vector
          tmpV = plot_cube(1,:,3*M+j)
          tmpPsi = plot_cube(1,:,3*M+N_disks+j)
          !! unnormalize the row
          do kk = 1, output_parameters % num_axial_points
             tmpV(kk) = tmpV(kk) * derived_qtys % u0
             tmpPsi(kk) = tmpPsi(kk)*w0
          end do

          !! compute fn at z=0. if < 0 then disk starts out below the
          !! \psi = \beta z - \omega t line and will never intersect it
          !! in this case add 2pi to its phase trajectory to include it.
          fn = real(tmpPsi(1)) + omega0 * T
          if (fn < 0.0) then
             do kk = 1, output_parameters % num_axial_points
                tmpPsi(kk) = tmpPsi(kk) + 2.0*pi
             end do
          end if

          !! loop on z
          z = dz
          part_plotted = .false.
          kk = 2
          do while (kk <= output_parameters % num_axial_points &
               .and. .not. part_plotted)
             fn = real(tmpPsi(kk)) - betae * z + omega0 * T
             !print*, 'j=',j,'kk=',kk, fn
             !pause
             !! when fn goes negative interpolate between these two points
             if (fn < 0.0) then
                dPsi = real(tmpPsi(kk)) - real(tmpPsi(kk-1))
                z_part = ((dPsi/dz)*(z-dz) - real(tmpPsi(kk-1)) - omega0*T) &
                     /(dPsi/dz - betae)

                dv = real(tmpV(kk)) - real(tmpV(kk-1))
                v_part = (dv/dz)*(z_part - (z - dz)) + real(tmpV(kk-1))

                if (units_structure % length_cm) then
                   z_part = z_part * 100.0
                   v_part = v_part * 100.0
                end if
                write(1,fmt="(en20.8, en20.8)") z_part, v_part
                part_plotted = .true.
             end if
             z = z + dz
             kk = kk + 1
          end do
       end do
       T = T + dt
    end do
    close(1)
  end subroutine phase_space


  subroutine create_file(i,j)
    integer, intent(in) :: i,j  !i = namelist #, j = frame #
    character(40) :: file_type
    character(10) :: dir_name
    character(70) :: file_name

    if (i <= 9) then
       write(dir_name,fmt="(A6,I1)")'movie_',i
    elseif (i >= 10 .and. i <= 99) then
       write(dir_name,fmt="(A6,I2)")'movie_',i
    endif

    !! setup file_type to match movie_type in movies.nml
    if (run_parameters % movie_data_array(i) % movie_type == 1) then   
       if (i <= 9) then
          write(file_type,fmt="(A10,A7,A20)")'./outputs/',&
                  dir_name,'/circ_voltage_frame.'
       elseif (i >= 10 .and. i <= 99) then
          write(file_type,fmt="(A10,A8,A20)")'./outputs/',&
                  dir_name,'/circ_voltage_frame.'
       endif
   
    elseif (run_parameters % movie_data_array(i) % movie_type == 2) then
       if (i <= 9) then
          write(file_type,fmt="(A10,A7,A20)")'./outputs/',&
                  dir_name,'/circ_current_frame.'
       elseif (i >= 10 .and. i <= 99) then
          write(file_type,fmt="(A10,A8,A20)")'./outputs/',&
                  dir_name,'/circ_current_frame.'
       endif

    elseif (run_parameters % movie_data_array(i) % movie_type == 3) then
       if (i <= 9) then
          write(file_type,fmt="(A10,A7,A20)")'./outputs/',&
                  dir_name,'/space_charge_frame.'
       elseif (i >= 10 .and. i <= 99) then
          write(file_type,fmt="(A10,A8,A20)")'./outputs/',&
                  dir_name,'/space_charge_frame.'
       endif

    elseif (run_parameters % movie_data_array(i) % movie_type == 4) then
       if (i <= 9) then
          write(file_type,fmt="(A10,A7,A20)")'./outputs/',&
                  dir_name,'/beam_velocty_frame.'
       elseif (i >= 10 .and. i <= 99) then
          write(file_type,fmt="(A10,A8,A20)")'./outputs/',&
                  dir_name,'/beam_velocty_frame.'
       endif

    elseif (run_parameters % movie_data_array(i) % movie_type == 5) then
       if (i <= 9) then
          write(file_type,fmt="(A10,A7,A20)")'./outputs/',&
                  dir_name,'/beam_density_frame.'
       elseif (i >= 10 .and. i <= 99) then
          write(file_type,fmt="(A10,A8,A20)")'./outputs/',&
                  dir_name,'/beam_density_frame.'
       endif

    elseif (run_parameters % movie_data_array(i) % movie_type == 6) then
       if (i <= 9) then
          write(file_type,fmt="(A10,A7,A20)")'./outputs/',&
                  dir_name,'/beam_current_frame.'
       elseif (i >= 10 .and. i <= 99) then
          write(file_type,fmt="(A10,A8,A20)")'./outputs/',&
                  dir_name,'/beam_current_frame.'
       endif
    endif
    
    ! setup frame file string
    if (i < 10) then
       if (j < 10) then
          write(file_name,fmt="(A37, I1, A4)") file_type,j,'.dat'     
       elseif (j > 9 .and. j < 100) then
          write(file_name,fmt="(A37, I2, A4)") file_type,j,'.dat'
       elseif (j > 99 .and. j < 1000) then
          write(file_name,fmt="(A37, I3, A4)") file_type,j,'.dat'
       endif
    elseif (i > 9 .and. i < 100) then
       if (j < 10) then
          write(file_name,fmt="(A38, I1, A4)") file_type,j,'.dat'     
       elseif (j > 9 .and. j < 100) then
          write(file_name,fmt="(A38, I2, A4)") file_type,j,'.dat'
       elseif (j > 99 .and. j < 1000) then
          write(file_name,fmt="(A38, I3, A4)") file_type,j,'.dat'
       endif
    endif

    !! open the file
    open(1,file=file_name,action='write')
  end subroutine create_file
end subroutine make_frames
