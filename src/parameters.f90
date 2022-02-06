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
!**************************************************************************!

module parameters
!**************************************************************!
!********************  modules to use  ************************!
!**************************************************************!
  use precision95
  use scan_data_class
  use movie_data_class

!** no implicit variable naming **!
  implicit none

  character(len=256) :: path_to_inputs_directory
  character(len=256) :: path_to_lmsuite_nml

!*********************************************************!
!********************  constants  ************************!
!*********************************************************!
  real(dp), parameter :: pi = 3.1415927
  real(dp), parameter :: me = 9.11e-31
  real(dp), parameter :: echarge = 1.6e-19
  real(dp), parameter :: eps0 = 8.85e-12
  real(dp), parameter :: mu0 = 1.25664e-6
  real(dp), parameter :: c = 2.99863e8     ! m/s
  complex(dp), parameter :: smo = (0.0,1.0)

!*** even in fortran 95, still need fix some array lengths :( *************!
!*** e.g. allocatable arrays cannot be used in derived type definitions ***!
  integer, parameter :: max_freqs = 500
  !integer, parameter :: max_freqs = 100
  !integer, parameter :: max_freqs = 20
  integer, parameter :: max_ckt_sections = 20
  integer, parameter :: max_scans = 10
  integer, parameter :: max_movies = 99  !! if this is 3 digits, must add code

!*** local variables ***!
  integer :: i, j


!*********************************************************!
!**** data structures with code, run, TWT, dispersion ****!
!***********  and frequency/power parameters  ************!
!*********************************************************!
  type interface_parameters_S
     logical :: echo_user_interface
     logical :: echo_namelists
     logical :: echo_initialization
     logical :: echo_runtime_one
     logical :: echo_runtime_two
     logical :: echo_runtime_three
     logical :: with_pauses
  end type interface_parameters_S

  type units_structure_S
     logical :: echo_units
     logical :: input_power_dBm
     logical :: output_power_dBm
     logical :: phase_degrees
     logical :: length_cm
     logical :: nepers
  end type units_structure_S

  type run_parameters_S
     logical :: echo_run
     character(3) :: select_code
     logical :: compute_reflections
     logical :: svea
     logical :: dc_constant
     !**** for scan input ****!
     integer :: num_scan_namelists
     type (scan_data), dimension(max_scans) :: scan_data_array
     integer :: num_movie_namelists
     type (movie_data), dimension(max_movies) :: movie_data_array
  end type run_parameters_S

  !!"qty" integers use the following legend
  !!1=circuit voltage, 2=circuit current, 3=space charge field,
  !!4 =beam velocity, 5=beam charge density, 6=beam current
  type output_parameters_S
     logical :: echo_output
     logical :: clean_outputs
     logical :: file_headers
     logical :: plot_dispersion
     logical :: power_out_vs_power_in
     logical :: gain_vs_power_in
     logical :: power_out_vs_phase_in
     logical :: phase_out_vs_power_in
     logical :: phase_out_vs_phase_in
     logical :: power_out_vs_freq
     logical :: gain_vs_freq
     logical :: phase_vs_freq
     integer :: num_axial_points
     logical :: circuit_power_vs_z
     logical :: magnitude_vs_z
     integer :: mvz_qty
     logical :: phase_vs_z
     integer :: pvz_qty
     logical :: phase_differences_vs_z
     integer :: pd_qty_1
     integer :: pd_qty_2
     logical :: hot_phase_velocity_vs_z
     logical :: disk_orbits_vs_z
     logical :: conserved_quantity_vs_z
     logical :: beam_energy_vs_z
     logical :: dc_beam_vel_vs_z
     logical :: phase_space
  end type output_parameters_S

  type frequency_parameters_S
     logical :: echo_frequencies
     integer :: num_input_freqs
     real(dp) :: base_frequency
     logical :: read_from_file
     integer, dimension(max_freqs) :: frequency_integer
     real(dp), dimension(max_freqs) :: power_input
     real(dp), dimension(max_freqs) :: phase_input
     integer :: highest_order_IMP
     integer :: min_ckt_freq
     integer :: max_ckt_freq
     integer :: min_space_charge_freq
     integer :: max_space_charge_freq
  end type frequency_parameters_S

  !! this is not complete. if we ever really wish to go to a dispersion
  !! solver, we need info about vane dimensions, dielectric rods, smeared, etc.
  type circuit_parameters_S
     logical :: echo_circuit
     real(dp) :: circuit_length
     integer :: number_ckt_sections
     logical :: read_from_file
     real(dp), dimension(max_ckt_sections) :: ckt_section_location
     ! first "section" is always at 0.0
     real(dp), dimension(max_ckt_sections) :: tape_width
     real(dp), dimension(max_ckt_sections) :: helix_radius
     real(dp), dimension(max_ckt_sections) :: helix_pitch
     real(dp), dimension(max_ckt_sections) :: vane_radius
     real(dp), dimension(max_ckt_sections) :: wall_radius
     real(dp) :: input_TL_impedance
     real(dp) :: output_TL_impedance
  end type circuit_parameters_S

  type dispersion_parameters_S
     logical :: echo_dispersion
     logical :: vph_over_c
     integer :: num_dispersion_freqs
     logical :: read_from_file
     integer, dimension(max_ckt_sections,max_freqs) :: dispersion_freq_integer
     real(dp), dimension(max_ckt_sections,max_freqs) :: phase_velocity
     real(dp), dimension(max_ckt_sections,max_freqs) :: impedance
     real(dp), dimension(max_ckt_sections,max_freqs) :: space_charge_redux
     logical :: intrplt_sections
     integer, dimension(max_ckt_sections - 1) :: interpol_sects_list
     logical :: intrplt_over_length
     real(dp) :: interpol_length
     logical :: use_sheath_model
     logical :: use_tape_model
     logical :: use_antonsen_formula
     logical :: use_correction
     logical :: use_beta_c
  end type dispersion_parameters_S

  type loss_parameters_S
     logical :: echo_losses
     integer :: number_loss_locations
     integer :: number_loss_freqs
     logical :: read_from_file
     real(dp), dimension(max_ckt_sections) :: loss_location
     integer, dimension(max_ckt_sections,max_freqs) :: loss_freq_integer
     real(dp), dimension(max_ckt_sections,max_freqs) :: loss
     logical :: intrplt_btwn_points
     logical :: intrplt_over_length
     real(dp) :: interpol_length
     logical :: use_loss_model
  end type loss_parameters_S

  type beam_parameters_S
     logical :: echo_beam
     logical :: relativistic
     logical :: DC_depression
     real(dp) :: voltage
     real(dp) :: current
     real(dp) :: outer_radius
     real(dp) :: inner_radius
     logical :: specify_velocity
     real(dp) :: beam_velocity
  end type beam_parameters_S

  type numerical_parameters_S
     logical :: echo_numerical
     integer :: num_grid_pts
     integer :: base_number_disks
     real(dp) :: fnt0
     logical :: adaptive_step
     real(dp) :: tolerance
  end type numerical_parameters_S

!*********************************************************!
!********  make variables of the types above  ************!
!*********************************************************!

  type (interface_parameters_S) interface_parameters
  type (units_structure_S) units_structure
  type (run_parameters_S) run_parameters
  type (output_parameters_S) output_parameters
  type (frequency_parameters_S) frequency_parameters
  type (circuit_parameters_S) circuit_parameters
  type (dispersion_parameters_S) dispersion_parameters
  type (loss_parameters_S) loss_parameters
  type (beam_parameters_S) beam_parameters
  type (numerical_parameters_S) numerical_parameters


contains
  !*********************************************************!
  !********  subroutine to set default parameters  *********!
  !*********************************************************!
  subroutine set_default_parameters
    !*** interface_parameters ***!
    interface_parameters % echo_user_interface = .true.
    interface_parameters % echo_namelists = .true.
    interface_parameters % echo_initialization = .true.
    interface_parameters % echo_runtime_one = .true.
    interface_parameters % echo_runtime_two = .true.
    interface_parameters % echo_runtime_three = .true.
    interface_parameters % with_pauses = .true.

    !*** units_structure ***!
    units_structure % echo_units = .true.
    units_structure % input_power_dBm = .false.
    units_structure % output_power_dBm = .true.
    units_structure % phase_degrees = .false.
    units_structure % length_cm = .false.
    units_structure % nepers = .true.

    !*** run_parameters ***!
    run_parameters % echo_run = .true.
    run_parameters % select_code = 'LMS'
    run_parameters % compute_reflections = .false.
    run_parameters % svea = .false.
    run_parameters % dc_constant = .false.
    run_parameters % num_scan_namelists = 0
    run_parameters % num_movie_namelists = 0

    !*** output_parameters ***!
    output_parameters % echo_output = .true.
    output_parameters % clean_outputs = .true.
    output_parameters % file_headers = .true.
    output_parameters % plot_dispersion = .true.
    output_parameters % power_out_vs_power_in = .false.
    output_parameters % gain_vs_power_in = .false.
    output_parameters % power_out_vs_phase_in = .false.
    output_parameters % phase_out_vs_power_in = .false.
    output_parameters % phase_out_vs_phase_in = .false.
    output_parameters % power_out_vs_freq = .true.
    output_parameters % gain_vs_freq = .false.
    output_parameters % phase_vs_freq = .true.
    output_parameters % num_axial_points = 100
    output_parameters % circuit_power_vs_z = .true.
    output_parameters % magnitude_vs_z = .true.
    output_parameters % mvz_qty = 5
    output_parameters % phase_vs_z = .true.
    output_parameters % pvz_qty = 1
    output_parameters % phase_differences_vs_z = .true.
    output_parameters % pd_qty_1 = 1
    output_parameters % pd_qty_2 = 5
    output_parameters % hot_phase_velocity_vs_z = .false.
    output_parameters % disk_orbits_vs_z = .false.
    output_parameters % conserved_quantity_vs_z = .false.
    output_parameters % beam_energy_vs_z = .false.
    output_parameters % dc_beam_vel_vs_z = .false.
    output_parameters % phase_space = .false.

    !*** frequency_parameters ***!
    frequency_parameters % echo_frequencies = .true.
    frequency_parameters % num_input_freqs = 1
    frequency_parameters % base_frequency = 1.0e6
    frequency_parameters % read_from_file = .false.
    frequency_parameters % frequency_integer = (/1600,(0, i=1, max_freqs-1)/)
    frequency_parameters % power_input = &
         (/1.0e-5, (0.0, i=1, max_freqs-1)/)
    frequency_parameters % phase_input = (/0.0, (0.0, i=1, max_freqs-1)/)
    frequency_parameters % highest_order_IMP = 2
    frequency_parameters % min_ckt_freq = 1600
    frequency_parameters % max_ckt_freq = 3200
    frequency_parameters % min_space_charge_freq = 0
    frequency_parameters % max_space_charge_freq = 3200
    
    !*** circuit_parameters ***!
    circuit_parameters % echo_circuit = .true.
    circuit_parameters % circuit_length = 0.42
    circuit_parameters % number_ckt_sections = 1
    circuit_parameters % read_from_file = .false.
    circuit_parameters % ckt_section_location = &
         (/ 0.0, (0.0, i=1, max_ckt_sections-1) /)
    ! first "section" is always at 0.0
    circuit_parameters % tape_width = &
         (/ 0.0, (0.0, i=1, max_ckt_sections-1) /)
    circuit_parameters % helix_radius = &
         (/ 0.2353e-2, (0.0, i=1, max_ckt_sections-1) /)
    circuit_parameters % helix_pitch = &
         (/ 0.0, (0.0, i=1, max_ckt_sections-1) /)
    circuit_parameters % vane_radius = &
         (/ 0.0, (0.0, i=1, max_ckt_sections-1) /)
    circuit_parameters % wall_radius = &
         (/ 0.0, (0.0, i=1, max_ckt_sections-1) /)
    circuit_parameters % input_TL_impedance = 50.0
    circuit_parameters % output_TL_impedance = 50.0
    
    !*** dispersion_parameters ***!
    dispersion_parameters % echo_dispersion = .true.
    dispersion_parameters % vph_over_c = .true.
    dispersion_parameters % num_dispersion_freqs = 2
    dispersion_parameters % read_from_file = .false.
    dispersion_parameters % dispersion_freq_integer&
         = reshape(source = (/1600, (0,i=1,max_ckt_sections-1), &
         3200, (0,i=1,max_ckt_sections-1), &
         (0,i=1,max_ckt_sections*(max_freqs-2))/), &
         shape = (/max_ckt_sections, max_freqs/))
    dispersion_parameters % phase_velocity = &
         reshape(source = (/0.0987921, (0.0,i=1,max_ckt_sections-1), &
         0.0879657, (0.0,i=1,max_ckt_sections-1), &
         (0.0,i=1,max_ckt_sections*(max_freqs-2))/), &
         shape = (/max_ckt_sections, max_freqs/))
    dispersion_parameters % impedance = &
         reshape(source = (/236.881, (0.0,i=1,max_ckt_sections-1), &
         39.7373, (0.0,i=1,max_ckt_sections-1), &
         (0.0,i=1,max_ckt_sections*(max_freqs-2))/), &
         shape = (/max_ckt_sections, max_freqs/))
    dispersion_parameters % space_charge_redux = &
         reshape(source = (/0.0, (0.0,i=1,max_ckt_sections-1), &
         0.0, (0.0,i=1,max_ckt_sections-1), &
         (0.0,i=1,max_ckt_sections*(max_freqs-2))/), &
         shape = (/max_ckt_sections, max_freqs/))
    dispersion_parameters % intrplt_over_length = .false.
    dispersion_parameters % interpol_sects_list = &
         (/ (-1, i=1,max_ckt_sections-1)/)
    dispersion_parameters % intrplt_over_length = .true.
    dispersion_parameters % interpol_length = 0.01
    dispersion_parameters % use_sheath_model = .false.
    dispersion_parameters % use_tape_model = .false.
    dispersion_parameters % use_antonsen_formula = .true.
    dispersion_parameters % use_correction = .false.
    dispersion_parameters % use_beta_c = .false.

    do i = 1, max_ckt_sections
       do j = 1, max_freqs
          dispersion_parameters % phase_velocity(i,j) = &
            c * dispersion_parameters % phase_velocity(i,j)
       end do
    end do

    !*** loss_parameters ***!
    loss_parameters % echo_losses = .true.
    loss_parameters % number_loss_locations = 1
    loss_parameters % number_loss_freqs = 2
    loss_parameters % read_from_file = .false.
    ! first "loss_location" is always at 0.0
    loss_parameters % loss_location = &
         (/ 0.0, (0.0, i=1, max_ckt_sections-1) /)
    loss_parameters % loss_freq_integer = &
         reshape(source = (/1600, (0,i=1,max_ckt_sections-1), &
         3200, (0,i=1,max_ckt_sections-1), &
         (0,i=1,max_ckt_sections*(max_freqs-2))/), &
         shape = (/max_ckt_sections, max_freqs/))
    loss_parameters % loss = &
         reshape(source = (/0.1, (0.0,i=1,max_ckt_sections-1), &
         0.1, (0.0,i=1,max_ckt_sections-1), &
         (0.0,i=1,max_ckt_sections*(max_freqs-2))/), &
         shape = (/max_ckt_sections, max_freqs/))
    loss_parameters % intrplt_btwn_points = .false.
    loss_parameters % intrplt_over_length = .true.
    loss_parameters % interpol_length = 0.01
    loss_parameters % use_loss_model = .false.

    !*** beam_parameters ***!
    beam_parameters % echo_beam = .true.
    beam_parameters % relativistic = .false.
    beam_parameters % DC_depression = .false.
    beam_parameters % voltage = 3100
    beam_parameters % current = 0.0655
    beam_parameters % outer_radius = 0.09652e-2   ! in meters
    beam_parameters % inner_radius = 0.0   ! in meters
    beam_parameters % specify_velocity = .false.
    beam_parameters % beam_velocity = 1.0e9

    !*** numerical_parameters ***!
    numerical_parameters % echo_numerical = .true.
    numerical_parameters % num_grid_pts = 500
    numerical_parameters % base_number_disks = 19 !! LATTE only
    numerical_parameters % fnt0 = 1.5 !!copies christine parameter, LATTE only
    numerical_parameters % adaptive_step = .false.
    numerical_parameters % tolerance = 1.0e-6
    
  end subroutine set_default_parameters

  !*********************************************************!
  !******  subroutine to read parameters from file  ********!
  !*********************************************************!
  subroutine read_parameters

    !*********************************!
    !**** namelist user_interface ****!
    !*********************************!
    logical echo_user_interface, &
         echo_namelists, &
         echo_initialization, &
         echo_runtime_one, &
         echo_runtime_two, &
         echo_runtime_three, &
         with_pauses
    namelist /user_interface/ echo_user_interface, &
         echo_namelists, &
         echo_initialization, &
         echo_runtime_one, &
         echo_runtime_two, &
         echo_runtime_three, &
         with_pauses
    !************************!
    !**** namelist units ****!
    !************************!
    logical echo_units, &
         input_power_dBm, &
         output_power_dBm, &
         phase_degrees, &
         length_cm, &
         nepers
    namelist /units/ echo_units, &
         input_power_dBm, &
         output_power_dBm, &
         phase_degrees, &
         length_cm, &
         nepers
    !**********************!
    !**** namelist run ****!
    !**********************!
    character(3) select_code
    logical echo_run, &
         compute_reflections, &
         svea, &
         dc_constant
    integer num_scan_namelists, &
         num_movie_namelists
    namelist /run/ echo_run, &
         select_code, &
         compute_reflections, &
         svea, &
         dc_constant, &
         num_scan_namelists, &
         num_movie_namelists
    !***********************!
    !**** namelist scan ****!
    !***********************!
    logical two_parameter, distribute_log, logic_param_1
    integer scanID, &
         num_points, &
         int_param_1, &
         int_param_2, &
         int_param_3
    real(dp) min, max, real_param_1, real_param_2
    namelist /scan/ scanID, &
         two_parameter, &
         min, &
         max, &
         num_points, &
         distribute_log, &
         int_param_1, &
         int_param_2, &
         int_param_3, &
         real_param_1, &
         real_param_2, &
         logic_param_1
    !************************!
    !**** namelist movie ****!
    !************************!
    integer :: movie_type, num_frames
    real(dp) :: time, scale
    namelist /movie/ movie_type, &
         time, &
         num_frames, &
         scale
    !******************************!
    !**** namelist output_data ****!
    !******************************!
    logical :: echo_output, &
         clean_outputs, &
         file_headers, &
         plot_dispersion, &
         power_out_vs_power_in, &
         gain_vs_power_in, &
         power_out_vs_phase_in, &
         phase_out_vs_power_in, &
         phase_out_vs_phase_in, &
         power_out_vs_freq, &
         gain_vs_freq, &
         phase_vs_freq, &
         circuit_power_vs_z, magnitude_vs_z, phase_vs_z, &
         phase_differences_vs_z, hot_phase_velocity_vs_z, disk_orbits_vs_z, &
         conserved_quantity_vs_z, beam_energy_vs_z, dc_beam_vel_vs_z, &
         phase_space
    integer :: num_axial_points, mvz_qty, pvz_qty, pd_qty_1, pd_qty_2
    namelist /output_data/ echo_output, &
         clean_outputs, &
         file_headers, &
         plot_dispersion, &
         power_out_vs_power_in, &
         gain_vs_power_in, &
         power_out_vs_phase_in, &
         phase_out_vs_power_in, &
         phase_out_vs_phase_in, &
         power_out_vs_freq, &
         gain_vs_freq, &
         phase_vs_freq, &
         num_axial_points, &
         circuit_power_vs_z, &
         magnitude_vs_z, &
         mvz_qty, &
         phase_vs_z, &
         pvz_qty, &
         phase_differences_vs_z, &
         pd_qty_1, &
         pd_qty_2, &
         hot_phase_velocity_vs_z, &
         disk_orbits_vs_z, &
         conserved_quantity_vs_z, &
         beam_energy_vs_z, &
         dc_beam_vel_vs_z, &
         phase_space
    !*********************************!
    !**** namelist frequency_list ****!
    !*********************************!
    logical echo_frequencies, &
         read_from_file
    integer num_input_freqs, &
         frequency_integer(max_freqs), &
         highest_order_IMP, &
         min_ckt_freq, &
         max_ckt_freq, &
         min_space_charge_freq, &
         max_space_charge_freq
    real(dp) base_frequency, &
         power_input(max_freqs), &
         phase_input(max_freqs)
    namelist /frequency_list/ echo_frequencies, &
         num_input_freqs, &
         base_frequency, &
         read_from_file, &
         frequency_integer, &
         power_input, &
         phase_input, &
         highest_order_IMP, &
         min_ckt_freq, &
         max_ckt_freq, &
         min_space_charge_freq, &
         max_space_charge_freq
    !**************************!
    !**** namelist circuit ****!
    !**************************!
    logical echo_circuit
         !read_from_file, &
    integer number_ckt_sections
    real(dp) circuit_length, &
         ckt_section_location(max_ckt_sections), &
         tape_width(max_ckt_sections), &
         helix_radius(max_ckt_sections), &
         helix_pitch(max_ckt_sections), &
         vane_radius(max_ckt_sections), &
         wall_radius(max_ckt_sections), &
         input_TL_impedance, &
         output_TL_impedance
    namelist /circuit/ echo_circuit, &
         circuit_length, &
         number_ckt_sections, &
         read_from_file, &
         ckt_section_location, &
         tape_width, &
         helix_radius, &
         helix_pitch, &
         vane_radius, &
         wall_radius, &
         input_TL_impedance, &
         output_TL_impedance
    !*****************************!
    !**** namelist dispersion ****!
    !*****************************!
    logical echo_dispersion, &
         vph_over_c, &
         !read_from_file, &
         intrplt_sections, &
         intrplt_over_length, &
         use_sheath_model, &
         use_tape_model, &
         use_antonsen_formula, &
         use_correction, &
         use_beta_c
    integer num_dispersion_freqs, &
         dispersion_freq_integer(max_ckt_sections * max_freqs), &
         interpol_sects_list(max_ckt_sections - 1)
    !these arrays are 1-dimensional, will make them 2-d upon reading
    real(dp) phase_velocity(max_ckt_sections * max_freqs), &
         impedance(max_ckt_sections * max_freqs), &
         space_charge_redux(max_ckt_sections * max_freqs), &
         interpol_length
    namelist /dispersion/ echo_dispersion, &
         vph_over_c, &
         num_dispersion_freqs, &
         read_from_file, &
         dispersion_freq_integer, &
         phase_velocity, &
         impedance, &
         space_charge_redux, &
         intrplt_sections, &
         interpol_sects_list, &
         intrplt_over_length, &
         interpol_length, &
         use_sheath_model, &
         use_tape_model, &
         use_antonsen_formula, &
         use_correction, &
         use_beta_c
    !*************************!
    !**** namelist losses ****!
    !*************************!
    logical echo_losses, &
         !read_from_file, &
         intrplt_btwn_points, &
         !intrplt_over_length, &
         use_loss_model
    integer number_loss_locations, &
         number_loss_freqs, &
         loss_freq_integer(max_ckt_sections * max_freqs)
    real(dp) loss_location(max_ckt_sections), &
         loss(max_ckt_sections * max_freqs)
         !interpol_length
    namelist /losses/ echo_losses, &
         number_loss_locations, &
         number_loss_freqs, &
         read_from_file, &
         loss_location, &
         loss_freq_integer, &
         loss, &
         intrplt_btwn_points, &
         intrplt_over_length, &
         interpol_length, &
         use_loss_model
    !***********************!
    !**** namelist beam ****!
    !***********************!
    logical echo_beam, relativistic, DC_depression, specify_velocity
    real(dp) voltage, current, outer_radius, inner_radius, beam_velocity
    namelist /beam/ echo_beam, relativistic, &
         DC_depression, voltage, current, outer_radius, &
         inner_radius, specify_velocity, beam_velocity
    !*********************************!
    !**** namelist numerical_data ****!
    !*********************************!
    logical echo_numerical, adaptive_step
    integer num_grid_pts, base_number_disks
    real(dp) fnt0, tolerance
    namelist /numerical_data/ echo_numerical, &
         num_grid_pts, base_number_disks, fnt0, adaptive_step, &
         tolerance

    !**** local variables ****!
    logical :: file_stat !variable to check existence of input files
    real(dp) tmp

!***************  start reading parameters  ****************!
10  format(/,/,/,'Reading parameters from lmsuite.nml',/)
    print 10

    open(1,file=path_to_lmsuite_nml)

!*****************  read namelist user_interface  *****************!
    print*, 'Reading namelist user_interface'
    read(1,user_interface)
    interface_parameters % echo_user_interface = echo_user_interface
    interface_parameters % echo_namelists = echo_namelists
    interface_parameters % echo_initialization = echo_initialization
    interface_parameters % echo_runtime_one = echo_runtime_one
    interface_parameters % echo_runtime_two = echo_runtime_two
    interface_parameters % echo_runtime_three = echo_runtime_three
    interface_parameters % with_pauses = with_pauses

!*****************  read namelist units  *****************!
    print*, 'Reading namelist units'
    read(1,units)
    units_structure % echo_units = echo_units
    units_structure % input_power_dBm = input_power_dBm
    units_structure % output_power_dBm = output_power_dBm
    units_structure % phase_degrees = phase_degrees
    units_structure % length_cm = length_cm
    units_structure % nepers = nepers

!*****************  read namelist run  *****************!
    print*, 'Reading namelist run'
    read(1,run)
    run_parameters % echo_run = echo_run
    run_parameters % select_code = select_code
    run_parameters % compute_reflections = compute_reflections
    run_parameters % svea = svea
    run_parameters % dc_constant = dc_constant
    run_parameters % num_scan_namelists = num_scan_namelists

    !**** read namelist scan if there are scans ****!
    if (run_parameters % num_scan_namelists > 0) then
       print "(t5,'reading scan_data.nml')"
       inquire(file='./inputs/scan_data.nml', exist=file_stat)
       if (.not. file_stat) then
          print*, ''
          print*, 'File inputs/scan_data.nml does not exist. Stopping.'
          print*, ''
          stop
       else if (file_stat) then
          open(2,file='./inputs/scan_data.nml')
          do i = 1, run_parameters % num_scan_namelists
             read(2,scan)
             call init_scan_data(run_parameters % scan_data_array(i), &
                  scanID, two_parameter, min, max, num_points, &
                  distribute_log, int_param_1, int_param_2, int_param_3, &
                  real_param_1, real_param_2, logic_param_1)
          end do
          close(2)
       end if
    end if

    run_parameters % num_movie_namelists = num_movie_namelists

    !**** read namelist movies if there are movies ****!
    if (run_parameters % num_movie_namelists > 0) then
       print "(t5,'reading movies.nml')"
       inquire(file=trim(path_to_inputs_directory)//"movies.nml", &
            exist=file_stat)
       if (.not. file_stat) then
          print*, ''
          print*, 'File inputs/movies.nml does not exist. Stopping.'
          print*, ''
          stop
       else if (file_stat) then
          open(2,file=trim(path_to_inputs_directory)//"movies.nml")
          do i = 1, run_parameters % num_movie_namelists
             read(2,movie)
             call init_movie_data(run_parameters % movie_data_array(i), &
                  movie_type, time, num_frames, scale)
          end do
          close(2)
       end if
    end if

!*****************  read namelist output_data  *****************!
    print*, 'Reading namelist output_data'
    read(1,output_data) 
    output_parameters % echo_output = echo_output
    output_parameters % clean_outputs = clean_outputs
    output_parameters % file_headers = file_headers
    output_parameters % plot_dispersion = plot_dispersion
    output_parameters % power_out_vs_power_in = power_out_vs_power_in
    output_parameters % gain_vs_power_in = gain_vs_power_in
    output_parameters % power_out_vs_phase_in = power_out_vs_phase_in
    output_parameters % phase_out_vs_power_in = phase_out_vs_power_in
    output_parameters % phase_out_vs_phase_in = phase_out_vs_phase_in
    output_parameters % power_out_vs_freq = power_out_vs_freq
    output_parameters % gain_vs_freq = gain_vs_freq
    output_parameters % phase_vs_freq = phase_vs_freq
    output_parameters % num_axial_points = num_axial_points
    output_parameters % circuit_power_vs_z = circuit_power_vs_z
    output_parameters % magnitude_vs_z = magnitude_vs_z
    output_parameters % mvz_qty = mvz_qty
    if (mvz_qty < 1 .or. mvz_qty > 6) then
       print*, 'Magnitude vs z quantity (mvz_qty in &output_data) is not between 1 and 6.'
       print*, 'Stopping.'
       stop
    end if
    output_parameters % phase_vs_z = phase_vs_z
    output_parameters % pvz_qty = pvz_qty
    if (pvz_qty < 1 .or. pvz_qty > 6) then
       print*, 'Phase vs z quantity (pvz_qty in &output_data) is not between 1 and 6.'
       print*, 'Stopping.'
       stop
    end if
    output_parameters % phase_differences_vs_z = phase_differences_vs_z
    output_parameters % pd_qty_1 = pd_qty_1
    output_parameters % pd_qty_2 = pd_qty_2
    if (pd_qty_1 < 1 .or. pd_qty_1 > 6 .or. &
         pd_qty_2 < 1 .or. pd_qty_2 > 6) then
       print*, 'Phase difference vs z quantity (pd_qty_1 or pd_qty_2 in &output_data)'
       print*, 'is not between 1 and 6.'
       print*, 'Stopping.'
       stop
    end if
    output_parameters % hot_phase_velocity_vs_z = hot_phase_velocity_vs_z
    output_parameters % disk_orbits_vs_z = disk_orbits_vs_z
    output_parameters % conserved_quantity_vs_z = conserved_quantity_vs_z
    output_parameters % beam_energy_vs_z = beam_energy_vs_z
    output_parameters % dc_beam_vel_vs_z = dc_beam_vel_vs_z
    output_parameters % phase_space = phase_space

!*****************  read namelist frequency_list ******************!
    print*, 'Reading namelist frequency_list'
    read(1,frequency_list)
    frequency_parameters % echo_frequencies = echo_frequencies
    frequency_parameters % num_input_freqs = num_input_freqs
    if (num_input_freqs > max_freqs) then
       print*, ''
       print*, 'num_input_freqs is larger than max_freqs.'
       print*, 'To run this many input frequencies, change max_freqs in'
       print*, 'source file parameters.f90 and recompile.'
       print*, 'Stopping.'
       stop
    end if
    frequency_parameters % base_frequency = base_frequency
    frequency_parameters % read_from_file = read_from_file
    !if read_from_file, overwrite data read from namelist
    if (frequency_parameters % read_from_file) then
       print "(t5,'reading file frequencies.in')"
       inquire(file='./inputs/frequencies.in', exist=file_stat)
       if (.not.file_stat) then
          print*, ''
          print*, 'File inputs/frequencies.in does not exist. Stopping.'
          print*, ''
          stop
       else if (file_stat) then
          open(2,file='./inputs/frequencies.in')
          do i = 1, frequency_parameters % num_input_freqs
             read(2, *) frequency_integer(i), &
                  power_input(i), &
                  phase_input(i)
          end do
          close(2)
       end if
    end if
    frequency_parameters % frequency_integer = frequency_integer
    !assign powers and phases, converting if necessary
    do i = 1, frequency_parameters % num_input_freqs
       if (units_structure % input_power_dBm) then
          frequency_parameters % power_input(i) = &
               0.001 * (10.0 ** (power_input(i)/10.0))
       else
          frequency_parameters % power_input(i) = power_input(i)
       end if
       if (units_structure % phase_degrees) then
          frequency_parameters % phase_input(i) = &
               (pi/180.0)*phase_input(i)
       else
          frequency_parameters % phase_input(i) = phase_input(i)
       end if
    end do
    frequency_parameters % highest_order_IMP = highest_order_IMP
    frequency_parameters % min_ckt_freq = min_ckt_freq
    frequency_parameters % max_ckt_freq = max_ckt_freq
    frequency_parameters % min_space_charge_freq = min_space_charge_freq
    frequency_parameters % max_space_charge_freq = max_space_charge_freq

!*****************  read namelist circuit ******************!
    print*, 'Reading namelist circuit'
    read(1,circuit)
    circuit_parameters % echo_circuit = echo_circuit
    circuit_parameters % number_ckt_sections = number_ckt_sections
    circuit_parameters % read_from_file = read_from_file
    !if read_from_file, overwrite data read into namelist
    if (circuit_parameters % read_from_file) then
       print "(t5,'reading file circuit.in')"
       inquire(file=trim(path_to_inputs_directory)//"/circuit.in", &
            exist=file_stat)
       if (.not. file_stat) then
          print*, ''
          print*, 'File inputs/circuit.in does not exist. Stopping.'
          print*, ''
          stop
       else if (file_stat) then
          open(2,file=trim(path_to_inputs_directory)//"/circuit.in")
          do i = 1, circuit_parameters % number_ckt_sections
             read(2, *) ckt_section_location(i), &
                  tape_width(i), &
                  helix_radius(i), &
                  helix_pitch(i), &
                  vane_radius(i), &
                  wall_radius(i)
          end do
          close(2)
       end if
    end if
    !assign circuit section information, converting to m if necessary
    do i = 1, circuit_parameters % number_ckt_sections
       if (units_structure % length_cm) then
          circuit_parameters % circuit_length = circuit_length/100.0
          circuit_parameters % ckt_section_location(i) = &
               ckt_section_location(i)/100.0
          circuit_parameters % tape_width(i) = tape_width(i)/100.0
          circuit_parameters % helix_radius(i) = helix_radius(i)/100.0
          circuit_parameters % helix_pitch(i) = helix_pitch(i)/100.0
          circuit_parameters % vane_radius(i) = vane_radius(i)/100.0
          circuit_parameters % wall_radius(i) = wall_radius(i)/100.0
       else
          circuit_parameters % circuit_length = circuit_length
          circuit_parameters % ckt_section_location(i) = &
               ckt_section_location(i)
          circuit_parameters % tape_width(i) = tape_width(i)
          circuit_parameters % helix_radius(i) = helix_radius(i)
          circuit_parameters % helix_pitch(i) = helix_pitch(i)
          circuit_parameters % vane_radius(i) = vane_radius(i)
          circuit_parameters % wall_radius(i) = wall_radius(i)
       end if
    end do
    circuit_parameters % input_TL_impedance = input_TL_impedance
    circuit_parameters % output_TL_impedance = output_TL_impedance

    !*****************  read namelist dispersion ******************!
    print*, 'Reading namelist dispersion'

    read(1,dispersion)
    dispersion_parameters % echo_dispersion = echo_dispersion
    dispersion_parameters % vph_over_c = vph_over_c
    dispersion_parameters % num_dispersion_freqs = num_dispersion_freqs
    dispersion_parameters % read_from_file = read_from_file

    if (dispersion_parameters % read_from_file) then
       print "(t5,'reading file dispersion.in')"
       inquire(file=trim(path_to_inputs_directory)//'/dispersion.in', &
            exist=file_stat)
       if (.not. file_stat) then
          print*, ''
          print*, 'File inputs/dispersion.in does not exist. Stopping.'
          print*, ''
          stop
       else if (file_stat) then
          ! if read_from_file, then read data from dispersion.in directly into
          ! dispersion_parameters % dispersion_freq_integer
          ! dispersion_parameters % phase_velocity
          ! dispersion_parameters % impedance
          ! dispersion_parameters % space_charge_redux
          ! if reads into namelist goes into array. then this array is put into
          ! matrices. this affects conversions. in the read_from_file case
          ! when a conversion is done one is overwriting data. in the namelist
          ! case one doesn't overwrite the array data, so parameters that
          ! need multiple coversions have to be processed differently.
          open(2,file=trim(path_to_inputs_directory)//'/dispersion.in')
          do i = 1, circuit_parameters % number_ckt_sections
             do j = 1, dispersion_parameters % num_dispersion_freqs
                read(2,*) dispersion_parameters % &
                     dispersion_freq_integer(i,j), &
                     dispersion_parameters % phase_velocity(i,j), &
                     dispersion_parameters % impedance(i,j), &
                     dispersion_parameters % space_charge_redux(i,j)
                !! convert phase velocity if necessary
                if (dispersion_parameters % vph_over_c) then
                   dispersion_parameters % phase_velocity(i,j) = &
                        dispersion_parameters % phase_velocity(i,j) * c
                else if (.not. dispersion_parameters % vph_over_c .and. &
                     units_structure % length_cm) then
                   dispersion_parameters % phase_velocity(i,j) = &
                        dispersion_parameters % phase_velocity(i,j) / 100.0
                end if
             end do
          end do
          close(2)
       end if
    else !not reading from dispersion.in
       do i = 1, circuit_parameters % number_ckt_sections
          do j = 1, dispersion_parameters % num_dispersion_freqs
             dispersion_parameters % dispersion_freq_integer(i,j) = &
                  dispersion_freq_integer(dispersion_parameters % &
                  num_dispersion_freqs * (i-1)+j)
             !! this allows multiple assignments, ending with the highest
             !! priority.
             !assign as if m/s
             dispersion_parameters % phase_velocity(i,j) = &
                  phase_velocity(dispersion_parameters % &
                  num_dispersion_freqs  * (i-1)+j)
             ! if cm_per_s convert cm/s to m/s.
             if (units_structure % length_cm) then
                dispersion_parameters % phase_velocity(i,j) = &
                     phase_velocity(dispersion_parameters % &
                     num_dispersion_freqs * (i-1)+j) / 100.0
             end if
             ! highest priority: if vph_over_c convert to m/s
             if (dispersion_parameters % vph_over_c) then
                dispersion_parameters % phase_velocity(i,j) = &
                     phase_velocity(dispersion_parameters % &
                     num_dispersion_freqs * (i-1)+j) * c
             end if
             dispersion_parameters % impedance(i,j) = &
                  impedance(dispersion_parameters % num_dispersion_freqs &
                  * (i-1)+j)
             dispersion_parameters % space_charge_redux(i,j) = &
                  space_charge_redux(dispersion_parameters % &
                  num_dispersion_freqs * (i-1)+j)
          end do
       end do
    end if
    dispersion_parameters % intrplt_sections = &
         intrplt_sections

    do i = 1, circuit_parameters % number_ckt_sections - 1
       if (interpol_sects_list(i) /= 0) then
          dispersion_parameters % interpol_sects_list(i) = &
               interpol_sects_list(i)
       else
          dispersion_parameters % interpol_sects_list(i) = -1
       end if
    end do

    dispersion_parameters % intrplt_over_length = &
         intrplt_over_length
    if (units_structure % length_cm) then
       dispersion_parameters % interpol_length = interpol_length/100.0
    else
       dispersion_parameters % interpol_length = interpol_length
    end if
    dispersion_parameters % use_sheath_model = use_sheath_model
    dispersion_parameters % use_tape_model = use_tape_model
    dispersion_parameters % use_antonsen_formula = use_antonsen_formula
    dispersion_parameters % use_correction = use_correction
    dispersion_parameters % use_beta_c = use_beta_c

!*****************  read namelist losses ******************!
    print*, 'Reading namelist losses'
    read(1,losses)
    loss_parameters % echo_losses = echo_losses
    loss_parameters % number_loss_locations = number_loss_locations
    loss_parameters % number_loss_freqs = number_loss_freqs
    loss_parameters % read_from_file = read_from_file

    if (loss_parameters % read_from_file) then
       print "(t5,'reading file losses.in')"
       inquire(file='./inputs/losses.in', exist=file_stat)
       if (.not. file_stat) then
          print*, ''
          print*, 'File inputs/losses.in does not exist. Stopping.'
          print*, ''
          stop
       else if (file_stat) then
          ! if read_from_file, then read data from dispersion.in directly into
          ! loss_parameters % loss_location
          ! loss_parameters % loss_freq_integer
          ! loss_parameters % loss
          ! if reads into namelist goes into array. then this array is put into
          ! matrices. this affects conversions. in the read_from_file case
          ! when a conversion is done one is overwriting data. in the namelist
          ! case one doesn't overwrite the array data, so parameters that
          ! need multiple coversions have to be processed differently.
          open(2,file='./inputs/losses.in')
          ! the first row are the loss locations
          read(2,*) (loss_parameters % loss_location(i), &
               i=1, loss_parameters % number_loss_locations)
          !convert loss locations to m if necessary
          do i = 1, loss_parameters % number_loss_locations
             if (units_structure % length_cm) then
                loss_parameters % loss_location(i) = &
                     loss_parameters % loss_location(i)/100.0
             end if
          end do
          do i = 1, loss_parameters % number_loss_locations
             do j = 1, loss_parameters % number_loss_freqs
                read(2,*) loss_parameters % loss_freq_integer(i,j), &
                     loss_parameters % loss(i,j)
                !! convert loss if necessary. first convert dB to nepers, then
                !! convert length
                if (.not. units_structure % nepers) then
                   loss_parameters % loss(i,j) = &
                        loss_parameters % loss(i,j) / 8.68589
                end if
                if (units_structure % length_cm) then
                   loss_parameters % loss(i,j) = &
                        loss_parameters % loss(i,j) * 100.0
                end if
             end do
          end do
          close(2)
       end if
    else !not reading from losses.in
       do i = 1, loss_parameters % number_loss_locations
             !convert loss locations to m if necessary
             if (units_structure % length_cm) then
                loss_parameters % loss_location(i) = loss_location(i)/100.0
             else
                loss_parameters % loss_location(i) = loss_location(i)
             end if
             do j = 1, loss_parameters % number_loss_freqs
                loss_parameters % loss_freq_integer(i,j) = &
                     loss_freq_integer(loss_parameters % &
                     number_loss_freqs * (i-1)+j)
                !! convert loss if necessary. this coversion has three cases
                !! since reading from an array rather than overwriting matrix.
                ! given dB/cm
                if (.not. units_structure % nepers .and. &
                     units_structure % length_cm) then
                   loss_parameters % loss(i,j) = &
                        loss(loss_parameters % number_loss_freqs &
                        * (i-1)+j) / 8.68589 * 100.0
                   ! given dB/m
                else if (.not. units_structure % nepers .and. &
                     .not. units_structure % length_cm) then
                   loss_parameters % loss(i,j) = &
                        loss(loss_parameters % number_loss_freqs &
                        * (i-1)+j) / 8.68589
                   ! given Np/cm
                else if (units_structure % nepers .and. &
                     units_structure % length_cm) then
                   loss_parameters % loss(i,j) = &
                        loss(loss_parameters % number_loss_freqs &
                        * (i-1)+j) * 100.0
                end if
             end do
       end do
    end if

    loss_parameters % intrplt_btwn_points = &
         intrplt_btwn_points
    loss_parameters % intrplt_over_length = &
         intrplt_over_length
    if (units_structure % length_cm) then
       loss_parameters % interpol_length = interpol_length/100.0
    else
       loss_parameters % interpol_length = interpol_length
    end if
    loss_parameters % use_loss_model = use_loss_model

!*****************  read namelist beam ******************!
    print*, 'Reading namelist beam'
    read(1,beam)
    beam_parameters % echo_beam = echo_beam
    beam_parameters % relativistic = relativistic
    beam_parameters % DC_depression = DC_depression
    beam_parameters % voltage = voltage
    beam_parameters % current = current
    if (units_structure % length_cm) then
       beam_parameters % outer_radius = outer_radius/100.0
       beam_parameters % inner_radius = inner_radius/100.0
       beam_parameters % beam_velocity = beam_velocity/100.0
    else
       beam_parameters % outer_radius = outer_radius
       beam_parameters % inner_radius = inner_radius
       beam_parameters % beam_velocity = beam_velocity
    end if
    beam_parameters % specify_velocity = specify_velocity
    
!*****************  read namelist numerical_data ******************!
    print*, 'Reading namelist numerical_data'
    read(1,numerical_data)
    numerical_parameters % echo_numerical = echo_numerical
    numerical_parameters % num_grid_pts = num_grid_pts
    if (num_grid_pts < num_axial_points) then
       print*, ''
       print*, 'num_grid_pts is less than num_axial_pts.'
       print*, 'Stopping.'
       stop
    end if
    numerical_parameters % base_number_disks = base_number_disks
    numerical_parameters % fnt0 = fnt0
    numerical_parameters % adaptive_step = adaptive_step
    numerical_parameters % tolerance = tolerance
    close(1)
  end subroutine read_parameters

!*********************************************************!
!******  subroutine to echo parameters to screen  ********!
!*********************************************************!
  subroutine echo_parameters
    character junk
    integer i, j
    real(dp) tmpPower, tmpPhase, tmpCSL, tmpTW, tmpHR, tmpHP, &
         tmpVR, tmpWR, tmpVPH, tmpLOC, tmpAlpha

20  format (/)
30  format(A)
    print 20

!*****************  echo namelist user_interface ******************!
    if (interface_parameters % echo_user_interface) then
       print*, 'namelist user_interface:'
       print*, '------------------------'
       print "(' echo_user_interface = ', L1)", &
            interface_parameters % echo_user_interface
       print "(' echo_namelists = ', L1)", &
            interface_parameters % echo_namelists
       print "(' echo_initialization = ', L1)", &
            interface_parameters % echo_initialization
       print "(' echo_runtime_one = ', L1)", &
            interface_parameters % echo_runtime_one
       print "(' echo_runtime_two = ', L1)", &
            interface_parameters % echo_runtime_two
       print "(' echo_runtime_three = ', L1)", &
            interface_parameters % echo_runtime_three
       print "(' with_pauses = ', L1)", interface_parameters % with_pauses
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist units  ******************!
    if (units_structure % echo_units) then
       print*, 'namelist units:'
       print*, '---------------'
       print "(' echo_units = ', L1)", units_structure % echo_units
       print "(' input_power_dBm = ', L1)", units_structure % input_power_dBm
       print "(' output_power_dBm = ', L1)", units_structure % output_power_dBm
       print "(' phase_degrees = ', L1)", units_structure % phase_degrees
       print "(' length_cm = ', L1)", units_structure % length_cm
       print "(' nepers = ', L1)", units_structure % nepers
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist run  ******************!
    if (run_parameters % echo_run) then
       print*, 'namelist run:'
       print*, '-------------'
       print "(' echo_run = ', L1)", run_parameters % echo_run
       print "(' select_code = ', A3)", run_parameters % select_code
       print "(' compute_reflections = ', L1)", &
            run_parameters % compute_reflections
       print "(' svea = ', L1)", run_parameters % svea
       print "(' dc_constant = ', L1)", run_parameters % dc_constant
       print "(' num_scan_namelists = ', I3)", &
            run_parameters % num_scan_namelists
       print "(' num_movie_namelists = ', I3)", &
            run_parameters % num_movie_namelists

       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if

       do i = 1, run_parameters % num_scan_namelists
          print "(' scan# = ', I2)", i
          print "(' ----------')"
          print "(' scanID = ', I2)", &
               run_parameters % scan_data_array(i) % scanID
          print "(' two_parameter = ', L1)", &
               run_parameters % scan_data_array(i) % two_parameter
          print "(' min = ', en11.2)", &
               run_parameters % scan_data_array(i) % min
          print "(' max = ', en11.2)", &
               run_parameters % scan_data_array(i) % max
          print "(' num_points = ', I4)", &
               run_parameters % scan_data_array(i) % num_points
          print "(' distribute_log = ', L1)", &
               run_parameters % scan_data_array(i) % distribute_log
          print "(' int_param_1 = ', I3)", &
               run_parameters % scan_data_array(i) % int_param_1
          print "(' int_param_2 = ', I3)", &
               run_parameters % scan_data_array(i) % int_param_2
          print "(' int_param_3 = ', I3)", &
               run_parameters % scan_data_array(i) % int_param_3
          print "(' real_param_1 = ', en11.2)", &
               run_parameters % scan_data_array(i) % real_param_1
          print "(' real_param_2 = ', en11.2)", &
               run_parameters % scan_data_array(i) % real_param_2
          print "(' logic_param_1 = ', L1)", &
               run_parameters % scan_data_array(i) % logic_param_1
          print 20
          if (interface_parameters % with_pauses) then
             print*, '(press enter to continue)'
             read 30, junk
          end if
       end do

       do i = 1, run_parameters % num_movie_namelists
          print "(' movie# = ', I2)", i
          print "(' ----------')"
          print "(' movie_type = ', I2)", &
               run_parameters % movie_data_array(i) % movie_type
          print "(' time = ', en15.8)", &
               run_parameters % movie_data_array(i) % time
          print "(' num_frames = ', I5)", &
               run_parameters % movie_data_array(i) % num_frames
          print "(' scale = ', en15.8)", &
               run_parameters % movie_data_array(i) % scale
          print 20
          if (interface_parameters % with_pauses) then
             print*, '(press enter to continue)'
             read 30, junk
          end if
       end do
    end if

!*****************  echo namelist output_data  ******************!
    if (output_parameters % echo_output) then
       print*, 'namelist output_data:'
       print*, '---------------------'
       print "(' echo_output = ', L1)", output_parameters % echo_output
       print "(' clean_outputs = ', L1)", output_parameters % clean_outputs
       print "(' file_headers = ', L1)", output_parameters % file_headers
       print "(' plot_dispersion = ', L1)", output_parameters % plot_dispersion
       print "(' power_out_vs_power_in = ', L1)", &
            output_parameters % power_out_vs_power_in
       print "(' gain_vs_power_in = ', L1)", &
            output_parameters % gain_vs_power_in
       print "(' power_out_vs_phase_in = ', L1)", &
            output_parameters % power_out_vs_phase_in
       print "(' phase_out_vs_power_in = ', L1)", &
            output_parameters % phase_out_vs_power_in
       print "(' phase_out_vs_phase_in = ', L1)", &
            output_parameters % phase_out_vs_phase_in
       print "(' power_out_vs_freq = ', L1)", &
            output_parameters % power_out_vs_freq
       print "(' gain_vs_freq = ', L1)", &
            output_parameters % gain_vs_freq
       print "(' phase_vs_freq = ', L1)", &
            output_parameters % phase_vs_freq
       print "(' num_axial_points =' , I4)", &
            output_parameters % num_axial_points
       print "(' circuit_power_vs_z = ', L1)", &
            output_parameters % circuit_power_vs_z
       print "(' magnitude_vs_z = ', L1)", output_parameters % magnitude_vs_z
       print "(' mvz_qty = ', I1)", output_parameters % mvz_qty
       print "(' phase_vs_z = ', L1)", output_parameters % phase_vs_z
       print "(' pvz_qty = ', I1)", output_parameters % pvz_qty
       print "(' phase_differences_vs_z = ', L1)", &
            output_parameters % phase_differences_vs_z
       print "(' pd_qty_1 = ', I1)", output_parameters % pd_qty_1
       print "(' pd_qty_2 = ', I1)", output_parameters % pd_qty_2
       print "(' hot_phase_velocity_vs_z = ', L1)", &
            output_parameters % hot_phase_velocity_vs_z
       print "(' disk_orbits_vs_z = ', L1)", &
            output_parameters % disk_orbits_vs_z
       print "(' conserved_quantity_vs_z = ', L1)", &
            output_parameters % conserved_quantity_vs_z
       print "(' beam_energy_vs_z = ', L1)", &
            output_parameters % beam_energy_vs_z
       print "(' dc_beam_vel_vs_z = ', L1)", &
            output_parameters % dc_beam_vel_vs_z
       print "(' phase_space = ', L1)", &
            output_parameters % phase_space

       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist frequency_list  ******************!
    if (frequency_parameters % echo_frequencies) then
       print*, 'namelist frequency_list:'
       print*, '------------------------'
       print "(' echo_frequencies = ', L1)", &
            frequency_parameters % echo_frequencies
       print "(' num_input_freqs = ', I3)", &
            frequency_parameters % num_input_freqs
       print "(' base_frequency = ', en11.2)", &
            frequency_parameters % base_frequency
       print "(' read_from_file = ', L1)", &
            frequency_parameters % read_from_file
       print "(' frequency_integer, power_input, phase_input:')"
       do i = 1, frequency_parameters % num_input_freqs
          if (units_structure % input_power_dBm) then
             tmpPower = 10.0*log10(frequency_parameters % power_input(i)/0.001)
          else
             tmpPower = frequency_parameters % power_input(i)
          end if
          if (units_structure % phase_degrees) then
             tmpPhase = (180.0/pi)*frequency_parameters % phase_input(i)
          else
             tmpPhase = frequency_parameters % phase_input(i)
          end if
          print "(t7, I5, t20, en11.2, t33, en11.2)", &
               frequency_parameters % frequency_integer(i), &
               tmpPower, &
               tmpPhase
       end do
       print "(' highest_order_IMP = ', I3)", &
            frequency_parameters % highest_order_IMP
       print "(' min_ckt_freq = ', I5)", &
            frequency_parameters % min_ckt_freq
       print "(' max_ckt_freq = ', I5)", &
            frequency_parameters % max_ckt_freq
       print "(' min_space_charge_freq = ', I5)", &
            frequency_parameters % min_space_charge_freq
       print "(' max_space_charge_freq = ', I5)", &
            frequency_parameters % max_space_charge_freq
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist circuit  ******************!
    if (circuit_parameters % echo_circuit) then
       print*, 'namelist circuit:'
       print*, '-----------------'
       print "(' echo_circuit = ', L1)", circuit_parameters % echo_circuit
       if (units_structure % length_cm) then
          print "(' circuit_length = ', en11.2)", &
               100.0*circuit_parameters % circuit_length
       else
          print "(' circuit_length = ', en11.2)", &
               circuit_parameters % circuit_length
       end if
       print "(' number_ckt_sections = ', I2)", &
            circuit_parameters % number_ckt_sections
       print "(' read_from_file = ', L1)", circuit_parameters % read_from_file
       print "(' #,    location, tape_width, helix_radius, helix_pitch, vane_radius, wall_radius')"
       do i= 1, circuit_parameters % number_ckt_sections
          if (units_structure % length_cm) then
             tmpCSL = circuit_parameters % ckt_section_location(i)*100.0
             tmpTW = circuit_parameters % tape_width(i)*100.0
             tmpHR = circuit_parameters % helix_radius(i)*100.0
             tmpHP = circuit_parameters % helix_pitch(i)*100.0
             tmpVR = circuit_parameters % vane_radius(i)*100.0
             tmpWR = circuit_parameters % wall_radius(i)*100.0
          else
             tmpCSL = circuit_parameters % ckt_section_location(i)
             tmpTW = circuit_parameters % tape_width(i)
             tmpHR = circuit_parameters % helix_radius(i)
             tmpHP = circuit_parameters % helix_pitch(i)
             tmpVR = circuit_parameters % vane_radius(i)
             tmpWR = circuit_parameters % wall_radius(i)
          end if

          print "(i2, t5, en11.2, t17, en11.2, t31, en11.2, t44, en11.2, t57, en11.2, t70, en11.2)", &
               i, &
               tmpCSL, &
               tmpTW, &
               tmpHR, &
               tmpHP, &
               tmpVR, &
               tmpWR
       end do
       print "(' input_TL_impedance = ', en11.2)", &
            circuit_parameters % input_TL_impedance
       print "(' output_TL_impedance = ', en11.2)", &
            circuit_parameters % output_TL_impedance
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist dispersion  ******************!
    if (dispersion_parameters % echo_dispersion) then
       print*, 'namelist dispersion:'
       print*, '--------------------'
       print "(' echo_dispersion = ', L1)", &
            dispersion_parameters % echo_dispersion
       print "(' vph_over_c = ', L1)", &
            dispersion_parameters % vph_over_c
       print "(' num_dispersion_freqs = ', I3)", &
            dispersion_parameters % num_dispersion_freqs
       print "(' read_from_file = ', L1)", &
            dispersion_parameters % read_from_file
       do i = 1, circuit_parameters % number_ckt_sections
          print "(' circuit section #', I2)", i
          print "(t2,'disp_freq_int, phase_velocity, impedance, space_charge_redux')"
          do j = 1, dispersion_parameters % num_dispersion_freqs
             tmpVPH = dispersion_parameters % phase_velocity(i,j)
             if (units_structure % length_cm) then
                tmpVPH = dispersion_parameters % phase_velocity(i,j) * 100.0
             end if
             if (dispersion_parameters % vph_over_c) then
                tmpVPH = dispersion_parameters % phase_velocity(i,j) / c
             end if
             !convert loss, have to check three cases.             

             print "(t5,i5,t17,en11.2,t31,en11.2,t43,en11.2)", &
                  dispersion_parameters % dispersion_freq_integer(i,j), &
                  tmpVPH, &
                  dispersion_parameters % impedance(i,j), &
                  dispersion_parameters % space_charge_redux(i,j)
          end do
          if (dispersion_parameters % num_dispersion_freqs >= 3 &
               .and. interface_parameters % with_pauses) then
             print*, '(press enter to continue)'
             read 30, junk
          end if
       end do
       print "(' intrplt_sections = ', L1)", &
            dispersion_parameters % intrplt_sections
       do i = 1, circuit_parameters % number_ckt_sections - 1
          if (dispersion_parameters % interpol_sects_list(i) /= -1) then
             print "('   interpol_sects_list(', I2,') =',I3)", i, &
                  dispersion_parameters % interpol_sects_list(i)
          end if
       end do

       print "(' intrplt_over_length = ', L1)", &
            dispersion_parameters % intrplt_over_length
       if (units_structure % length_cm) then
          print "(' interpol_length = ', en11.2)", &
            dispersion_parameters % interpol_length * 100.0
       else
          print "(' interpol_length = ', en11.2)", &
            dispersion_parameters % interpol_length
       end if
       print "(' use_sheath_model = ', L1)", &
            dispersion_parameters % use_sheath_model
       print "(' use_tape_model = ', L1)", &
            dispersion_parameters % use_tape_model
       print "(' use_antonsen_formula = ', L1)", &
            dispersion_parameters % use_antonsen_formula
       print "(' use_correction = ', L1)", &
            dispersion_parameters % use_correction
       print "(' use_beta_c = ', L1)", &
            dispersion_parameters % use_beta_c
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist loss  ******************!
    if (loss_parameters % echo_losses) then
       print*, 'namelist losses:'
       print*, '----------------'
       print "(' echo_losses = ', L1)", loss_parameters % echo_losses
       print "(' number_loss_locations = ', I3)", &
            loss_parameters % number_loss_locations
       print "(' number_loss_freqs = ', I3)", &
            loss_parameters % number_loss_freqs
       print "(' read_from_file = ', L1)", &
            loss_parameters % read_from_file

       do i = 1, loss_parameters % number_loss_locations
          tmpLOC = loss_parameters % loss_location(i)
          if (units_structure % length_cm) then
             tmpLOC = tmpLOC * 100.0
          end if
          print "(' loss_location = ', en11.2)", &
               tmpLOC
         print "(t2,'loss_freq_int',t22,'loss')"
          do j = 1, loss_parameters % number_loss_freqs
             !convert loss, have to check three cases.             
             ! input is Np/m
             tmpAlpha = loss_parameters % loss(i,j)
             ! input is dB/cm
             if (.not. units_structure % nepers .and. &
                  units_structure % length_cm) then
                tmpAlpha = loss_parameters % loss(i,j) * 8.68589 / 100.0
             ! input is dB/m
             else if (.not. units_structure % nepers .and. &
                  .not. units_structure % length_cm) then
                tmpAlpha = loss_parameters % loss(i,j) * 8.68589
             ! input is Np/cm
             else if (units_structure % nepers .and. &
                  units_structure % length_cm) then
                tmpAlpha = loss_parameters % loss(i,j) / 100.0
             end if

             print "(t5,i5,t17,en11.2,t31,en11.2,t43,en11.2,t60,en11.2)", &
                  loss_parameters % loss_freq_integer(i,j), &
                  tmpAlpha
          end do
          if (loss_parameters % number_loss_freqs >= 3 &
               .and. interface_parameters % with_pauses) then
             print*, '(press enter to continue)'
             read 30, junk
          end if
       end do

       print "(' intrplt_btwn_points = ', L1)", &
            loss_parameters % intrplt_btwn_points
       print "(' intrplt_over_length = ', L1)", &
            loss_parameters % intrplt_over_length

       if (units_structure % length_cm) then
          print "(' interpol_length = ', en11.2)", &
               loss_parameters % interpol_length * 100.0
       else
          print "(' interpol_length = ', en11.2)", &
               loss_parameters % interpol_length
       end if
       print "(' use_loss_model = ', L1)", &
            loss_parameters % use_loss_model

       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist beam  ******************!
    if (beam_parameters % echo_beam) then
       print*, 'namelist beam:'
       print*, '--------------'
       print "(' echo_beam = ', L1)", beam_parameters % echo_beam
       print "(' relativistic = ', L1)", beam_parameters % relativistic
       print "(' DC_depression = ', L1)", beam_parameters % DC_depression
       print "(' voltage = ', en11.2)", beam_parameters % voltage
       print "(' current = ', en11.2)", beam_parameters % current
       if (units_structure % length_cm) then
          print "(' outer_radius = ', en11.2)", &
               beam_parameters % outer_radius * 100.0
          print "(' inner_radius = ', en11.2)", &
               beam_parameters % inner_radius * 100.0
       else
          print "(' outer_radius = ', en11.2)", beam_parameters % outer_radius
          print "(' inner_radius = ', en11.2)", beam_parameters % inner_radius
       end if
       print "(' specify_velocity = ', L1)", &
            beam_parameters % specify_velocity
       if (units_structure % length_cm) then
          print "(' beam_velocity = ', en11.2)", &
               beam_parameters % beam_velocity * 100.0
       else
          print "(' beam_velocity = ', en11.2)", &
               beam_parameters % beam_velocity
       end if
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

!*****************  echo namelist numerical_data  ******************!
    if (numerical_parameters % echo_numerical) then
       print*, 'namelist numerical_data:'
       print*, '------------------------'
       print "(' num_grid_pts = ', I4)", numerical_parameters % num_grid_pts
       print "(' base_number_disks = ', I4)", &
            numerical_parameters % base_number_disks
       print "(' fnt0 = ', en11.2)", numerical_parameters % fnt0
       print "(' adaptive_step = ', L1)", numerical_parameters % adaptive_step
       print "(' tolerance = ', en11.2)", numerical_parameters % tolerance
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if
  end subroutine echo_parameters

end module
