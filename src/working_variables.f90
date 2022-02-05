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

module working_variables
  use precision95
  use parameters

  implicit none

  !! derived qtys for beam ... and more?
  type derived_qtys_S
     real(dp) :: u0
     real(dp) :: rho0
     real(dp) :: beam_area
  end type derived_qtys_S
  
  !*** variable definitions ***!
  type (derived_qtys_S) derived_qtys
  
  !normalized fundamental frequency
  real(dp) :: w0
  !number of disks
  integer :: N_disks

  !frequency sets OmegaAlpha and integer array fl
  integer, dimension(:,:), allocatable :: OmegaAlpha
  integer, dimension(:), allocatable :: fl
  integer, dimension(:,:), allocatable :: pair_matrix !for MUSE, S-MUSE

  !dispersion matrices
  real(dp), dimension(:,:), allocatable :: vph_matrix
  real(dp), dimension(:,:), allocatable :: nvph_matrix !normalized vph
  real(dp), dimension(:,:), allocatable :: K_matrix
  real(dp), dimension(:,:), allocatable :: alpha_matrix
  real(dp), dimension(:,:), allocatable :: scrf_matrix
  real(dp), dimension(:,:), allocatable :: nscrf_matrix !normalized scrf
  !ckt parameter matrices, normalized
  real(dp), dimension(:,:), allocatable :: R_matrix
  real(dp), dimension(:,:), allocatable :: L_matrix
  real(dp), dimension(:,:), allocatable :: G_matrix
  real(dp), dimension(:,:), allocatable :: C_matrix
  real(dp), dimension(:,:), allocatable :: pC_matrix !! pierce C
  !complex phase factor, Z_\ell(z)
  complex(dp), dimension(:,:), allocatable :: Z_factor_matrix

  !3-d array to keep plot data in
  complex(dp), dimension(:,:,:), allocatable :: plot_cube

  !vectors used in RK integrator
  complex(dp), dimension(:), allocatable :: yk_1
  complex(dp), dimension(:), allocatable :: yk_2
  complex(dp), dimension(:), allocatable :: yk_3
  complex(dp), dimension(:), allocatable :: k_1
  complex(dp), dimension(:), allocatable :: k_2
  complex(dp), dimension(:), allocatable :: k_3
  complex(dp), dimension(:), allocatable :: k_4
  complex(dp), dimension(:), allocatable :: y_nminus1
  !! for dormand prince (4,5)
  complex(dp), dimension(:,:), allocatable :: k_45
  real(dp), dimension(7) :: c_45
  real(dp), dimension(7) :: b_45
  real(dp), dimension(7) :: bhat_45
  real(dp), dimension(7,7) :: a_45
  !vectors used in F()
  complex(dp), dimension(:), allocatable :: netE
  complex(dp), dimension(:), allocatable :: rhol
  real(dp), dimension(:), allocatable :: exp_sum_real
  real(dp), dimension(:), allocatable :: exp_sum_imag
  real(dp), dimension(:,:), allocatable :: exp_real
  real(dp), dimension(:,:), allocatable :: exp_imag
  real(dp), dimension(:,:), allocatable :: disk_matrix
  !MUSE LHS matrix and RHS vector
  complex(dp), dimension(:,:), allocatable :: thesis_V
  complex(dp), dimension(:), allocatable :: thesis_w
  integer, dimension(:), allocatable :: ipiv

  logical :: vs_z_plot = .false. !make true if any of the vs_z's true
  logical :: vs_freq_plot = .false. !make true if any of the vs_freqs's true
  logical :: vs_powin_plot = .false. !make true if any of the vs_powin's true

  logical :: dcconst = .false. !flags ODE routines for LATTE dc_constant

  logical :: run_one = .true. !true if any runtime_one, _two, or _three true
  logical :: run_two = .false. !true if any runtime_two, or _three true
  logical :: run_three = .false. !true if any runtime_three true

  !some other variables
  integer :: M !total number of positive frequencies
  integer :: num_ckt_freqs !number of circuit frequencies
  integer :: num_spch_freqs !number of circuit frequencies
  integer :: first_ckt_index !index of lowest ckt freq
  integer :: last_ckt_index !index of highest ckt freq
  integer :: first_spch_index !index of lowest space charge freq
  integer :: last_spch_index !index of highest space charge freq
  character :: modelID !flag used in F() to know which derivative to compute
  

contains
  !******************************************************************!
  !********************  initialization routines  *******************!
  !******************************************************************!
  
  !*******************  compute derived quantities  *****************!
  subroutine compute_derived_qtys
    implicit none
    integer tmpDisk
    character junk

    !! beam velocity. if beam velocity specified assign it, else compute it
    if (beam_parameters % specify_velocity) then
       derived_qtys % u0 = beam_parameters % beam_velocity
       !! need to compute voltage to go along with it so that pierce
       !! parameter is consistent
       beam_parameters % voltage = 0.5 * me * (derived_qtys % u0)**2 &
            / echarge
    else
       derived_qtys % u0 = sqrt(2.0 * echarge * beam_parameters % voltage &
            / me)
    end if

    !! beam area
    derived_qtys % beam_area = pi * ((beam_parameters % outer_radius)**2 &
         - (beam_parameters % inner_radius)**2)
    !! dc charge density
    derived_qtys % rho0 = beam_parameters % current &
         / (derived_qtys % u0 * derived_qtys % beam_area)
    

    !compute normalized fundamental frequency
    w0 = 2*pi*frequency_parameters % base_frequency &
         * circuit_parameters % circuit_length / derived_qtys % u0
    ! compute N_disks
    ! formula is fnt0 * num_input_freqs * number_harmonics 
    ! * base_number_disks 
    ! where number_harmonics is estimated by
    ! max_space_charge_freq / 1st drive frequency
    N_disks = numerical_parameters % fnt0 * num_spch_freqs &
         * numerical_parameters % base_number_disks &
         * (frequency_parameters % max_space_charge_freq &
         / frequency_parameters % frequency_integer(1))

    ! check if N_disks is even or odd and change it to odd
    tmpDisk = N_disks
    do while (tmpDisk > 0)
       tmpDisk = tmpDisk - 2
    end do
    if (tmpDisk == 0) then
       N_disks = N_disks - 1
    end if
    
    if (interface_parameters % echo_initialization) then
       print*, 'Computing derived quantities:'
       print "(' w0 = ', f8.3, '  (normalized base frequency)')", w0
       print "(' N_disks = ', I6, '  (total number of disks)')", N_disks
       if (units_structure % length_cm) then
          print "(' u0 =',en11.3,' (cm/s)')", derived_qtys % u0 * 100.0
          print "(' beam_area = ',en11.3,' (cm^2)')", &
               derived_qtys % beam_area * (100.0)**2
          print "(' rho0 = ',en11.3,' (C/cm^3)')", &
               derived_qtys % rho0 / (100.0**3)
       else
          print "(' u0 =',en11.3,' (m/s)')", derived_qtys % u0
          print "(' beam_area = ',en11.3,' (m^2)')", derived_qtys % beam_area
          print "(' rho0 = ',en11.3,' (C/m^3)')", derived_qtys % rho0
       end if
20     format (/)
30     format(A)
       print 20
       if (interface_parameters % with_pauses) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if
  end subroutine compute_derived_qtys
  
  !*******************  create frequency array fl  ******************!
  subroutine create_frequency_array(echo)
    implicit none
    logical, intent(in) :: echo

    integer, dimension(:,:), allocatable :: tmpMatrix
    integer, dimension(:,:), allocatable :: tmpArray !make 2-d so can reshape
    integer OA_size, tM_rows, tM_columns, HOI !!highest_order_IMP
    integer i,j,k
    character junk
    
    if (interface_parameters % echo_initialization .and. echo) then
       print*, ''
       print*, 'Creating frequency array'
    end if
    
    HOI = frequency_parameters % highest_order_IMP
    
    !*** compute OmegaAlphas ***!
    ! allocate OmegaAlpha
    allocate(OmegaAlpha(HOI, -max_freqs:max_freqs))
    ! set OmegaAlpha to 0's
    do i = 1, HOI
       do j = -max_freqs, max_freqs
          OmegaAlpha(i,j) = 0
       end do
    end do
    
    !fill in OmegaAlpha(1,:) with freq_integer
    do j = 1, frequency_parameters % num_input_freqs
       OmegaAlpha(1,j) = frequency_parameters % frequency_integer(j)
       OmegaAlpha(1,-j) = -frequency_parameters % frequency_integer(j)
    end do
    
    ! OA_size keeps track of maximum number of nonzero
    ! frequencies in OmegaAlpha, including negative frequencies
    OA_size = 2 * frequency_parameters % num_input_freqs
    
    !! compute the OmegaAlpha
    do i = 2, HOI
       !tM_rows is #rows in tmpMatrix. originally written for i-1, but
       !it seems i/2 is better since uses symmetry.
       tM_rows = i/2
       !tM_rows = i-1
       !tM_columns is #columns in tmpMatrix. lowest upper bound I have been
       !able to figure is OA_size**2. should be something smaller using
       !structure of frequency sets, but I haven't figured it.
       tM_columns = OA_size**2
       ! allocate tmpMatrix
       allocate(tmpMatrix(tM_rows,tM_columns))
       !make tmpMatrix all zeros
       tmpMatrix = &
            reshape(source=(/(0, k = 1, tM_rows * tM_columns)/), &
            shape = (/tM_rows , tM_columns/))
       
       do j = 1, tM_rows !! j loop fills tmpMatrix, with OmegaAlpha as source
          !really only need to go up to (i-1)/2 since j's above this are
          !repeats...but fix that later
          call outer_product(tmpMatrix,OmegaAlpha,i,j)
       end do
       
       !flatten tmpMatrix into tmpArray
       allocate(tmpArray(1,tM_rows*tM_columns))
       tmpArray = reshape(source=tmpMatrix, shape=(/1, tM_rows*tM_columns/))
       
       !! sort tmpArray.
       call sort_tmpArray(tmpArray, tM_rows * tM_columns)
       
       !! make new OmegaAlpha from tmpArray
       call new_OmegaAlpha(tmpArray,tM_rows*tM_columns,OmegaAlpha,OA_size,i)
       
       !! kill tmpMatrix and tmpArray
       deallocate(tmpMatrix)
       deallocate(tmpArray)
    end do !! end construction of OmegaAlpha
    
    !*** make fl from OmegaAlpha ***!
    call make_fl(OmegaAlpha, M)
    
    !*** count circuit frequencies and space charge freqs ***!
    call count_frequencies

    if (interface_parameters % echo_initialization .and. echo) then
       print*, ''
       print "(' A total of ', I3, ' frequencies (',I3,' positive) for', I3, ' drive frequencies')", &
            2*M+1, M, frequency_parameters % num_input_freqs
       print "(' with highest order IMP =', I3)", HOI
       
       print*, ''
       print "(' List of frequencies including all IMPs up to order', I3,':')", HOI
       print*,(fl(i),i=-M,M)
       print 20
       print "(' min_space_charge_freq = ', I5, t44,'min_ckt_freq = ', I5)", &
            frequency_parameters % min_space_charge_freq, &
            frequency_parameters % min_ckt_freq
       print "(t11,'max_ckt_freq = ', I5, t35,'max_space_charge_freq = ', I5)", &
            frequency_parameters % max_ckt_freq, &
            frequency_parameters % max_space_charge_freq
       if (interface_parameters % echo_initialization .and. echo) then
          print "(t10,'num_ckt_freqs = ', I5,t42,'num_spch_freqs = ', I5)", &
               num_ckt_freqs, num_spch_freqs
       end if
20     format (/)
30     format(A)
       print 20
       if (interface_parameters % with_pauses .and. echo) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if

    deallocate(OmegaAlpha)

  contains
    subroutine outer_product(tmpMatrix,OmegaAlpha,i,j)
      implicit none
      integer, dimension(:,:), intent(inout) :: tmpMatrix
      integer, dimension(HOI,-max_freqs:max_freqs), intent(in) :: &
           OmegaAlpha
      integer, intent(in) :: i, j
      integer k, l, m, r, s, tmpFreq
      logical not_in
      
      m = 1 !! m tracks index going into tmpArray
      do k = -max_freqs, max_freqs
         do l = -max_freqs, max_freqs
            if (OmegaAlpha(j,k) == 0 .or. OmegaAlpha(i-j,l) == 0) then
               tmpFreq = 0
            else
               tmpFreq = OmegaAlpha(j,k) + OmegaAlpha(i-j,l)
            end if
            !! first check that it isn't already in tmpArray
            not_in = .true.
            do r = 1, j
               do s = 1, tM_columns
                  if (tmpMatrix(r,s) == tmpFreq) then
                     not_in = .false.
                  end if
               end do
            end do
            if (tmpFreq /= 0 .and. not_in) then
               tmpMatrix(j,m) = tmpFreq
               m = m + 1
            end if
         end do
      end do
    end subroutine outer_product
    
    subroutine new_OmegaAlpha(tmpArray, N, OmegaAlpha, OA_size, i)
      integer, dimension(:,:), intent(in) :: tmpArray
      integer, intent(in) :: N
      integer, dimension(HOI,-max_freqs:max_freqs), intent(inout) :: &
           OmegaAlpha
      integer, intent(inout) :: OA_size
      integer, intent(in) :: i
      integer j, k
      
      !! find location in tmpArray of first greater than zero entry
      k = 1
      do while (tmpArray(1,k) <= 0)
         k = k + 1
      end do
      
      !! from k start putting elements of tmpArray into OmegaAlpha
      j = 1
      do while (k <= N)
         if (j > max_freqs) then
            print*, ''
            print*, 'Too many frequencies. Reduce highest_order_IMP or increase max_freqs.'
            print*, 'Stopping.'
            stop
         end if
         OmegaAlpha(i,j) = tmpArray(1,k)
         OmegaAlpha(i,-j) = -tmpArray(1,k)
         k = k + 1
         j = j + 1
      end do
      
      !set the new size of OA_size
      if (2*(j-1) > OA_size) then
         OA_size = 2*(j-1)
      end if
    end subroutine new_OmegaAlpha
    
    subroutine make_fl(OmegaAlpha, M)
      integer, dimension(HOI,-max_freqs:max_freqs), intent(in) :: &
           OmegaAlpha
      integer, intent(inout) :: M
      integer i,j,k,l,tmpFreq
      logical not_in
      
      !!flatten all of OmegaAlpha into tmpArray, don't use reshape b/c
      !!want to get rid of dupes
      allocate(tmpArray(1,HOI*(2*max_freqs+1)))
      !!initialize tmpArray to zeros
      tmpArray = reshape((/(0,i=1,HOI*(2*max_freqs+1))/),&
           shape=(/1,HOI*(2*max_freqs+1)/))
      k = 1
      do i = 1, HOI
         do j = -max_freqs, max_freqs
            tmpFreq = OmegaAlpha(i,j)
            not_in = .true.
            do l = 1, k
               if (tmpArray(1,l) == tmpFreq) then
                  not_in = .false.
               end if
            end do
            if (not_in) then
               tmpArray(1,k) = tmpFreq
               k = k + 1
            end if
         end do
      end do
      M = (k-1)/2 !! total number of positive frequencies, global variable!
      
      !!sort tmp_Array
      call sort_tmpArray(tmpArray,HOI*(2*max_freqs+1))
      
      !find first greater than zero entry
      k = 1
      do while (tmpArray(1,k) <= 0)
         k = k + 1
      end do
      
      !allocate space for fl
      allocate(fl(-M:M))
      
      do i = 0, M
         fl(i) = 0
         fl(-i) = 0
      end do

      !fill in fl
      j = 1
      do while (k <= HOI*(2*max_freqs+1))
         if (j > max_freqs) then
            print*, ''
            print*, 'Too many frequencies. Reduce highest_order_IMP or increase max_freqs.'
            print*, 'Stopping.'
            stop
         end if
         fl(j) = tmpArray(1,k)
         fl(-j) = -tmpArray(1,k)
         k = k + 1
         j = j + 1
      end do
    end subroutine make_fl
    
    
    !! bubble sort taken from Hahn, F90 for scientists and engineers
    subroutine sort_tmpArray(tmpArray, N)
      integer, dimension(:,:), intent(inout) :: tmpArray
      integer, intent(in) :: N
      integer :: tmp, j, k
      logical :: sorted
      
      sorted = .false.
      k = 0
      
      do while (.not. sorted)
         sorted = .true.
         k = k + 1
         do j = 1, N-k
            if (tmpArray(1,j) > tmpArray(1,j+1)) then
               tmp = tmpArray(1,j)
               tmpArray(1,j) = tmpArray(1,j+1)
               tmpArray(1,j+1) = tmp
               sorted = .false.
            end if
         end do
      end do
    end subroutine sort_tmpArray
  end subroutine create_frequency_array

  subroutine create_harmonic_array(echo)
    implicit none
    logical, intent(in) :: echo

    integer :: i
    character :: junk

    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, 'USING CREATE_HARMONIC_ARRAY'
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    print*, ''
    pause

    allocate(fl(-frequency_parameters % highest_order_IMP : &
         frequency_parameters % highest_order_IMP))
    
    M = frequency_parameters % highest_order_IMP

    fl(0) = 0
    do i = 1, frequency_parameters % highest_order_IMP
       fl(i) = i * frequency_parameters % frequency_integer(1)
       fl(-i) = -fl(i)
    end do

    !*** count circuit frequencies and space charge freqs ***!
    call count_frequencies

    if (interface_parameters % echo_initialization .and. echo) then
       print*, ''
       print "(' A total of ', I3, ' frequencies (',I3,' positive) for', I3, ' drive frequencies')", &
            2*M+1, M, frequency_parameters % num_input_freqs
       print "(' with highest order IMP =', I3)", &
            frequency_parameters % highest_order_IMP
       
       print*, ''
       print "(' List of frequencies including all IMPs up to order', I3,':')",frequency_parameters % highest_order_IMP
       print*,(fl(i),i=-M,M)
       print 20
       print "(' min_space_charge_freq = ', I5, t44,'min_ckt_freq = ', I5)", &
            frequency_parameters % min_space_charge_freq, &
            frequency_parameters % min_ckt_freq
       print "(t11,'max_ckt_freq = ', I5, t35,'max_space_charge_freq = ', I5)", &
            frequency_parameters % max_ckt_freq, &
            frequency_parameters % max_space_charge_freq
       if (interface_parameters % echo_initialization .and. echo) then
          print "(t10,'num_ckt_freqs = ', I5,t42,'num_spch_freqs = ', I5)", &
               num_ckt_freqs, num_spch_freqs
       end if
20     format (/)
30     format(A)
       print 20
       if (interface_parameters % with_pauses .and. echo) then
          print*, '(press enter to continue)'
          read 30, junk
       end if
    end if
  end subroutine create_harmonic_array

  subroutine count_frequencies
    integer i
    num_ckt_freqs = 0
    num_spch_freqs = 0
      
    !! count the number of circuit and space charge frequencies
    !! and assign the first_ckt_index and first_spch_index
    do i = 1, M
       if (frequency_parameters % min_ckt_freq <= fl(i) &
            .and. fl(i) <= frequency_parameters % max_ckt_freq) then
          num_ckt_freqs = num_ckt_freqs + 1
          if (num_ckt_freqs == 1) then
             first_ckt_index = i
          end if
       end if
       if (frequency_parameters % min_space_charge_freq <= fl(i) &
            .and. fl(i) <= frequency_parameters % max_space_charge_freq) then
          num_spch_freqs = num_spch_freqs + 1
          if (num_spch_freqs == 1) then
             first_spch_index = i
          end if
       end if
    end do
    
    !! assign the last_ckt_index and last_spch_index
    last_ckt_index = first_ckt_index + num_ckt_freqs - 1
    last_spch_index = first_spch_index + num_spch_freqs - 1
  end subroutine count_frequencies

  subroutine make_pair_matrix
    integer :: i, j, kk, num_pairs

    allocate(pair_matrix(-M:M,1:2*(2*M+1)+2))
    do i = -M, M
       do j = 1, 2*(M+2)
          pair_matrix(i,j) = 0
       end do
    end do

    do i = -M, M
       pair_matrix(i,1) = i
       num_pairs = 0
       do j = -M, M
          kk = find_freq(fl(i) - fl(j))
          if (fl(kk)+fl(j)==fl(i)) then
             num_pairs = num_pairs+1
             pair_matrix(i,2*num_pairs+1) = j
             pair_matrix(i,2*num_pairs+2) = kk
             !print*, 2*num_pairs+1, 2*num_pairs+2
             !print*, 'i=',i,'j=',j,'kk=',kk, 'np=', num_pairs
          end if
       end do
       pair_matrix(i,2) = num_pairs
    end do

    do i = -M, M
       !print*, 'i=',i,'num_pairs=',pair_matrix(i,2)
       do j = 1, pair_matrix(i,2)
          !print*, 'fl(i)=',fl(i),'fl(j)=',fl(pair_matrix(i,2*j+1)), &
          !'fl(kk)=',fl(pair_matrix(i,2*j+2))
       end do
    end do
  end subroutine make_pair_matrix

  function find_freq(freq_int)
    integer find_freq
    integer, intent(in) :: freq_int

    integer :: i

    find_freq = 0
    do i = -M, M
       if (fl(i) == freq_int) then
          find_freq = i
       end if
    end do
  end function find_freq
  !*******************  initialize plot arrays  *********************!
  
  !****************  initialize dispersion arrays  ******************!
  subroutine initialize_dispersion_arrays
    character junk
    
    !! allocate dispersion arrays of size
    !!(num_ckt_sections, num_input_freqss, num_grid_pts)
    !! second 0 index corresponds to z = 0.0
    allocate(vph_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(nvph_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(K_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(alpha_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(scrf_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(nscrf_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(Z_factor_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    
    !! compute initial dispersion
    if (interface_parameters % echo_initialization) then
       print*, 'Computing vph, K, alpha, and scrf matrices'
    end if
    call compute_dispersion_arrays
    
    !! plot dispersion if requested
    if (output_parameters % plot_dispersion) then
       if (interface_parameters % echo_initialization) then
          print*, 'Plotting vph, K, alpha and scrf data'
       end if
       call plot_dispersion_data
    end if
    
    !! now get derived circuit quantities
    allocate(R_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(L_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(G_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(C_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    allocate(pC_matrix(-M:M, 0:numerical_parameters % num_grid_pts))
    
    if (interface_parameters % echo_initialization) then
       print*, 'Computing derived circuit quantities R, L, G, C, and pierce C'
    end if
    call derived_ckt_qtys(R_matrix, L_matrix, G_matrix, C_matrix,&
         pC_matrix)
    
20  format (/)
30  format(A)
    print 20
    if (interface_parameters % with_pauses) then
       print*, '(press enter to continue)'
       read 30, junk
    end if
  end subroutine initialize_dispersion_arrays

  !! initialize A matrix, S matrix, jacobian, frequency_pair matrix
  
  !********************  dispersion routines  *******************!
  subroutine compute_dispersion_arrays
    real(dp) z_start, z_end, tmpz, tmpvph, tmpK, tmp_scrf, tmp_alpha
    integer arr_start, arr_end, i, j, k
    
    
    !***************** vph, K, R matrices ******************!
    !!set the zero frequency parameters to zero
    do k = 0, numerical_parameters % num_grid_pts
       vph_matrix(0,k) = 0.0
       nvph_matrix(0,k) = 0.0
       K_matrix(0,k) = 0.0
       alpha_matrix(0,k) = 0.0
       scrf_matrix(0,k) = 0.0
       nscrf_matrix(0,k) = 0.0
       Z_factor_matrix(0,k) = 0.0
    end do
    
    !! loop on circuit sections
    do i = 1, circuit_parameters % number_ckt_sections
       !! get length indexing
       z_start = circuit_parameters % ckt_section_location(i)
       arr_start = numerical_parameters % num_grid_pts &
            * (z_start / circuit_parameters % circuit_length)
       !! stop if first section does not start at index == 0
       if (i == 1 .and. arr_start /= 0) then
          print*, 'dispersion:'
          print*, 'First ckt_section_location is not equal to 0.0.'
          print*, 'Stopping.'
          stop
       end if
       !!if last section set end to circuit_length
       if (i == circuit_parameters % number_ckt_sections) then
          z_end = circuit_parameters % circuit_length
          arr_end = numerical_parameters % num_grid_pts
       else !!else look to next section
          z_end = circuit_parameters % ckt_section_location(i+1)
          arr_end = numerical_parameters % num_grid_pts &
               * (z_end / circuit_parameters % circuit_length)
       end if

       !!if arr_start and/or arr_end are larger than num_grid_points reset them
       if (arr_start > numerical_parameters % num_grid_pts) then
          arr_start = numerical_parameters % num_grid_pts
       end if
       if (arr_end > numerical_parameters % num_grid_pts) then
          arr_end = numerical_parameters % num_grid_pts
       end if
       
       do j = 1, M !! loop on all frequencies
          !! get the values of vph, K for this section and frequency
          if (dispersion_parameters % use_tape_model) then
             !call tape_helix(i,j,tmpvph,tmpK)
             if (j==1) then
                print*, ''
                print*, '***********************************************'
                print*, 'Tape model not implemented, using sheath model.'
                print*, '***********************************************'
                print*, ''
             end if
             call sheath_helix(i,j,tmpvph,tmpK)
          else if (dispersion_parameters % use_sheath_model) then
             if (j==1) then
                print*, ''
                print*, '***************************************************&
                     &*******'
                print*, 'WARNING: The sheath helix model has not been &
                     &extensively'
                print*, 'tested. Furthermore, since the model does not account'
                print*, 'for a barrel, it should be used for qualitatative &
                     &purposes'
                print*, 'only.'
                print*, '***************************************************&
                     &*******'
                print*, ''
             end if
             call sheath_helix(i,j,tmpvph,tmpK)
          else
             call data_helix(i,j,tmpvph,tmpK)
          end if
          
          !! get value for R for this section and frequency
          if (dispersion_parameters % use_antonsen_formula) then
             call antonsen_formula(i,j,tmp_scrf, tmpvph)
          else
             call data_scrf(i,j,tmp_scrf)
          end if
          
          !!fill in vph_matrix, K_matrix, scrf_matrix
          do k = arr_start, arr_end
             !!plus and minus frequencies
             vph_matrix(j,k) = tmpvph
             vph_matrix(-j,k) = tmpvph
             K_matrix(j,k) = tmpK
             K_matrix(-j,k) = tmpK
             scrf_matrix(j,k) = tmp_scrf
             scrf_matrix(-j,k) = tmp_scrf
          end do
       end do !! end loop on frequencies
    end do !! end loop on ckt sections
    
    !! interpolate results
    if (dispersion_parameters % intrplt_sections .or. &
         dispersion_parameters % intrplt_over_length) then
       call interpolate_dispersion(vph_matrix, K_matrix, scrf_matrix)
    end if
    
    !! fill in nvph_matrix, normalized phase velocity
    !! fill in nscrf_matrix. the purpose of this matrix is so that calls
    !! to nscrf() return normalized values without having to compute
    !! normalizations each time. however, want to leave the scrf_matrix
    !! in tact for plotting, diagnositics, etc.
    do i = 1, M
       do j = 0, numerical_parameters % num_grid_pts
          nvph_matrix(i,j) = vph_matrix(i,j) / derived_qtys % u0
          nvph_matrix(-i,j) = vph_matrix(-i,j) / derived_qtys % u0

          nscrf_matrix(i,j) = scrf_matrix(i,j) &
               * echarge * derived_qtys % rho0 &
               * (circuit_parameters % circuit_length)**2 &
               / ((derived_qtys % u0)**2 * me * eps0)
          nscrf_matrix(-i,j) = scrf_matrix(-i,j) &
               * echarge * derived_qtys % rho0 &
               * (circuit_parameters % circuit_length)**2 &
               / ((derived_qtys % u0)**2 * me * eps0)

          tmpz = real(j)/real(numerical_parameters % num_grid_pts)
          if (i >= first_ckt_index .and. i <= last_ckt_index) then
             Z_factor_matrix(i,j) = &
                  exp(smo*fl(i)*w0*(1.0 - 1.0/nvph_matrix(i,j)) * tmpz)
             Z_factor_matrix(-i,j) = Z_factor_matrix(i,j)
          else
             Z_factor_matrix(i,j) = 0.0
             Z_factor_matrix(-i,j) = 0.0
          end if
       end do
    end do
    
    !***************** alpha matrix ******************!
    
    !! loop on loss_location 
    do i = 1, loss_parameters % number_loss_locations
       !! get length indexing
       z_start = loss_parameters % loss_location(i)
       arr_start = numerical_parameters % num_grid_pts &
            * (z_start / circuit_parameters % circuit_length)
       !! stop if first section does not start at index == 0
       if (i == 1 .and. arr_start /= 0) then
          print*, 'loss:'
          print*, 'First loss_location is not equal to 0.0.'
          print*, 'Stopping.'
          stop
       end if
       !!if last section set end to circuit_length
       if (i == loss_parameters % number_loss_locations) then
          z_end = circuit_parameters % circuit_length
          arr_end = numerical_parameters % num_grid_pts
       else !!else look to next section
          z_end = loss_parameters % loss_location(i+1)
          arr_end = numerical_parameters % num_grid_pts &
               * (z_end / circuit_parameters % circuit_length)
       end if

       !!if arr_start and/or arr_end are larger than num_grid_points reset them
       if (arr_start > numerical_parameters % num_grid_pts) then
          arr_start = numerical_parameters % num_grid_pts
       end if
       if (arr_end > numerical_parameters % num_grid_pts) then
          arr_end = numerical_parameters % num_grid_pts
       end if
       
       !print*, arr_start, arr_end
       
       do j = 1, M !! loop on all frequencies
          !! get the value of alpha for this location and frequency
          if (loss_parameters % use_loss_model) then
             !call loss_model(i,j,tmp_alpha)
          end if !put "else" here if loss_model implemented
          !else
          call data_alpha(i,j,tmp_alpha)
          !end if
          
          !! fill in alpha_matrix with tmp_alpha
          !! if arr_start == arr_end then the last loss_location is equal to
          !! circuit length. in this case don't want to record loss data. if
          !! last loss_location is greater than circuit length, non-issue
          if (arr_start /= arr_end) then
             do k = arr_start, arr_end
                alpha_matrix(j,k) = tmp_alpha
                alpha_matrix(-j,k) = tmp_alpha
             end do
          end if
       end do !! end loop on frequencies
    end do !! end loop on loss locations
    
    !! interpolate results
    if (loss_parameters % intrplt_btwn_points &
         .or. loss_parameters % intrplt_over_length) then
       call interpolate_loss(alpha_matrix)
    end if
  end subroutine compute_dispersion_arrays
  subroutine plot_dispersion_data
    integer :: i, j
    real(dp) :: scale
    
    if (units_structure % length_cm) then
       scale = 100.0
    else
       scale = 1.0
    end if
    
    open(1,file='./outputs/vph.dat', action='write')
    open(2,file='./outputs/impedance.dat', action='write')
    open(3,file='./outputs/scrf.dat', action='write')
    open(4,file='./outputs/loss.dat', action='write')
    
    !! write some file headers
    if (output_parameters % file_headers) then
       write(1, fmt="(A43)") 'cold circuit phase velocity versus distance'
       write(1, fmt="(A18)") 'column = frequency'
       write(2, fmt="(A45)") 'circuit interaction impedance versus distance'
       write(2, fmt="(A18)") 'column = frequency'
       write(3, fmt="(A45)") 'space charge reduction factor versus distance'
       write(3, fmt="(A18)") 'column = frequency'
       write(4, fmt="(A35)") 'circuit attenuation versus distance'
       write(4, fmt="(A18)") 'column = frequency'
       write(1, fmt="(A11)", advance='no') 'distance'
       write(2, fmt="(A11)", advance='no') 'distance'
       write(3, fmt="(A11)", advance='no') 'distance'
       write(4, fmt="(A11)", advance='no') 'distance'
       do j = 1 , M
          if (fl(j) >= frequency_parameters % min_ckt_freq &
               .and. fl(j) <= frequency_parameters % max_ckt_freq) then
             write(1, fmt="(I15)", advance='no') fl(j)
             write(2, fmt="(I15)", advance='no') fl(j)
             write(4, fmt="(I15)", advance='no') fl(j)
          end if
          if (fl(j) >= frequency_parameters % min_space_charge_freq .and. &
               fl(j) <= frequency_parameters % max_space_charge_freq) then
             write(3, fmt="(I15)", advance='no') fl(j)
          end if
       end do
       write(1, fmt="(A1)", advance='yes')
       write(2, fmt="(A1)", advance='yes')
       write(3, fmt="(A1)", advance='yes')
       write(4, fmt="(A1)", advance='yes')
    end if

    !! write the data
    do i = 0, numerical_parameters % num_grid_pts
       write(1, fmt="(en11.2)", advance='no') i &
            * circuit_parameters % circuit_length * scale &
            /numerical_parameters % num_grid_pts
       write(2, fmt="(en11.2)", advance='no') i &
            * circuit_parameters % circuit_length * scale &
            /numerical_parameters % num_grid_pts
       write(3, fmt="(en11.2)", advance='no') i &
            * circuit_parameters % circuit_length * scale &
            /numerical_parameters % num_grid_pts
       write(4, fmt="(en11.2)", advance='no') i &
            * circuit_parameters % circuit_length * scale &
            /numerical_parameters % num_grid_pts
       !!loop on frequencies, picking out circuit frequencies
       do j = 1, M
          if (fl(j) >= frequency_parameters % min_ckt_freq &
               .and. fl(j) <= frequency_parameters % max_ckt_freq) then
             !!write phase velocity
             if (units_structure % length_cm .and. .not. &
                  dispersion_parameters % vph_over_c) then
                !cm/s
                write(1, fmt="(en20.8)", advance='no') &
                     vph_matrix(j,i) * 100.0
             else if (dispersion_parameters % vph_over_c) then
                !vph/c
                write(1, fmt="(en20.8)", advance='no') vph_matrix(j,i) / c
             else
                !m/s
                write(1, fmt="(en20.8)", advance='no') vph_matrix(j,i)
             end if
             !!write impedance
             write(2, fmt="(en20.8)", advance='no') K_matrix(j,i)
             !!write loss
             if (units_structure % length_cm .and. &
                  .not. units_structure % nepers) then
                !dB/cm
                write(4, fmt="(en20.8)", advance='no') &
                     alpha_matrix(j,i) * 8.68589 / 100.0
             else if (units_structure % length_cm .and. &
                  units_structure % nepers) then
                !Np/cm
                write(4, fmt="(en20.8)", advance='no') &
                     alpha_matrix(j,i) / 100.0
             else if (.not. units_structure % length_cm .and. &
                  .not. units_structure % nepers) then
                !dB/m
                write(4, fmt="(en20.8)", advance='no') &
                     alpha_matrix(j,i) * 8.68589
             else
                !Np/m
                write(4, fmt="(en20.8)", advance='no') &
                     alpha_matrix(j,i)
             end if
          end if
          !!write space charge reduction factor if a space charge freq
          if (fl(j) >= frequency_parameters % min_space_charge_freq .and. &
               fl(j) <= frequency_parameters % max_space_charge_freq) then
             write(3, fmt="(en20.8)", advance='no') scrf_matrix(j,i)
          end if
       end do
       write(1,fmt="(A1)")
       write(2,fmt="(A1)")
       write(3,fmt="(A1)")
       write(4,fmt="(A1)")
    end do
    close(1)
    close(2)
    close(3)
    close(4)
  end subroutine plot_dispersion_data
  
  !! compute normalized derived ckt qtys
  subroutine derived_ckt_qtys(R_matrix, L_matrix, G_matrix, C_matrix,&
       pC_matrix)
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: R_matrix
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: L_matrix
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: G_matrix
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: C_matrix
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: pC_matrix
    
    integer i, j
    complex(dp) Gamma, complex_K, X, Y
    real(dp) a
    
    !! loop on frequency
    do i = 1, M
       !! loop on distance
       do j = 0, numerical_parameters % num_grid_pts
          if (fl(i) >= frequency_parameters % min_ckt_freq &
               .and. fl(i) <= frequency_parameters % max_ckt_freq) then
             Gamma = cmplx(alpha_matrix(i,j), &
                  fl(i) * frequency_parameters % base_frequency &
                  / vph_matrix(i,j))
             complex_K = cmplx(K_matrix(i,j))
             
             !! compute X, and Y
             X = conjg(Gamma * complex_K)
             Y = conjg(Gamma / complex_K)
             
             !!assign values. these are normalized!!
             R_matrix(i,j) = real(X) &
                  * (circuit_parameters % circuit_length / K_matrix(i,j))
             L_matrix(i,j) = -(1.0/(fl(i) &
                  * frequency_parameters % base_frequency))*aimag(X) &
                  * (derived_qtys % u0 / K_matrix(i,j))
             G_matrix(i,j) = real(Y) &
                  * (circuit_parameters % circuit_length * K_matrix(i,j))
             C_matrix(i,j) = -(1.0/(fl(i) &
                  * frequency_parameters % base_frequency))*aimag(Y) &
                  * (derived_qtys % u0 * K_matrix(i,j))
             
             !!pierce parameter
             pC_matrix(i,j) = (K_matrix(i,j)*beam_parameters % current &
                  /(4.0*beam_parameters % voltage)) ** 0.333333333
             
             !print*, R_matrix(i,j), L_matrix(i,j)
             !print*, G_matrix(i,j), C_matrix(i,j)
             !print*, pC_matrix(i,j)
          else
             R_matrix(i,j) = 0.0
             L_matrix(i,j) = 0.0
             G_matrix(i,j) = 0.0
             C_matrix(i,j) = 0.0
             pC_matrix(i,j) = 0.0
          end if
       end do !! end loop on distance
    end do !! end loop on frequency
  end subroutine derived_ckt_qtys
  
  !*** functions of (z,f) from dispersion arrays ***!
  !! phase velocity function
  function vph(z,i)
    real(dp) vph
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance
    j = numerical_parameters % num_grid_pts &
         * (z / circuit_parameters % circuit_length)
    
    vph = vph_matrix(i,j)
  end function vph
  
  !! impedance function
  function K(z,i)
    real(dp) K
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance
    j = numerical_parameters % num_grid_pts &
         * (z / circuit_parameters % circuit_length)
    
    K = K_matrix(i,j)
  end function K
  
  !! space charge reduction factor function
  function scrf(z,i)
    real(dp) scrf
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance
    j = numerical_parameters % num_grid_pts &
         * (z / circuit_parameters % circuit_length)
    
    scrf = scrf_matrix(i,j)
  end function scrf
  
  !! loss Np/m function
  function alpha(z,i)
    real(dp) alpha
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance
    j = numerical_parameters % num_grid_pts &
         * (z / circuit_parameters % circuit_length)
    
    alpha = alpha_matrix(i,j)
  end function alpha
  
  !! normalized resistance R function
  function R(z,i)
    real(dp) R
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer

    !!first get integer for distance, z is normalized to 1.0
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if
    
    R = R_matrix(i,j)
  end function R
  
  !! normalized inductance L function
  function L(z,i)
    real(dp) L
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance, z normalized to 1.0
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if

    L = L_matrix(i,j)
  end function L
  
  !! normalized conductance G function
  function G(z,i)
    real(dp) G
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance, z normalized to 1.0
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if
    
    G = G_matrix(i,j)
  end function G
  
  !! normalized capacitance Ca function
  function Ca(z,i)
    real(dp) Ca
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if
    
    Ca = C_matrix(i,j)
  end function Ca
  
  !! normalized space charge reduction factor function
  function nscrf(z,i)
    real(dp) nscrf
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance, z normalized to 1.0
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if
    
    nscrf = nscrf_matrix(i,j)
  end function nscrf

  !! normalized cold circuit phase velocity function
  function nvph(z,i)
    real(dp) nvph
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance, z normalized to 1.0
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if
    
    nvph = nvph_matrix(i,j)
  end function nvph
  
  !! pierce parameter C function
  function pC(z,i)
    real(dp) pC
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance, z normalized to 1.0
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if
    
    pC = pC_matrix(i,j)
  end function pC

  !! complex phase factor Z_\ell(z)
  function Z_factor(z,i)
    complex(dp) Z_factor
    real(dp), intent(in) :: z !! distance
    integer, intent(in) :: i !! frequency index
    
    integer j !!distance integer
    
    !!first get integer for distance, z normalized to 1.0
    j = numerical_parameters % num_grid_pts * z
         !* (z / circuit_parameters % circuit_length)
    !!if z beyond length of TWT, set j to last entry
    if (j > numerical_parameters % num_grid_pts) then
       j = numerical_parameters % num_grid_pts
    end if
    
    Z_factor = Z_factor_matrix(i,j)
  end function Z_factor
  
  
  
  !*** compute dispersion parameters from models ***!
  subroutine tape_helix(i,j,tmpvph,tmpK)
    integer, intent(in) :: i, j !! i is section index and j is freq index
    real(dp), intent(inout) :: tmpvph, tmpK
    
  end subroutine tape_helix

  subroutine sheath_helix(i,j,tmpvph,tmpK)
    integer, intent(in) :: i, j !! i is section index and j is freq index
    real(dp), intent(inout) :: tmpvph, tmpK
    
    real(dp) :: a, gamma, cotanpsi, x1, x2
    !! bessel functions from mini_slatec library
    real(dp) :: dbesi0, dbesi1, dbesk0, dbesk1

    !print*, 'computing sheath helix for i=', i, 'j=', j

    a = circuit_parameters % helix_radius(i)
    cotanpsi = 2.0 * pi * a / circuit_parameters % helix_pitch(i)
    !cotanpsi = circuit_parameters % helix_pitch(i) / (2.0 * a)

    !!set initial guesses around k
    x1 = fl(j) * 2.0 * pi * frequency_parameters % base_frequency * &
         sqrt(mu0*eps0)
    x2 = x1*5.0

    call bracket(x1,x2,i,j,cotanpsi)

    gamma = findRoot(x1,x2,i,j,cotanpsi)

    tmpvph = c/cotanpsi * sqrt(dbesi0(gamma*a)*dbesk0(gamma*a)/ &
         (dbesi1(gamma*a)*dbesk1(gamma*a)))

    tmpK = sqrt(mu0/eps0) * cotanpsi / (2.0*pi) * &
         sqrt(dbesi0(gamma*a)*dbesk0(gamma*a) * &
         (dbesi1(gamma*a)*dbesk1(gamma*a)))

    !print*, tmpvph/c, tmpK
  end subroutine sheath_helix

  !! given x1, x2, adjust them until they surround a root
  !! converting Numerical recipes in C zbrac() routine to Java
  subroutine bracket(x1, x2, i, j,cotanpsi)
    real(dp), intent(inout) :: x1, x2
    integer, intent(in) :: i, j
    real(dp), intent(in) :: cotanpsi

    integer :: NTRY = 50
    integer :: kk
    real(dp) :: FACTOR = 1.6
    real(dp) :: f1, f2

    !!print*, 'In Bracket Function'
	
    if (x1 == x2) then
       print*, 'bad initial range in bracket'
    end if

    f1 = dispersionFunction(x1,i,j,cotanpsi)
    f2 = dispersionFunction(x2,i,j,cotanpsi)

    !!print*, 'initial: f(x1) = " + f1 + "\t f(x2) = " + f2'

    do kk = 1, NTRY
       if (f1*f2 < 0.0) then
          return
       else if (abs(f1) < abs(f2)) then
          x1 = x1 + FACTOR*(x1-x2)
          f1 = dispersionFunction(x1,i,j,cotanpsi)
       else
          x2 = x2 + FACTOR*(x2-x1)
          f2 = dispersionFunction(x2,i,j,cotanpsi)
       end if
       !!print*, 'j = ', j, 'x1 = ', x1, ' f(x1) = ', &
       !!& f1, 'x2 = ', x2, ' f(x2) = ', f2
    end do
  end subroutine bracket

  function findRoot(x1, x2, i, j,cotanpsi)
    real(dp) :: findRoot
    real(dp), intent(in) :: x1, x2
    integer, intent(in) :: i, j
    real(dp), intent(in) :: cotanpsi


    integer :: KKMAX = 100
    integer :: kk
    real(dp) :: xacc, dx, f, fmid, xmid, rtb

    findRoot = 0.0

    xacc = 1.0E-15 * (abs(x1) + abs(x2))/2.0

    f = dispersionFunction(x1,i,j,cotanpsi)
    fmid = dispersionFunction(x2,i,j,cotanpsi)

    if (f*fmid >= 0.0) then
       print*, '*****************************************************'
       print*, 'In sheath helix: root must be bracketed for bisection'
       print*, 'Stopping.'
       print*, '*****************************************************'
       stop
    end if

    if (f < 0.0) then
       rtb = x1
       dx = x2 - x1
    else
       rtb = x2
       dx = x1 - x2
    end if
	
    do kk = 1,  KKMAX
       dx = dx * 0.5
       xmid = rtb + dx
       fmid = dispersionFunction(xmid,i,j,cotanpsi)
       if (fmid <= 0.0) then
          rtb = xmid
       end if
       if (abs(dx) < xacc .or. fmid == 0.0) then
          findRoot = rtb
          return
       end if
    end do

    print*,'Too many bisections'
    return

  end function findRoot

  function dispersionFunction(g,i,j,cotanpsi)
    real(dp) :: dispersionFunction
    real(dp), intent(in) :: g, cotanpsi
    integer, intent(in) :: i,j

    real(dp) :: k, a
    !! bessel functions from mini_slatec library
    real(dp) :: dbesi0, dbesi1, dbesk0, dbesk1

    k = fl(j)*2.0*pi* frequency_parameters % base_frequency*sqrt(mu0*eps0)
    a = circuit_parameters % helix_radius(i)

    dispersionFunction = g**2 - k**2 * cotanpsi**2 &
         * dbesi1(g*a) * dbesk1(g*a)/(dbesi0(g*a)*dbesk0(g*a))
  end function dispersionFunction

  subroutine data_helix(i,j,tmpvph,tmpK)
    integer, intent(in) :: i, j !! i is section index and j is freq index
    real(dp), intent(inout) :: tmpvph, tmpK
    
    integer k, tmpLow, tmpHigh
    real(dp) vp2, vp1, K2, K1
    
    !!determine if frequency is a ckt freq
    if (fl(j) < frequency_parameters % min_ckt_freq &
         .or. fl(j) > frequency_parameters % max_ckt_freq) then
       !print "('f(',I3,') =',I3,' outside of ckt_freqs')", j, fl(j)
       tmpvph = 0.0
       tmpK = 0.0
    else
       !!find neighboring frequencies in disp_freq_integer
       k = 1
       tmpHigh = dispersion_parameters % dispersion_freq_integer(i,1)
       !find low-side frequency
       do while (tmpHigh < fl(j) &
            .and. k+1 <= dispersion_parameters % num_dispersion_freqs)
          k = k + 1
          tmpHigh = dispersion_parameters % dispersion_freq_integer(i,k)
       end do
       !!assign tmpLow
       if (k-1 >= 1) then
          tmpLow = dispersion_parameters % dispersion_freq_integer(i,k-1)
       else
          tmpLow = tmpHigh
       end if
       
       !!assign the answer
       if (tmpHigh == fl(j)) then
          tmpvph = dispersion_parameters % phase_velocity(i,k)
          tmpK = dispersion_parameters % impedance(i,k)
       else if (tmpLow < fl(j) .and. fl(j) < tmpHigh) then
          !!interpolate their values
          vp2 = dispersion_parameters % phase_velocity(i,k)
          vp1 = dispersion_parameters % phase_velocity(i,k-1)
          tmpvph = (vp2-vp1)/(tmpHigh-tmpLow) * (fl(j)-tmpLow) + vp1
          K2 = dispersion_parameters % impedance(i,k)
          K1 = dispersion_parameters % impedance(i,k-1)
          tmpK = (K2-K1)/(tmpHigh-tmpLow) * (fl(j)-tmpLow) + K1
       else
          !print*, tmpLow, fl(j), tmpHigh
          print*, 'Range of frequencies in dispersion namelist is not equal to or larger'
          print*, 'than range of computed circuit frequencies.'
          print*, 'Stopping.'
          stop
       end if
       !print*, tmpLow, fl(j), tmpHigh
       !print*, vp1, tmpvph, vp2
       !print*, K1, tmpK, K2
    end if
  end subroutine data_helix

  subroutine antonsen_formula(i,j,tmp_scrf, tmpvph)
    integer, intent(in) :: i, j !! i is section index and j is freq index
    real(dp), intent(inout) :: tmp_scrf
    real(dp), intent(in) :: tmpvph !use this to use beta_c

    !! bessel functions from mini_slatec library
    real(dp) :: dbesi0, dbesi1, dbesk0, dbesk1
    !! local variables
    real(dp) :: kappa
    real(dp) :: rbo, rbi, rh
    real(dp) :: term1, term2, term3

    !! assign local variables
    if (dispersion_parameters % use_beta_c .and. tmpvph /= 0.0) then
       ! use beta_c
       kappa = sqrt((2.0*pi*fl(j)*frequency_parameters % base_frequency)**2 &
            *(1.0/(tmpvph**2) - 1.0/(c**2)))
    else
       ! use beta_e
       kappa = sqrt((2.0*pi*fl(j)*frequency_parameters % base_frequency)**2 &
            *(1.0/(derived_qtys % u0**2) - 1.0/(c**2)))
    end if
    rbo = beam_parameters % outer_radius
    rbi = beam_parameters % inner_radius
    rh = circuit_parameters % helix_radius(i)

    if (rh == 0.0) then
       print*, 'Helix radius must be nonzero to compute space charge reduction factor.'
       print*, 'Stopping.'
       stop
    end if
    
    if (rbi /= 0.0) then !!annular beam
       term1 = rbo*dbesi1(kappa*rbo) - rbi*dbesi1(kappa*rbi)
       term2 = -rbo*dbesk1(kappa*rbo) &
            - (dbesk0(kappa*rh)/dbesi0(kappa*rh))*term1
       term3 = rbi*dbesi1(kappa*rbi)*(rbi*dbesk1(kappa*rbi) &
            - rbo*dbesk1(kappa*rbo))
    else !! pencil beam
       term1 = rbo*dbesi1(kappa*rbo)
       term2 = -rbo*dbesk1(kappa*rbo) &
            - (dbesk0(kappa*rh)/dbesi0(kappa*rh))*term1
       term3 = 0.0
    end if

    !! assign the answer
    tmp_scrf = 1 + 2.0/(rbo**2 - rbi**2) * (term1*term2 - term3)

  end subroutine antonsen_formula

  subroutine data_scrf(i,j,tmp_scrf)
    integer, intent(in) :: i, j !! i is section index and j is freq index
    real(dp), intent(inout) :: tmp_scrf
    
    integer k, tmpLow, tmpHigh
    real(dp) scrf_2, scrf_1
    
    !!determine if frequency is a space_charge_freq
    if (fl(j) < frequency_parameters % min_space_charge_freq &
         .or. fl(j) > frequency_parameters % max_space_charge_freq) then
       !print "('f(',I3,') =',I3,' outside of ckt_freqs')", j, fl(j)
       tmp_scrf = 0.0
    else
       !!find neighboring frequencies in disp_freq_integer
       k = 1
       tmpHigh = dispersion_parameters % dispersion_freq_integer(i,1)
       !find low-side frequency
       do while (tmpHigh < fl(j) &
            .and. k+1 <= dispersion_parameters % num_dispersion_freqs)
          k = k + 1
          tmpHigh = dispersion_parameters % dispersion_freq_integer(i,k)
       end do
       !!assign tmpLow
       if (k-1 >= 1) then
          tmpLow = dispersion_parameters % dispersion_freq_integer(i,k-1)
       else
          tmpLow = tmpHigh
       end if
       
       !!assign the answer
       if (tmpHigh == fl(j)) then
          tmp_scrf = dispersion_parameters % space_charge_redux(i,k)
       else if (tmpLow < fl(j) .and. fl(j) < tmpHigh) then
          !!interpolate their values
          scrf_2 = dispersion_parameters % space_charge_redux(i,k)
          scrf_1 = dispersion_parameters % space_charge_redux(i,k-1)
          tmp_scrf = (scrf_2-scrf_1)/(tmpHigh-tmpLow) * (fl(j)-tmpLow) &
               + scrf_1
       else
          print*, 'Range of frequencies in dispersion namelist is not equal to or larger'
          print*, 'than range of computed space charge frequencies.'
          print*, 'Stopping.'
          stop
       end if
       !print*, tmpLow, fl(j), tmpHigh
       !print*, scrf_1, tmp_scrf, scrf_2
    end if
  end subroutine data_scrf
  
  subroutine loss_model(i,j,tmp_alpha)
    integer, intent(in) :: i, j !! i is section index and j is freq index
    real(dp), intent(inout) :: tmp_alpha
    
    print*, 'loss_model not implemented, using data'
  end subroutine loss_model
  
  subroutine data_alpha(i,j,tmp_alpha)
    integer, intent(in) :: i, j !! i is section index and j is freq index
    real(dp), intent(inout) :: tmp_alpha
    
    integer k, tmpLow, tmpHigh
    real(dp) alpha_2, alpha_1
    
    !!determine if frequency is a ckt freq
    if (fl(j) < frequency_parameters % min_ckt_freq &
         .or. fl(j) > frequency_parameters % max_ckt_freq) then
       !print "('f(',I3,') =',I3,' outside of ckt_freqs')", j, fl(j)
       tmp_alpha = 0.0
    else
       !!find neighboring frequencies in loss_freq_integer
       k = 1
       tmpHigh = loss_parameters % loss_freq_integer(i,1)
       !find low-side frequency
       do while (tmpHigh < fl(j) &
            .and. k+1 <= loss_parameters % number_loss_freqs)
          k = k + 1
          tmpHigh = loss_parameters % loss_freq_integer(i,k)
       end do
       !!assign tmpLow
       if (k-1 >= 1) then
          tmpLow = loss_parameters % loss_freq_integer(i,k-1)
       else
          tmpLow = tmpHigh
       end if
       
       !!assign the answer
       if (tmpHigh == fl(j)) then
          tmp_alpha = loss_parameters % loss(i,k)
       else if (tmpLow < fl(j) .and. fl(j) < tmpHigh) then
          !!interpolate their values
          alpha_2 = loss_parameters % loss(i,k)
          alpha_1 = loss_parameters % loss(i,k-1)
          tmp_alpha = (alpha_2-alpha_1)/(tmpHigh-tmpLow) &
               * (fl(j)-tmpLow) + alpha_1
       else
          print*, 'Range of frequencies in losses namelist is not equal to or larger'
          print*, 'than range of computed circuit frequencies.'
          print*, 'Stopping.'
          stop
       end if
       !print*, tmpLow, fl(j), tmpHigh
       !print*, alpha_1, tmp_alpha, alpha_2
    end if
  end subroutine data_alpha
  
  subroutine interpolate_dispersion(vph_matrix, K_matrix, scrf_matrix)
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: vph_matrix
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: K_matrix
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: scrf_matrix
    
    integer :: i, j, k, arr_start, arr_end, section
    real(dp) :: z_start, z_end, vp1, vp2, K1, K2, scrf_1, scrf_2

    !! intrplt_over_length takes precedence, so do it first
    if (dispersion_parameters % intrplt_over_length) then

       !! loop on circuit sections
       do i = 1, circuit_parameters % number_ckt_sections - 1
          !! get length indexing
          z_start = circuit_parameters % ckt_section_location(i+1) &
               - dispersion_parameters % interpol_length/2.0
          z_end = z_start + dispersion_parameters % interpol_length
          !! check that z_start and z_end are well defined, i.e.
          !! 1/2 interpol_length in an internal section extends no
          !! more than half section length, and is not longer than
          !! the first and last sections
          if (z_start - dispersion_parameters % interpol_length/2.0 &
               < 0.0 .and. i == 1) then
             print*, 'interpolate dispersion:'
             print "(' (1/2 * interpol_length) greater than first section length.')"
             print*, 'Stopping.'
             stop
          else if (z_start - dispersion_parameters % interpol_length/2.0 &
               < circuit_parameters % ckt_section_location(i) &
               .and. i > 1) then
             print*, 'interpolate dispersion:'
             print "(' interpol_length greater than section length for section #',I2)", i
             print*, 'Stopping.'
             stop
          end if
          if (z_end + dispersion_parameters % interpol_length/2.0 &
               > circuit_parameters % ckt_section_location(i+2) &
               .and. i+2 <= circuit_parameters % number_ckt_sections) then
             print*, 'interpolate dispersion:'
             print "(' interpol_length greater than section length for section #',I2)", i+1
             print*, 'Stopping.'
             stop
          else if (z_end > circuit_parameters % circuit_length &
               .and. i+2 > circuit_parameters % number_ckt_sections) then
             print*, 'interpolate dispersion:'
             print "(' (1/2 * interpol_length) greater than final section length, section #',I2)", i+1
             print*, 'Stopping.'
             stop
          end if
       
          !! get indices for these z_start and z_end
          arr_start = numerical_parameters % num_grid_pts &
               * (z_start / circuit_parameters % circuit_length)
          arr_end = numerical_parameters % num_grid_pts &
               * (z_end / circuit_parameters % circuit_length)
       
          !print*, z_start, z_end
          !print*, arr_start, arr_end
       
          do j = 1, M !! loop on all frequencies
             vp1 = vph_matrix(j,arr_start)
             vp2 = vph_matrix(j,arr_end)
             K1 = K_matrix(j,arr_start)
             K2 = K_matrix(j,arr_end)
             scrf_1 = scrf_matrix(j,arr_start)
             scrf_2 = scrf_matrix(j,arr_end)
             do k = arr_start, arr_end
                vph_matrix(j,k) = (vp2-vp1)/(arr_end - arr_start) &
                     * (k - arr_start) + vp1
                vph_matrix(-j,k) = (vp2-vp1)/(arr_end - arr_start) &
                     * (k - arr_start) + vp1
                K_matrix(j,k) = (K2-K1)/(arr_end - arr_start) &
                     * (k - arr_start) + K1
                K_matrix(-j,k) = (K2-K1)/(arr_end - arr_start) &
                     * (k - arr_start) + K1
                scrf_matrix(j,k) = (scrf_2-scrf_1)/(arr_end - arr_start) &
                     * (k - arr_start) + scrf_1
                scrf_matrix(-j,k) = (scrf_2-scrf_1)/(arr_end - arr_start) &
                     * (k - arr_start) + scrf_1
             end do
          end do
       end do
    else if (dispersion_parameters % intrplt_sections) then
       !! first check that interpol_sects_list is in proper order
       do i = 1, circuit_parameters % number_ckt_sections - 1
          section = dispersion_parameters % interpol_sects_list(i)

          if (i == 1 .and. section == 1) then
             print*, ''
             print*, '*******************************************************'
             print*, 'WARNING: 1 is not a valid entry in interpol_sects_list.'
             print*, 'Skipping interpolation.'
             print*, '*******************************************************'
             print*, ''
             return
          end if


          if (section /= -1 .and. &
               section <= circuit_parameters % number_ckt_sections) then
             !! do the interpolation for this section

             !! get length indexing
             z_start = circuit_parameters % ckt_section_location(section)
             if (section <= circuit_parameters % number_ckt_sections - 1) then
                z_end = circuit_parameters % ckt_section_location(section+1)
             else
                z_end = circuit_parameters % circuit_length
             end if
             !print*, z_start, z_end
             !! get indices for z_start and z_end
             arr_start = numerical_parameters % num_grid_pts &
                  * (z_start / circuit_parameters % circuit_length)
             arr_end = numerical_parameters % num_grid_pts &
                  * (z_end / circuit_parameters % circuit_length)
             !print*, arr_start, arr_end

             do j = 1, M !! loop on all frequencies
                !! pick left values from previous section, pick right values
                !! from end of current section.
                vp1 = vph_matrix(j,arr_start-1)
                vp2 = vph_matrix(j,arr_end-1)
                !print*, vp1, vp2
                K1 = K_matrix(j,arr_start-1)
                K2 = K_matrix(j,arr_end-1)
                scrf_1 = scrf_matrix(j,arr_start-1)
                scrf_2 = scrf_matrix(j,arr_end-1)
                do k = arr_start, arr_end
                   vph_matrix(j,k) = (vp2-vp1)/(arr_end - arr_start) &
                        * (k - arr_start) + vp1
                   vph_matrix(-j,k) = (vp2-vp1)/(arr_end - arr_start) &
                        * (k - arr_start) + vp1
                   K_matrix(j,k) = (K2-K1)/(arr_end - arr_start) &
                        * (k - arr_start) + K1
                   K_matrix(-j,k) = (K2-K1)/(arr_end - arr_start) &
                        * (k - arr_start) + K1
                   scrf_matrix(j,k) = (scrf_2-scrf_1)/(arr_end - arr_start) &
                        * (k - arr_start) + scrf_1
                   scrf_matrix(-j,k) = (scrf_2-scrf_1)/(arr_end - arr_start) &
                        * (k - arr_start) + scrf_1
                end do
             end do

          end if

          if (section > circuit_parameters % number_ckt_sections) then
             print*, ''
             print*, '******************************************************&
                  &********************'
             print*, 'WARNING: entry in interpol_sects_list is greater than &
                  &number_ckt_sections.'
             print*, 'Skipping entry.'
             print*, '******************************************************&
                  &********************'
             print*, ''
          end if
       end do
    end if
  end subroutine interpolate_dispersion
  
  subroutine interpolate_loss(alpha_matrix)
    real(dp), dimension(-M:M,0:numerical_parameters % num_grid_pts), &
         intent(inout) :: alpha_matrix
    
    integer i, j, k, arr_start, arr_end
    real(dp) z_start, z_end, alpha_1, alpha_2
    
    !! interpolate over length
    if (loss_parameters % intrplt_over_length) then
       !! loop on loss locations
       do i = 1, loss_parameters % number_loss_locations - 1
          !! get length indexing
          z_start = loss_parameters % loss_location(i+1) &
               - loss_parameters % interpol_length/2.0
          z_end = z_start + loss_parameters % interpol_length
          !! check that z_start and z_end are well defined, i.e.
          !! 1/2 interpol_length in an internal section extends no
          !! more than half section length, and is not longer than
          !! the first and last sections
          if (z_start - loss_parameters % interpol_length/2.0 &
               < 0.0 .and. i == 1) then
             print*, 'interpolate loss:'
             print "(' (1/2 * interpol_length) greater than first section length.')"
             print*, 'Stopping.'
             stop
          else if (z_start - loss_parameters % interpol_length/2.0 &
               < loss_parameters % loss_location(i) &
               .and. i > 1) then
             print*, 'interpolate loss:'
             print "(' interpol_length greater than section length for section #',I2)", i
             print*, 'Stopping.'
             stop
          end if
          if (z_end + loss_parameters % interpol_length/2.0 &
               > loss_parameters % loss_location(i+2) &
               .and. i+2 <= loss_parameters % number_loss_locations) then
             print*, 'interpolate loss:'
             print "(' interpol_length greater than section length for section #',I2)", i+1
             print*, 'Stopping.'
             stop
          else if (z_end > circuit_parameters % circuit_length &
               .and. i+2 > loss_parameters % number_loss_locations) then
             print*, 'interpolate loss:'
             print "(' (1/2 * interpol_length) greater than final section length, section #',I2)", i+1
             print*, 'Stopping.'
             stop
          end if
          
          !! get indices for these z_start and z_end
          arr_start = numerical_parameters % num_grid_pts &
               * (z_start / circuit_parameters % circuit_length)
          arr_end = numerical_parameters % num_grid_pts &
               * (z_end / circuit_parameters % circuit_length)
          
          do j = 1, M !! loop on all frequencies
             alpha_1 = alpha_matrix(j,arr_start)
             alpha_2 = alpha_matrix(j,arr_end)
             do k = arr_start, arr_end
                alpha_matrix(j,k) = (alpha_2-alpha_1)/(arr_end - arr_start) &
                     * (k - arr_start) + alpha_1
                alpha_matrix(-j,k) = (alpha_2-alpha_1)/(arr_end - arr_start)&
                     * (k - arr_start) + alpha_1
             end do
          end do
       end do !! end loop on loss locations
       !! interpolate between points
    else if (loss_parameters % intrplt_btwn_points) then
       if (loss_parameters % loss_location(loss_parameters % &
            number_loss_locations) < circuit_parameters % &
            circuit_length) then
          print*, 'Final loss_location is less than circuit_length.'
          print*, 'Stopping.'
          stop
       end if
       !! loop on loss locations
       do i = 1, loss_parameters % number_loss_locations - 1
          z_start = loss_parameters % loss_location(i)
          z_end = loss_parameters % loss_location(i+1)
          !! get indices for these z_start and z_end
          arr_start = numerical_parameters % num_grid_pts &
               * (z_start / circuit_parameters % circuit_length)
          arr_end = numerical_parameters % num_grid_pts &
               * (z_end / circuit_parameters % circuit_length)
          
          !print*, z_start, z_end
          !print*, arr_start, arr_end
          
          !! loop on all frequencies
          do j = 1, M
             alpha_1 = alpha_matrix(j,arr_start)
             alpha_2 = alpha_matrix(j,arr_end)
             !! loop on points in alpha_matrix
             do k = arr_start, arr_end
                alpha_matrix(j,k) = (alpha_2-alpha_1)/(arr_end - arr_start) &
                     * (k - arr_start) + alpha_1
                alpha_matrix(-j,k) = (alpha_2-alpha_1)/(arr_end - arr_start)&
                     * (k - arr_start) + alpha_1
             end do !! end loop on points in alpha_matrix
          end do !! end loop on frequencies
          
       end do !! end loop on loss locations
    end if
    
  end subroutine interpolate_loss
  
end module working_variables
