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
subroutine latte_dc_constant(n_eqns)
  use precision95
  use parameters, only : pi
  use parameters, only : smo
  use parameters, only : eps0
  use parameters, only : output_parameters
  use parameters, only : run_parameters
  use parameters, only : circuit_parameters
  use parameters, only : units_structure
  use parameters, only : numerical_parameters
  use parameters, only : beam_parameters
  use parameters, only : frequency_parameters
  use working_variables, only : derived_qtys
  use working_variables, only : w0
  use working_variables, only : modelID
  use working_variables, only : run_two
  use working_variables, only : k_1
  use working_variables, only : k_2
  use working_variables, only : k_3
  use working_variables, only : k_4
  use working_variables, only : yk_1
  use working_variables, only : yk_2
  use working_variables, only : yk_3
  use working_variables, only : y_nminus1
  use working_variables, only : k_45
  use working_variables, only : c_45
  use working_variables, only : b_45
  use working_variables, only : bhat_45
  use working_variables, only : a_45
  use working_variables, only : M
  use working_variables, only : fl
  use working_variables, only : first_ckt_index
  use working_variables, only : last_ckt_index
  use working_variables, only : first_spch_index
  use working_variables, only : last_spch_index
  use working_variables, only : vs_z_plot
  use working_variables, only : run_three
  use working_variables, only : N_disks
  use working_variables, only : K
  use working_variables, only : pC
  use working_variables, only : disk_matrix
  use working_variables, only : dcconst
  use working_variables, only : netE
  use working_variables, only : rhol
  use working_variables, only : exp_sum_real
  use working_variables, only : exp_sum_imag
  use working_variables, only : exp_real
  use working_variables, only : exp_imag
  use working_variables, only : thesis_V
  use working_variables, only : thesis_w
  use working_variables, only : ipiv

  integer, intent(in) :: n_eqns


  character(len=30) :: filename
  integer :: i, j, kk, l
  real(dp) :: dz, z, tmpdcvel, tmptmpdcvel, tmpmoddcvel, tmptmpmoddcvel, &
       tmpP, tmpPhi
  complex(dp) :: tmpV
  complex(dp), dimension(:,:), allocatable :: dc_plot_cube
  complex(dp), dimension(:), allocatable :: dc_tmpY, dc_tmpmodY
  complex(dp), dimension(:), allocatable :: y
    
  if (run_two) then
     print*, 'Running latte with dc_constant = true'
     print*, 'on first pass to compute average dc velocities'
  end if

  !! subroutines contained within this subroutine are basically reproductions
  !! of those in commonTWT.i90. This is the only way I could figure out how
  !! to fit into the current structure of the code without cluttering up
  !! commonTWT.i90 with stuff specific to running LATTE with dc_constant = true

  dcconst = .false.

  allocate(y(n_eqns))
  !********************************************************************!
  ! this is a repeat of the routine single_pass_only from commonTWT.i90

  call dc_init_latte_vector(y)
  call dc_norm_latte_vector(y)

  allocate(dc_plot_cube(output_parameters % num_axial_points, n_eqns))

  call dc_singlepass
  !********************************************************************!

  !call dc_mag_vs_z(dc_plot_cube)
  if (output_parameters % disk_orbits_vs_z) then
     call dc_disk_orbits(dc_plot_cube)
  end if

  if (run_parameters % num_scan_namelists /= 0) then
     !! write the voltage phase info to a file
     open(2,file='./outputs/phase_vs_pin.DCC.dat',action='write', &
          position='append')
     !! assuming that the first voltage i=1 is the one we want!
     tmpV = dc_plot_cube(1,1)*K(0.0d0,1)*beam_parameters % current &
          / pC(0.0d0,1)
     tmpP = 10.0*log10(real(2.0*tmpV*conjg(tmpV)/(K(0.0d0,1)))/0.001)
     !! then make tmpV the output voltage
     tmpV = dc_plot_cube(output_parameters % num_axial_points,1) &
          *K(0.0d0,1)*beam_parameters % current / pC(0.0d0,1) !! uniform params
     tmpPhi = 180.0*atan2(aimag(tmpV),real(tmpV))/pi
     write(2,fmt="(en20.8,en20.8)") tmpP, tmpPhi
     close(2)
  end if

  !! disk_matrix has the velocities taken from the plot cube
  if (.not. allocated(disk_matrix)) then
     allocate(disk_matrix(2*N_disks,output_parameters % num_axial_points))
  end if

  do i = 1, 2*N_disks
     do j = 1, output_parameters % num_axial_points
        disk_matrix(i,j) = dc_plot_cube(j,3*M+i)
     end do
  end do
  deallocate(dc_plot_cube)


  !! save dc velocity into a vector and put the reaveraged velocity
  !! into the velocity matrix. export dc velocity and reaveraged
  !! dc velocity into dc_beam_vs_z.L.sp.dat
  write(filename,fmt="(A29)") './outputs/dcbeam_vs_z.DCC.dat'
  open(1,file=filename,action='write')
  !! set up the z step information
  dz = circuit_parameters % circuit_length &
       / (output_parameters % num_axial_points - 1)
  !! multiply scale by z for length in output file
  if (units_structure % length_cm) then
     scale = 100.0
  else
     scale = 1.0
  end if
  z = 0.0

  allocate(dc_tmpY(n_eqns))
  do i = 1, 3*M
     dc_tmpY(i) = 0.0
  end do

  do i = 1, output_parameters % num_axial_points
     do j = 1, 2*N_disks
        dc_tmpY(3*M+j) = disk_matrix(j,i)
     end do
     call dc_unnorm_latte_vector(dc_tmpY)

     !! first compute average
     print*, 'computing dc: i=', i, 'of', output_parameters % num_axial_points
     tmpdcvel = 0.0
     do j = 1, N_disks
        tmptmpdcvel = 0.0
        !! first do the l == 0 contribution
        do kk = 1, N_disks
           tmptmpdcvel = tmptmpdcvel &
                + 1.0/real(dc_tmpY(3*M+kk))
        end do
        !! then do the l /= 0 contributions
        do l = first_spch_index, last_spch_index
           do kk = 1, N_disks
              tmptmpdcvel = tmptmpdcvel &
                   + 2.0*cos(fl(l)*(real(dc_tmpY(3*M+N_disks+j)) &
                   - real(dc_tmpY(3*M+N_disks+kk))))/real(dc_tmpY(3*M+kk))
           end do
        end do
        tmpdcvel = tmpdcvel + 1.0/tmptmpdcvel
     end do
     
     !! subtract average from each velocity, add u0 back on,
     !! then recompute average, do we get nearly u0?
     allocate(dc_tmpmodY(n_eqns))
     do j = 1, 3*M
        dc_tmpmodY(j) = 0.0
     end do
     !! subtract off the variable dc component, add the nonvariable
     do j = 1, 2*N_disks
        dc_tmpmodY(3*M+j) = dc_tmpY(3*M+j) - tmpdcvel + derived_qtys % u0
     end do
     
     !! redo the averaging
     tmpmoddcvel = 0.0
     do j = 1, N_disks
        tmptmpmoddcvel = 0.0
        !! first do the l == 0 contribution
        do kk = 1, N_disks
           tmptmpmoddcvel = tmptmpmoddcvel &
                + 1.0/real(dc_tmpmodY(3*M+kk))
        end do
        !! then do the l /= 0 contributions
        do l = first_spch_index, last_spch_index
           do kk = 1, N_disks
              tmptmpmoddcvel = tmptmpmoddcvel &
                   + 2.0*cos(fl(l)*(real(dc_tmpmodY(3*M+N_disks+j)) &
                   - real(dc_tmpmodY(3*M+N_disks+kk)))) &
                   /real(dc_tmpmodY(3*M+kk))
           end do
        end do
        tmpmoddcvel = tmpmoddcvel + 1.0/tmptmpmoddcvel
     end do
     deallocate(dc_tmpmodY)

     !! modify the disk_matrix values
     do j = 1, N_disks
        disk_matrix(j,i) = disk_matrix(j,i) - tmpdcvel/derived_qtys % u0 + 1.0
     end do

     if (units_structure % length_cm) then
        tmpdcvel = 100.0 * tmpdcvel
        tmpmoddcvel = 100.0 * tmpmoddcvel
     end if
     write(1, fmt="(en15.5, en20.8, en20.8)") z*scale, tmpdcvel, tmpmoddcvel

     z = z+dz
  end do

  !do i = 1, N_disks
  !do j = 1, output_parameters % num_axial_points
  !print*, disk_matrix(i,j)
  !end do
  !end do

  close(1)
  deallocate(dc_tmpY)
  deallocate(y)

  !! now the usual LATTE will run with dcconst = .true.
  dcconst = .true.

  if (run_two) then
     print*, 'finished first pass'
     print*, 'on second pass to correct for dc velocities'
  end if
  
contains

  !*** initialize latte vector ***!
  subroutine dc_init_latte_vector(dc_y)
    complex(dp), dimension(:), intent(inout) :: dc_y
    integer :: i, j

    !initialize y vector first to zeros
    do i = 1, n_eqns
       dc_y(i) = cmplx(0.0,0.0)
    end do

    !initialize circuit voltage
    j = 1
    do i = 1, M !!loop on frequencies, pick out drives and assign voltage
       if (fl(i) == frequency_parameters % frequency_integer(j)) then
          !positive frequency
          dc_y(i) = &
               sqrt(frequency_parameters % power_input(j) * K(0.0d0,i)/2.0) &
               * exp(smo * frequency_parameters % phase_input(j))
          !print*, log10(abs(dc_y(i)))
          j = j + 1
       end if
    end do

    !initialize circuit current
    !!loop on ckt frequencies
    !do i = first_ckt_index, last_ckt_index
    do i = first_ckt_index, last_ckt_index
       dc_y(i+M) = -dc_y(i)/K(0.0d0,i)
       !print*, 'i=',i,dc_y(i+M)
    end do
    
    !initialize disk velocities
    do i = 3*M + 1, 3*M + N_disks
       dc_y(i) = cmplx(derived_qtys % u0)
       !print*, 'i=',i,dc_y(i)
    end do
 
    !initialize disk phases
    j = 0
    do i = 3*M + N_disks + 1, 3*M + 2*N_disks
       !dc_y(i) = cmplx((2.0*pi)*j/N_disks)
       dc_y(i) = -cmplx((2.0*pi)*j/N_disks)
       !print*, 'i=',i,dc_y(i)
       j = j + 1
    end do
  end subroutine dc_init_latte_vector

  !*** normalize latte vector ***!
  subroutine dc_norm_latte_vector(dc_y)
    complex(dp), dimension(:), intent(inout) :: dc_y
    integer :: i,j

    !normalize voltage
    do i = 1, M !!loop on frequencies, pick out drives and assign voltage
       !positive frequency
       if (i >= first_ckt_index .and. i <= last_ckt_index) then
  	  dc_y(i) = dc_y(i) &
               * pC(0.0d0,i)/(K(0.0d0,i)*beam_parameters % current)
          !print*, log10(abs(dc_y(i)))
       end if
       !print*, 'i=',i,dc_y(i)
    end do
    
    !normalize circuit current
    do i = 1, M !!loop on frequencies
       if (i >= first_ckt_index .and. i <= last_ckt_index) then
	  dc_y(i+M) = dc_y(i+M) * pC(0.0d0,i)/ beam_parameters % current
       end if
       !print*, 'i=',i,dc_y(i+M)
    end do

    !assuming normalize always called with no space charge field?

    !normalize disk velocities
    do i = 3*M + 1, 3*M + N_disks
       dc_y(i) = dc_y(i)/derived_qtys % u0
       !print*, 'i=',i,dc_y(i)
    end do
 
    !normalize disk phases
    do i = 3*M + N_disks + 1, 3*M + 2*N_disks
       dc_y(i) = dc_y(i)/w0
       !print*, 'i=',i,dc_y(i)
    end do
  end subroutine dc_norm_latte_vector

  !*** unnormalize latte vector ***!
  subroutine dc_unnorm_latte_vector(dc_y)
    complex(dp), dimension(:), intent(inout) :: dc_y
    integer :: i

    !unnormalize voltage
    do i = 1, M !!loop on frequencies, pick out drives and assign voltage
       !positive frequency
       if (i >= first_ckt_index .and. i <= last_ckt_index) then
	  dc_y(i) = dc_y(i)*K(0.0d0,i)*beam_parameters % current / pC(0.0d0,i)
       end if
    end do
    
    !unnormalize circuit current
    do i = 1, M !!loop on frequencies
       if (i >= first_ckt_index .and. i <= last_ckt_index) then
          dc_y(i+M) = dc_y(i+M) * beam_parameters % current / pC(0.0d0,i)
       end if
    end do

    !unnormalize space charge field
    do i = 1, M !!loop on frequencies
       if (i >= first_spch_index .and. i <= last_spch_index) then
	  dc_y(i+2*M) = dc_y(i+2*M) * circuit_parameters % circuit_length &
               * derived_qtys % rho0 / eps0
       end if
    end do

    !unnormalize disk velocities
    do i = 3*M + 1, 3*M + N_disks
       dc_y(i) = dc_y(i) * derived_qtys % u0
    end do
    
    !unnormalize disk phases
    do i = 3*M + N_disks + 1, 3*M + 2*N_disks
       dc_y(i) = dc_y(i) * w0
    end do
  end subroutine dc_unnorm_latte_vector

  subroutine dc_mag_vs_z(dc_plot_cube)
    complex(dp), dimension(:,:), intent(in) :: dc_plot_cube
    character(len=40) :: filename
    integer :: i, j
    real(dp) :: z, dz, scale
    complex(dp), dimension(:), allocatable :: dc_tmpY

    if (run_two) then
       print*, '  generating mag_vs_z.DCC.dat file'
    end if

    !! create the file
    write(filename,fmt="(A26)") './outputs/mag_vs_z.DCC.dat'
    open(1,file=filename,action='write')

    !! set up the z step information
    dz = circuit_parameters % circuit_length &
         / (output_parameters % num_axial_points - 1)

    !! multiply scale by z for length in output file
    if (units_structure % length_cm) then
       scale = 100.0
    else
       scale = 1.0
    end if

    !! allocate dc_tmpY, the vector passed to the unnormalization routine
    allocate(dc_tmpY(n_eqns))

    z = 0.0

    !! outer loop. makes one row per axial point
    do i = 1, output_parameters % num_axial_points
       !! copy plot_cube(i,:) into dc_tmpY
       dc_tmpY = dc_plot_cube(i,:)
       call dc_unnorm_latte_vector(dc_tmpY)

       write(1, fmt="(en15.5)", advance='no') z*scale
       do j = first_ckt_index, last_ckt_index
          if (abs(dc_tmpY(j)) == 0.0) then
             dc_tmpY(j) = 1.0e-15
          end if
          write(1, fmt="(en15.5)", advance='no') log10(abs(dc_tmpY(j)))
       end do
       write(1, fmt="(A1)")
       z = z+dz
    end do

    close(1)

  end subroutine dc_mag_vs_z

  subroutine dc_disk_orbits(dc_plot_cube)
    complex(dp), dimension(:,:), intent(in) :: dc_plot_cube

    character(len=40) :: filename
    integer :: i, j
    real(dp) :: z, dz, scale
    complex(dp), dimension(:), allocatable :: dc_tmpY

    if (run_two) then
       print*, '  generating orbits_vs_z.DCC.dat file'
       print*, '  this file is WITHOUT dc correction'
    end if

    !! create the file
    write(filename,fmt="(A29)") './outputs/orbits_vs_z.DCC.dat'
    open(1,file=filename,action='write')

    !! set up the z step information
    dz = circuit_parameters % circuit_length &
         / (output_parameters % num_axial_points - 1)

    !! multiply scale by z for length in output file
    if (units_structure % length_cm) then
       scale = 100.0
    else
       scale = 1.0
    end if

    !! allocate dc_tmpY, the vector passed to the unnormalization routine
    allocate(dc_tmpY(n_eqns))

    z = 0.0

    !! outer loop. makes one row per axial point
    do i = 1, output_parameters % num_axial_points
       !! copy plot_cube(i,:) into dc_tmpY
       dc_tmpY = dc_plot_cube(i,:)
       call dc_unnorm_latte_vector(dc_tmpY)

       write(1, fmt="(en15.5)", advance='no') z*scale
       do j = 1, N_disks
          write(1, fmt="(en15.5)", advance='no') real(dc_tmpY(3*M+N_disks+j))
       end do
       write(1, fmt="(A1)")
       z = z+dz
    end do

    close(1)

  end subroutine dc_disk_orbits


  subroutine dc_singlepass
    implicit none

    external runge_kutta, & !runge kutta integrator
         !F, & !compute right hand side of ODE system
         JACOBN, & !compute jacobian matrix
         FA, & !compute A in A(y,t)*dy/dt = F(y,t)
         G, & !root finder, see cdriv3.f
         USERS !for special linear algebra routines
    
    real(dp) z, & !indep var, 1st call = 0.0, returns z at which sol'n given
         zout, &!point at which soln desired
         eps, & !requested relative acc., set=numerical_parameters % tolerance
         h, & !step size
         hmin, & !minimum magnitude of step-size.
         hmax !maximum magnitude of step-size.

    !! local variables
    integer i
    real(dp) dz, delta, scale

    !! values below are common to all routines
    eps = numerical_parameters % tolerance

    !! normalized length is 1.0, integrate up to this value
    dz = 1.0 / (output_parameters % num_axial_points - 1)
    !! compute delta, a small number to add to 1.0 for checking
    !! in integration kernel. it is the difference between 1 and
    !! adding dz, num_axial_points - 1 times. (not equal to
    !! dz*[num_axial_points-1] !!)
    z = 0.0
    do i = 1, output_parameters % num_axial_points - 1
       z = z + dz
    end do
    !delta = z - 1.0
    delta = z - 0.99999

    !! set initial step size
    h = 1.0 / numerical_parameters % num_grid_pts
    z = 0.0
    zout = dz

    !! allocate vectors used in RK, move this to a subroutine
    allocate(k_1(n_eqns))
    allocate(k_2(n_eqns))
    allocate(k_3(n_eqns))
    allocate(k_4(n_eqns))
    allocate(yk_1(n_eqns))
    allocate(yk_2(n_eqns))
    allocate(yk_3(n_eqns))
    allocate(y_nminus1(n_eqns))
    allocate(k_45(7,n_eqns))
    c_45(1) = 0.0
    c_45(2) = 1.0/5.0
    c_45(3) = 3.0/10.0
    c_45(4) = 7.0/10.0
    c_45(5) = 5.0/6.0
    c_45(6) = 1.0
    c_45(7) = 1.0
    b_45(1) = 59.0/630.0
    b_45(2) = 0.0
    b_45(3) = 125.0/288.0
    b_45(4) = 125.0/504.0
    b_45(5) = 27.0/160.0
    b_45(6) = 1.0/18.0
    b_45(7) = 0.0
    bhat_45(1) = 59.0/630.0 + 4.0/225.0 * 0.1
    bhat_45(2) = 0.0
    bhat_45(3) = -5.0/72.0 * 0.1 + 125.0/288.0
    bhat_45(4) = 5.0/18.0  * 0.1 + 125.0/504.0
    bhat_45(5) = -63.0/200.0 * 0.1 + 27.0/160.0
    bhat_45(6) = -41.0/45.0 * 0.1 + 1.0/18.0
    bhat_45(7) = 0.1
    a_45(2,1) = 1.0/5.0
    a_45(3,1) = 3.0/40.0
    a_45(4,1) = 203.0/360.0
    a_45(5,1) = -8285.0/13608.0
    a_45(6,1) = 136.0/315.0
    a_45(7,1) = 59.0/630.0
    a_45(3,2) = 9.0/40.0
    a_45(4,2) = -49.0/24.0
    a_45(5,2) = 1925.0/648.0
    a_45(6,2) = -5.0/3.0
    a_45(7,2) = 0.0
    a_45(4,3) = 98.0/45.0
    a_45(5,3) = -500.0/243.0
    a_45(6,3) = 575.0/288.0
    a_45(7,3) = 125.0/288.0
    a_45(5,4) = 100.0/189.0
    a_45(6,4) = -15.0/56.0
    a_45(7,4) = 125.0/504.0
    a_45(6,5) = 81.0/160.0
    a_45(7,5) = 27.0/160.0
    a_45(7,6) = 1.0/18.0
    allocate(netE(M))
    allocate(rhol(M))
    allocate(exp_sum_real(M))
    allocate(exp_sum_imag(M))
    allocate(exp_real(M,N_disks))
    allocate(exp_imag(M,N_disks))
    allocate(thesis_V(2*(2*M+1),2*(2*M+1)))
    allocate(thesis_w(2*(2*M+1)))
    allocate(ipiv(2*(2*M+1)))

    !! put the input vector into the plot_cube
    dc_plot_cube(1,:) = y

    i = 2
    !! integration kernel
    do while (i <= output_parameters % num_axial_points &
         .and. zout <= 1.0 + delta)
       if (run_three) then
          if (units_structure % length_cm) then
             scale = 100.0
          else
             scale = 1.0
          end if
          print*, 'in integration kernel. calling runge_kutta() at z = ', &
               zout * circuit_parameters % circuit_length * scale
       end if
       call runge_kutta(n_eqns,z,h,y,zout,eps,hmin,hmax)

       !! put the vector into the plot cube. if .not. vs_z_plot, then i
       !! never goes past 2
       dc_plot_cube(i,:) = y

       zout = zout + dz !advance the output point
       i = i + 1 !advance i
    end do

    !! release working vectors, move this to a subroutine
    deallocate(k_1)
    deallocate(k_2)
    deallocate(k_3)
    deallocate(k_4)
    deallocate(yk_1)
    deallocate(yk_2)
    deallocate(yk_3)
    deallocate(y_nminus1)
    deallocate(k_45)
    deallocate(netE)
    deallocate(rhol)
    deallocate(exp_sum_real)
    deallocate(exp_sum_imag)
    deallocate(exp_real)
    deallocate(exp_imag)
    deallocate(thesis_V)
    deallocate(thesis_w)
    deallocate(ipiv)

  end subroutine dc_singlepass
  
end subroutine latte_dc_constant
