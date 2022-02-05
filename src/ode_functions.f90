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


subroutine runge_kutta(n_eqns,z,h,y,zout,eps,hmin,hmax)
  use precision95
  use working_variables, only : y_nminus1
  integer, intent(in) :: n_eqns
  real(dp), intent(inout) :: z
  real(dp), intent(inout) :: h
  complex(dp), dimension(1:n_eqns), intent(inout) :: y
  real(dp), intent(in) :: zout
  real(dp), intent(in) :: eps
  real(dp), intent(in) :: hmin, hmax

  !************ fixed step routine ***************!

  do while (z <= zout)
     !! compute the next step
     !! answer is returned in y
     y_nminus1 = y ! this holds the last value
     call ynplus1_rk4(n_eqns,z,h,y)
     !call ynplus1_rk45(n_eqns,z,h,y)
     z = z + h
  end do

  !! interpolate back to zout
  do i = 1, n_eqns
     y(i) = (y(i)-y_nminus1(i))*(zout - z + h)/h + y_nminus1(i)
  end do
  !! set z equal to zout
  z = zout


  !************ end fixed step routine ***************!
  
  !! adaptive step size

contains
  subroutine ynplus1_rk4(n_eqns,z,h,y)
    use precision95
    use working_variables, only : k_1
    use working_variables, only : k_2
    use working_variables, only : k_3
    use working_variables, only : k_4
    use working_variables, only : yk_1
    use working_variables, only : yk_2
    use working_variables, only : yk_3

    integer, intent(in) :: n_eqns
    real(dp), intent(inout) :: z
    real(dp), intent(in) :: h
    complex(dp), dimension(1:n_eqns), intent(inout) :: y

    external F

    ! local variables
    integer i, j

    !!first call to F, k_1 is ydot()
    call F(n_eqns,z,y,k_1)
    !!scale k_1 by h
    do i = 1, n_eqns
       k_1(i) = h*k_1(i)
    end do
    
    !!second call to F, k_2 is ydot
    do i = 1, n_eqns
       yk_1(i) = y(i) + 0.5*k_1(i)
    end do
    call F(n_eqns,z+h/2,yk_1,k_2)
    !!scale k_2 by h
    do i = 1, n_eqns
       k_2(i) = h*k_2(i)
    end do


    !!third call to F, k_3 is ydot
    do i = 1, n_eqns
       yk_2(i) = y(i) + 0.5*k_2(i)
    end do
    call F(n_eqns,z+h/2,yk_2,k_3)
    !!scale k_3 by h
    do i = 1, n_eqns
       k_3(i) = h*k_3(i)
    end do
    
    !!fourth call to F, k_4 is ydot
    do i = 1, n_eqns
       !print*, y(i), k_3(i), y(i)+k_3(i)
       !print*, y(i), k_3(i), y(i)+1.0*k_3(i)
       yk_3(i) = y(i) + 1.0*k_3(i)
    end do
    call F(n_eqns,z+h,yk_3,k_4)
    !!scale k_4 by h
    do i = 1, n_eqns
       k_4(i) = h*k_4(i)
    end do
    
    !!assign the answer
    do i = 1, n_eqns
       y(i) = y(i) + k_1(i)/6.0 + k_2(i)/3.0 + k_3(i)/3.0 + k_4(i)/6.0
    end do
    
  end subroutine ynplus1_rk4

  subroutine ynplus1_rk45(n_eqns,z,h,y)
    use precision95
    use working_variables, only : k_45
    use working_variables, only : c_45
    use working_variables, only : b_45
    use working_variables, only : bhat_45
    use working_variables, only : a_45

    integer, intent(in) :: n_eqns
    real(dp), intent(inout) :: z
    real(dp), intent(in) :: h
    complex(dp), dimension(1:n_eqns), intent(inout) :: y

    ! local variables
    integer i, j, k
    complex(dp), dimension(1:n_eqns) :: F_arg

    external F

    call F(n_eqns,z+h*c_45(1),y,k_45(1,:))
    
    do i = 2, 7
       !! make the argument to send to F
       do j = 1, n_eqns
          F_arg(j) = 0.0
          do k = 1, i-1
             F_arg(j) = F_arg(j) + a_45(i,k)*k_45(k,j)
          end do
          F_arg(j) = y(j) + h*F_arg(j)
       end do
       !! call F
       call F(n_eqns,z+h*c_45(i),F_arg,k_45(i,:))
    end do

    !!make the solution
    do i = 1, n_eqns
       F_arg(i) = 0.0
       do j = 1, 7
          F_arg(i) = F_arg(i) + bhat_45(j)*k_45(j,i)
       end do
       y(i) = y(i) + h*F_arg(i)
    end do

  end subroutine ynplus1_rk45

end subroutine runge_kutta

subroutine F(n, z, y, ydot)
  use precision95
  use parameters, only : c
  use parameters, only : eps0
  use parameters, only : smo
  use parameters, only : run_parameters
  use parameters, only : output_parameters
  use parameters, only : circuit_parameters
  use working_variables, only : modelID
  use working_variables, only : M
  use working_variables, only : fl
  use working_variables, only : pair_matrix
  use working_variables, only : find_freq
  use working_variables, only : first_ckt_index
  use working_variables, only : last_ckt_index
  use working_variables, only : first_spch_index
  use working_variables, only : last_spch_index
  use working_variables, only : num_ckt_freqs
  use working_variables, only : num_spch_freqs
  use working_variables, only : w0
  use working_variables, only : N_disks
  use working_variables, only : K
  use working_variables, only : vph
  use working_variables, only : R
  use working_variables, only : L
  use working_variables, only : G
  use working_variables, only : Ca
  use working_variables, only : nscrf
  use working_variables, only : nvph
  use working_variables, only : pC
  use working_variables, only : Z_factor
  use working_variables, only : netE
  use working_variables, only : rhol
  use working_variables, only : exp_sum_real
  use working_variables, only : exp_sum_imag
  use working_variables, only : exp_real
  use working_variables, only : exp_imag
  use working_variables, only : thesis_V
  use working_variables, only : thesis_w
  use working_variables, only : ipiv
  use working_variables, only : derived_qtys
  use working_variables, only : disk_matrix
  use working_variables, only : dcconst

  implicit none

  integer, intent(in) :: n
  real(dp), intent(in) :: z
  complex(dp), dimension(1:n), intent(in) :: y
  complex(dp), dimension(1:n), intent(inout) :: ydot

  integer :: i, j, kk, p, q
  !! variables for LATTE
  complex(dp) :: i_freq, exp_sum, vel_rhs2
  real(dp) :: tmp_arg, tmp_phase, vel_rhs, tmp_oneoverv, tmp_real, tmp_imag
  real(dp) :: dz, dv
  real(dp), dimension(:), allocatable :: vel_vector

  !! variables for MUSE
  integer :: info
  !! variables for S-MUSE
  complex(dp) :: vv, Irho, Erho, vrho

  !print*, dcconst

  !************  for LATTE dc_constant = true  ****************!
  !! interpolate the velocities in disk_matrix
  if (dcconst) then
     if (.not. allocated(vel_vector)) then
        allocate(vel_vector(N_disks))
     end if
     dz = 1.0 / (output_parameters % num_axial_points - 1)
     i = z / dz + 1
     if (i >= output_parameters % num_axial_points) then
        i = output_parameters % num_axial_points - 1
     end if
     do j = 1, N_disks
        dv = disk_matrix(j,i+1) - disk_matrix(j,i)
        vel_vector(j) = (dv/dz)*(z - (i-1)*dz) + disk_matrix(j,i)
        !print*, 'j=',j,vel_vector(j)
     end do
  end if

  !************************************************************!


  !************************************************************!
  !************************************************************!
  !********************* LATTE derivatives ********************!
  !************************************************************!
  !************************************************************!

  !************************************************************!
  !****   LATTE: NO slowly varying envelope approximation  ****!
  !************************************************************!
  if (modelID=='L' .and. .not. run_parameters % svea) then
     !compute 1/v for later computations and store it in ydot(i)
     do i = 3*M + N_disks + 1, 3*M + 2*N_disks
        if (dcconst) then
           !! for dcconst, only have to modify here. whenever v is needed
           !! it takes it from this vector. but 
           ydot(i) = 1.0/vel_vector(i-3*M-N_disks)
        else
           ydot(i) = 1.0/real(y(i-N_disks))
        end if
     end do

     !! make netE vector that is used on RHS of all velocity equations
     do i = 1, M
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           netE(i) = 2.0*pC(z,i)**2*(smo*fl(i)*w0*L(z,i)-R(z,i))*y(i+M) &
                + nscrf(z,i)*y(i+2*M)
        else if (i >= first_spch_index .and. i <= last_spch_index) then
           netE(i) = nscrf(z,i)*y(i+2*M)
        else
           netE(i) = (0.0,0.0);
        end if
        exp_sum_real(i) = 0.0d0
        exp_sum_imag(i) = 0.0d0
     end do

     !! velocity derivatives
     do i = 3*M + 1, 3*M + N_disks
        vel_rhs = 0.0
        tmp_phase = real(y(i+N_disks))*w0
        tmp_oneoverv = real(ydot(i+N_disks))
        do j = first_spch_index, last_spch_index
           tmp_arg = fl(j)*tmp_phase
           tmp_real = tmp_oneoverv*cos(tmp_arg)
           tmp_imag = tmp_oneoverv*sin(tmp_arg)

           !compute velocity RHS
           vel_rhs = vel_rhs + real(netE(j))*tmp_real &
                - aimag(netE(j))*tmp_imag

           !use tmp_real/imag to build exp_sum_real/imag, remember to conjugate
           !and divide by N_disks when using!!
           !note looping on j, i.e. for fixed particle, doing all frequencies
           exp_sum_real(j) = exp_sum_real(j) + tmp_real
           exp_sum_imag(j) = exp_sum_imag(j) + tmp_imag
        end do
        ydot(i) = 2.0*vel_rhs

     end do

     !! phase derivatives
     do i = 3*M + N_disks + 1, 3*M + 2*N_disks
        !\Psi derivative
        ydot(i) = 1 - ydot(i)
     end do

     !circuit and space charge derivatives
     do i = 1, M
        i_freq = smo*fl(i)*w0
        exp_sum = cmplx(exp_sum_real(i),-exp_sum_imag(i))/real(N_disks)
        !! circuit derivatives
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           !! voltage derivatives
           ydot(i) = -i_freq*y(i) + (R(z,i)-i_freq*L(z,i))*y(i+M)

           !print*, z,'i=',i,abs(ydot(i))
           

           !! current derivatives
           ydot(i+M) = (G(z,i) - i_freq*Ca(z,i))*y(i) &
                - i_freq*y(i+M) + i_freq*pC(z,i) * exp_sum
        else
           ydot(i) = (0.0,0.0)
           ydot(i+M) = (0.0,0.0)
        end if

        !! space charge field derivative
        if (i >= first_spch_index .and. i <= last_spch_index) then
           !! space charge derivatives
           ydot(i+2*M) = -i_freq*y(i+2*M) + exp_sum
        else
           ydot(i+2*M) = (0.0,0.0)
        end if
     end do

  !**************************************************************!
  !****   LATTE: WITH slowly varying envelope approximation  ****!
  !**************************************************************!
  else if (modelID=='L' .and. run_parameters % svea) then
     !compute 1/v for later computations and store it in ydot(i)
     do i = 3*M + N_disks + 1, 3*M + 2*N_disks
        if (dcconst) then
           !! for dcconst, only have to modify here. whenever v is needed
           !! it takes it from this vector
           ydot(i) = 1.0/vel_vector(i-3*M-N_disks)
        else
           ydot(i) = 1.0/real(y(i-N_disks))
        end if
     end do

     do i = 1, M
        netE(i) = cmplx(0.0d0,0.0d0)
        rhol(i) = cmplx(0.0d0,0.0d0)
     end do

     !! compute complex exponentials once
     do i = 3*M + 1, 3*M + N_disks
        tmp_phase = w0*real(y(i+N_disks))
        tmp_oneoverv = real(ydot(i+N_disks))
        do j = first_spch_index, last_spch_index
           tmp_arg = fl(j)*tmp_phase
           tmp_real = cos(tmp_arg)
           tmp_imag = sin(tmp_arg)
           !! exp_real and exp_imag contain e^{i\fl\omega_0\Psi}
           exp_real(j,i-3*M) = tmp_real
           exp_imag(j,i-3*M) = tmp_imag
           rhol(j) = rhol(j) + tmp_oneoverv*cmplx(tmp_real,-tmp_imag)
        end do
     end do

     !! make rhol vector containing charge density Fourier coeff,
     !! Il vector for circuit frequencies, netE vector to appear
     !! in velocity equation
     do i = 1, M
        rhol(i) = (1.0/real(N_disks)) * rhol(i)
        i_freq = smo*fl(i)*w0
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           !! voltage derivative
           !ydot(i) = 0.5*(i_freq*(L(z,i)*Ca(z,i) - 1.0) &
           !- smo*R(z,i)*G(z,i)/(fl(i)*w0) &
           !- (R(z,i)*Ca(z,i) + G(z,i)*L(z,i)))*y(i) &
           !+ 0.5*(R(z,i)-i_freq*L(z,i))*pC(z,i)*rhol(i)
           ydot(i) = 0.5*nvph(z,i)*(i_freq*(L(z,i)*Ca(z,i) &
                - 1.0/(nvph(z,i)**2)) &
                - smo*R(z,i)*G(z,i)/(fl(i)*w0) &
                - (R(z,i)*Ca(z,i) + G(z,i)*L(z,i)))*y(i) &
                + 0.5*nvph(z,i)*(R(z,i) &
                -i_freq*L(z,i))*pC(z,i)*Z_factor(z,i)*rhol(i)
           !! current derivatives, no current computed with SVEA
           ydot(i+M) = (0.0,0.0)

           !! netE() without V derivative term
           !netE(i) = -2.0*pC(z,i)**2*i_freq*y(i) + nscrf(z,i)*y(i+2*M)
           netE(i) = -2.0*pC(z,i)**2*i_freq*conjg(Z_factor(z,i))*y(i) &
                /nvph(z,i) + nscrf(z,i)*y(i+2*M)
           !! netE() with V derivative term
           !netE(i) = -2.0*pC(z,i)**2*(ydot(i)+i_freq*y(i))+nscrf(z,i)*y(i+2*M)
        else if (i >= first_spch_index .and. i <= last_spch_index) then
           netE(i) = nscrf(z,i)*y(i+2*M)
        else
           ydot(i) = (0.0,0.0);
           ydot(i+M) = (0.0,0.0);

           netE(i) = (0.0,0.0);
        end if

        if (i >= first_spch_index .and. i <= last_spch_index) then
           !! space charge derivative
           ydot(i+2*M) = -i_freq*y(i+2*M) + rhol(i)
        else
           ydot(i+2*M) = (0.0,0.0);
        end if
     end do

     !! velocity derivatives
     do i = 3*M + 1, 3*M + N_disks
        vel_rhs = 0.0d0
        do j = first_spch_index, last_spch_index
           !compute velocity RHS. exp_real, exp_imag have 1/v in them
           vel_rhs = vel_rhs + real(netE(j))*exp_real(j,i-3*M) &
                - aimag(netE(j))*exp_imag(j,i-3*M)
        end do
        ydot(i) = ydot(i+N_disks)*2.0*vel_rhs
     end do

     !! phase derivatives
     do i = 3*M + N_disks + 1, 3*M + 2*N_disks
        !\Psi derivative
        ydot(i) = 1 - ydot(i)
     end do

  !***********************************************************!
  !***********************************************************!
  !********************* MUSE derivatives ********************!
  !***********************************************************!
  !***********************************************************!
  else if (modelID=='M') then
     !! put DC voltage, current and space charge field to zero
     ydot(M+1) = (0.0,0.0)
     ydot(2*(M+1)) = (0.0,0.0)
     ydot(3*(M+1)) = (0.0,0.0)

     !circuit and space charge derivatives
     do i = 1, M
        !! circuit derivatives
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           i_freq = smo*fl(i)*w0
           !! voltage derivatives
           ydot(i) = -i_freq*y(i) + (R(z,i)-i_freq*L(z,i))*y(i+M+1)

           !! current derivatives
           ydot(i+M+1) = (G(z,i) - i_freq*Ca(z,i))*y(i) &
                - i_freq*y(i+M+1) + i_freq*pC(z,i)*y(i+4*(M+1))
        else
           ydot(i) = (0.0,0.0)
           ydot(i+M+1) = (0.0,0.0)
        end if

        !! space charge field derivative
        if (i >= first_spch_index .and. i <= last_spch_index) then
           !! space charge derivatives
           ydot(i+2*(M+1)) = -i_freq*y(i+2*(M+1)) + y(i+4*(M+1))
        else
           ydot(i+2*(M+1)) = (0.0,0.0)
        end if
     end do

     !! make RHS vector for velocity and density
     do i = 1, M !! loop on frequency \ell
        vv = (0.0,0.0)
        vrho = (0.0,0.0)

        !pm comments starting with !pm were put in when inserting pair_matrix

        !pm kk = 0 !! set difference frequency to zero
        !pm do j = -last_spch_index, last_spch_index ! loop on "m" for products
        
        do p = 1, pair_matrix(i,2)  !pm
           j = pair_matrix(i,2*p+1) !pm
           kk= pair_matrix(i,2*p+2) !pm

           !!f_m /= 0 and f_n /= 0
           if (j /= i .and. j /= 0 .and. abs(j) >= first_spch_index) then

              !pm kk = find_freq(fl(i) - fl(j))

              if (abs(kk) >= first_spch_index) then
                 if (j < 0 .and. kk > 0) then !!f_m < 0 and f_n > 0
                    vv = vv &
                         + smo*fl(kk)*w0*conjg(y(abs(j)+3*(M+1)))*y(kk+3*(M+1))
                    vrho = vrho &
                         + smo*fl(i)*w0*conjg(y(abs(j)+3*(M+1)))*y(kk+4*(M+1))
                 else if (j > 0 .and. kk < 0) then !!f_m > 0 and f_n < 0
                    vv = vv &
                         + smo*fl(kk)*w0*y(j+3*(M+1))*conjg(y(abs(kk)+3*(M+1)))
                    vrho = vrho &
                         + smo*fl(i)*w0*y(j+3*(M+1))*conjg(y(abs(kk)+4*(M+1)))
                 else if (j > 0 .and. kk > 0) then !!f_m > 0, f_n > 0
                    vv = vv + smo*fl(kk)*w0*y(j+3*(M+1))*y(kk+3*(M+1))
                    vrho = vrho + smo*fl(i)*w0*y(j+3*(M+1))*y(kk+4*(M+1))
                 end if
              end if
           !!f_m = f_\ell /= 0 (since i /= M+1), f_n = 0
           else if (j == i .and. abs(j) >= first_spch_index) then
              !! vv gets nothing b/c f_n = 0
              vrho = vrho + smo*fl(i)*w0*y(j+3*(M+1))*y(5*(M+1))
           else if (j == 0) then !! f_m = 0, f_n = f_\ell /= 0 (since i /= M+1)
              vv = vv + smo*fl(i)*w0*y(4*(M+1))*y(i+3*(M+1))
              vrho = vrho + smo*fl(i)*w0*y(4*(M+1))*y(i+4*(M+1))
           end if
        end do

        !!make the rhs's.
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           !! vdot
           thesis_w(i+M+1) = 2.0*pC(z,i)**2*(smo*fl(i)*w0*L(z,i) &
                - R(z,i))*y(i+M+1) + nscrf(z,i)*y(i+2*(M+1)) &
                + smo*fl(i)*w0*y(i+3*(M+1)) - vv
           !! rhodot
           thesis_w(i+3*M+2) = smo*fl(i)*w0*y(i+4*(M+1)) - vrho
        else if (i >= first_spch_index .and. i <= last_spch_index) then
           !! vdot
           thesis_w(i+M+1) = nscrf(z,i)*y(i+2*(M+1)) &
                + smo*fl(i)*w0*y(i+3*(M+1)) - vv
           !! rhodot
           thesis_w(i+3*M+2) = smo*fl(i)*w0*y(i+4*(M+1)) - vrho
        else
           thesis_w(i+M+1) = (0.0,0.0)
           thesis_w(i+3*M+2) = (0.0,0.0)
        end if
     end do

     !! put the negative frequencies into w
     do i = 1, M
        !!rhs of vdot
        thesis_w(M+1-i) = conjg(thesis_w(i+M+1))
        !!rhs of rhodot
        thesis_w(3*M+2-i) = conjg(thesis_w(i+3*M+2))
     end do

     !! dc RHS of v, rho
     thesis_w(M+1) = (0.0,0.0)
     thesis_w(3*M+2) = (0.0,0.0)
     
     !! fill in thesis_V
     ! initialize to zeros
     do i = 1, 2*(2*M+1)
        do j = 1, 2*(2*M+1)
           thesis_V = (0.0,0.0)
        end do
     end do
     !! blocks. from thesis notation i => \ell, j => n
     do i = 1, 2*M+1

        !pm do j = 1, 2*M+1
        do p = 1, pair_matrix(i-M-1,2) !pm
           j = pair_matrix(i-M-1,2*p+1)+M+1 !pm
           kk = pair_matrix(i-M-1,2*p+2) !pm
           !pm kk = find_freq(fl(i-M-1)-fl(j-M-1))

           !pm if (fl(i-M-1) == fl(j-M-1) + fl(kk)) then
           if (kk < 0) then
              !print*, 'i-M-1=',i-M-1,'j-M-1=',j-M-1,'kk=',kk
              !block I
              thesis_V(i,j) = conjg(y(abs(kk)+3*(M+1))) 
              !block II
              thesis_V(i+2*M+1,j) = conjg(y(abs(kk)+4*(M+1))) 
              !block III
              thesis_V(i+2*M+1,j+2*M+1) = thesis_V(i,j)
           else if (kk > 0) then
              !print*, 'i-M-1=',i-M-1,'j-M-1=',j-M-1,'kk=',kk
              !block I
              thesis_V(i,j) = y(kk+3*(M+1))
              !block II
              thesis_V(i+2*M+1,j) = y(kk+4*(M+1))
              !block III
              thesis_V(i+2*M+1,j+2*M+1) = thesis_V(i,j)
           else if (kk == 0 .and. j == i) then
              !print*, 'i-M-1=',i-M-1,'j-M-1=',j-M-1,'kk=',kk
              !block I
              thesis_V(i,j) = y(4*(M+1))
              !block II
              thesis_V(i+2*M+1,j) = y(5*(M+1))
              !block III
              thesis_V(i+2*M+1,j+2*M+1) = thesis_V(i,j)
           end if
           !pm end if
        end do
     end do

     !! solve Ax=b
     call zgesv(2*(2*M+1),1,thesis_V,2*(2*M+1),ipiv,thesis_w,2*(2*M+1),info)

     !! set the derivative
     do i = 1, M
        ydot(i+3*(M+1)) = thesis_w(i+M+1)
        ydot(i+4*(M+1)) = thesis_w(i+3*M+2)
     end do
     ydot(4*(M+1)) = real(thesis_w(M+1))
     ydot(5*(M+1)) = real(thesis_w(3*M+2))

     !! if dc_constant, set the dc derivatives of the velocity and
     !! density to zero
     if (run_parameters % dc_constant) then
        ydot(4*(M+1)) = (0.0,0.0)
        ydot(5*(M+1)) = (0.0,0.0)
     end if

  !*************************************************************!
  !*************************************************************!
  !********************* S-MUSE derivatives ********************!
  !*************************************************************!
  !*************************************************************!
  else if (modelID=='S') then
     !! put DC to zero
     ydot(M+1) = (0.0,0.0)
     ydot(2*(M+1)) = (0.0,0.0)
     ydot(3*(M+1)) = (0.0,0.0)
     ydot(4*(M+1)) = (0.0,0.0)
     ydot(5*(M+1)) = (0.0,0.0)

     !circuit and space charge derivatives
     do i = 1, M
        !! circuit derivatives
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           i_freq = smo*fl(i)*w0
           !! voltage derivatives
           ydot(i) = -i_freq*y(i) + (R(z,i)-i_freq*L(z,i))*y(i+M+1)

           !! diagnostics to see if mathematica is using same values
           !print*, -i_freq / circuit_parameters % circuit_length
           !print*, -i_freq*K(z*circuit_parameters % circuit_length,i) &
           !*derived_qtys % u0 &
           !/(vph(z*circuit_parameters % circuit_length,i) &
           !*circuit_parameters % circuit_length)

           !! current derivatives
           ydot(i+M+1) = (G(z,i) - i_freq*Ca(z,i))*y(i) &
                - i_freq*y(i+M+1) + i_freq*pC(z,i)*y(i+4*(M+1))

           !! diagnostics to see if mathematica is using same values
           !print*, -i_freq*derived_qtys % u0 &
           !/(K(z*circuit_parameters % circuit_length,i) &
           !*vph(z*circuit_parameters % circuit_length,i) &
           !*circuit_parameters % circuit_length)
           !print*, i_freq* derived_qtys % u0 * derived_qtys % beam_area &
           !/circuit_parameters % circuit_length

        else
           ydot(i) = (0.0,0.0)
           ydot(i+M+1) = (0.0,0.0)
        end if

        !! space charge field derivative
        if (i >= first_spch_index .and. i <= last_spch_index) then
           !! space charge derivatives
           ydot(i+2*(M+1)) = -i_freq*y(i+2*(M+1)) + y(i+4*(M+1))
           !! diagnostics to see if mathematica is using same values
           !print*, 1/eps0
        else
           ydot(i+2*(M+1)) = (0.0,0.0)
        end if


        !! make nonlinear terms
        vv = (0.0,0.0)
        Irho = (0.0,0.0)
        Erho = (0.0,0.0)
        vrho = (0.0,0.0)

        !pm comments starting with pm were put in when inserting pair_matrix
        !pm p = 0
        !pm do q = -last_spch_index, last_spch_index
        do j = 1, pair_matrix(i,2)  !pm
           q = pair_matrix(i,2*j+1) !pm
           p = pair_matrix(i,2*j+2) !pm

           !pm p = find_freq(fl(i) - fl(q))

           if (p /= 0 .and. q /= 0 .and. abs(q) >= first_spch_index &
                .and. abs(p) >= first_spch_index &
                .and. abs(p) <= last_spch_index) then
              if (p < 0) then
                 vv = vv + smo*fl(p)*w0*conjg(y(abs(p)+3*(M+1)))*y(q+3*(M+1))
                 if (abs(p) >= first_ckt_index) then
                    Irho = Irho + 2*pC(z,p)**2*(R(z,p)-smo*fl(p)*w0*L(z,p)) &
                         *conjg(y(abs(p)+M+1))*y(q+4*(M+1))
                 end if
                 Erho = Erho + nscrf(z,p)*conjg(y(abs(p)+2*(M+1)))*y(q+4*(M+1))
                 vrho = vrho &
                      + smo*fl(i)*w0*conjg(y(abs(p)+3*(M+1)))*y(q+4*(M+1))
              else if (q < 0 ) then
                 vv = vv + smo*fl(p)*w0*y(p+3*(M+1))*conjg(y(abs(q)+3*(M+1)))
                 if (abs(p) >= first_ckt_index) then
                    Irho = Irho + 2*pC(z,p)**2*(R(z,p)-smo*fl(p)*w0*L(z,p)) &
                         *y(p+M+1)*conjg(y(abs(q)+4*(M+1)))
                 end if
                 Erho = Erho + nscrf(z,p)*y(p+2*(M+1))*conjg(y(abs(q)+4*(M+1)))
                 vrho = vrho &
                      + smo*fl(i)*w0*y(p+3*(M+1))*conjg(y(abs(q)+4*(M+1)))
              else
                 vv = vv + smo*fl(p)*w0*y(p+3*(M+1))*y(q+3*(M+1))
                 if (abs(p) >= first_ckt_index) then
                    Irho = Irho + 2*pC(z,p)**2*(R(z,p)-smo*fl(p)*w0*L(z,p)) &
                         *y(p+M+1)*y(q+4*(M+1))
                 end if
                 Erho = Erho + nscrf(z,p)*y(p+2*(M+1))*y(q+4*(M+1))
                 vrho = vrho + smo*fl(i)*w0*y(p+3*(M+1))*y(q+4*(M+1))
              end if
           end if
        end do

        !! velocity derivative
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           ydot(i+3*(M+1)) = 2.0*pC(z,i)**2*(smo*fl(i)*w0*L(z,i) &
                - R(z,i))*y(i+M+1) + nscrf(z,i)*y(i+2*(M+1)) - vv
        else if (i >= first_spch_index .and. i <= last_spch_index) then
           ydot(i+3*(M+1)) = nscrf(z,i)*y(i+2*(M+1)) - vv
        else
           ydot(i+3*(M+1)) = (0.0,0.0)
        end if

        !! density derivative
        if (i >= first_ckt_index .and. i <= last_ckt_index) then
           ydot(i+4*(M+1)) = -2.0*pC(z,i)**2*(smo*fl(i)*w0*L(z,i) &
                - R(z,i))*y(i+M+1) - nscrf(z,i)*y(i+2*(M+1)) &
                - smo*fl(i)*w0*y(i+3*(M+1)) + Irho - Erho + vv - vrho
        else if (i >= first_spch_index .and. i <= last_spch_index) then
           ydot(i+4*(M+1)) = -nscrf(z,i)*y(i+2*(M+1)) &
                - smo*fl(i)*w0*y(i+3*(M+1)) + Irho - Erho + vv - vrho
        else
           ydot(i+4*(M+1)) = (0.0,0.0)
        end if
     end do
  end if
  
end subroutine F

subroutine JACOBN(n, t, y, dfdy, matdim, ml, mu)
  use precision95
  !!jacobian for right hand sides
  use working_variables, only : modelID
  use working_variables, only : M

  implicit none

  integer, intent(in) :: n
  real(dp), intent(in) :: t
  complex(dp), dimension(1:n), intent(in) :: y
  complex(dp), dimension(1:n,1:n), intent(inout) :: dfdy
  integer, intent(in) :: matdim
  integer, intent(in) :: ml
  integer, intent(in) :: mu


  if (modelID=='L') then
     print*, 'compute LATTE jacobian'
  else if (modelID=='M') then
     print*, 'compute MUSE jacobian'
  else if (modelID=='S') then
     print*, 'compute S-MUSE jacobian'
  end if
end subroutine JACOBN

subroutine FA(n, t, y, A, matdim, ml, mu, nde)
  use precision95
  !!computes A matrix of A(y,t)*dy/dt = F(y,t)
  implicit none

  integer, intent(in) :: n
  real(dp), intent(in) :: t
  complex(dp), dimension(1:n), intent(in) :: y
  complex(dp), dimension(1:n,1:n), intent(inout) :: A
  integer, intent(in) :: matdim
  integer, intent(in) :: ml
  integer, intent(in) :: mu
  integer, intent(in) :: nde

  !!A assumed to be full, unless has banded structure?
  
end subroutine FA

subroutine G(n, t, y, iroot)
  use precision95
  !!subroutine for using root solver. can use root solver to stop at certain
  !!values of output, etc
  implicit none

  integer, intent(in) :: n
  real(dp), intent(in) :: t
  complex(dp), dimension(1:n), intent(in) :: y
  integer, intent(inout) :: iroot
  
end subroutine G

subroutine USERS(y, yh, ywt, save1, save2, t, h, el, impl, n, nde, iflag)
  use precision95
  !!called by cdriv3.f to solve certain linear systems
  !!e.g. to use sparse methods
  implicit none

  integer, intent(in) :: n
  complex(dp), dimension(1:n), intent(inout) :: y
  complex(dp), dimension(1:n), intent(inout) :: yh
  complex(dp), dimension(1:n), intent(inout) :: ywt
  complex(dp), dimension(1:n), intent(inout) :: save1
  complex(dp), dimension(1:n), intent(inout) :: save2
  real(dp), intent(inout) :: t
  real(dp), intent(inout) :: h
  real(dp), intent(inout) :: el
  integer, intent(in) :: impl
  integer, intent(in) :: nde
  integer, intent(in) :: iflag

end subroutine USERS
