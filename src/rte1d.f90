! File Name: rte1d.f90
! Description: Subprograms specific to 3D plane parallel RTE
! Created: Thu Jan 23, 2017 | 05:52pm EST
! Last Modified: Mon Jan 23, 2017 | 07:30pm EST

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!                           GNU GPL LICENSE                            !
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
!                                                                      !
! Copyright Oliver Evans 2017 <oliverevans96@gmail.com>                !
!                                                                      !
! This program is free software: you can redistribute it and/or modify !
! it under the terms of the GNU General Public License as published by !
! the Free Software Foundation, either version 3 of the License, or    !
! (at your option) any later version.                                  !
!                                                                      !
! This program is distributed in the hope that it will be useful,      !
! but WITHOUT ANY WARRANTY; without even the implied warranty of       !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         !
! GNU General Public License for more details.                         !
!                                                                      !
! You should have received a copy of the GNU General Public License    !
! along with this program. If not, see <http://www.gnu.org/licenses/>. !
!                                                                      !
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

module rte1d

use rte_core

contains

! Calculate the angle between two unit vectors in the directions
! (theta,phi) and (theta_p[rime],phi_p[rime])
function anglediff(theta,phi,theta_p,phi_p)
    implicit none

    ! INPUTS
    ! theta - azimuthal angle #1
    double precision theta
    ! phi - polar angle #1
    double precision phi
    ! theta - azimuthal angle #2
    double precision theta_p
    ! phi - polar angle #2
    double precision phi_p

    ! OUTPUT
    ! anglediff - angle between two vectors
    double precision anglediff

    anglediff = acos(sin(phi)*sin(phi_p)*cos(theta-theta_p) + cos(phi)*cos(phi_p))
end function


! Checkerboard successive over-relaxation for the 2D RTE
! Periodic boundary conditions in the x direction
subroutine sor(rad, aa, bb, beta, jmax, lmax, mmax, vmax&
                zmin, zmax, tol, maxiter, omega)
    implicit none

    ! INPUTS:
    ! rad - radiance array to be updated.
    ! should be already set with some initial guess
    double precision, dimension(imax,jmax,lmax) :: rad
    ! jmax, lmax, mmax - array size in z, theta, phi dimensions
    integer, intent(in) :: jmax, lmax, mmax
    ! vmax - size of VSF data array
    integer, intent(in) :: vmax
    ! aa - absorption coefficient over space
    ! bb - scattering coefficient over space
    double precision, dimension(jmax), intent(in) :: aa, bb
    ! beta - normalized volume scattering function evenly spaced array
    ! beta_1 = vsf(dphi); vsf(0) not meaningful
    double precision, dimension(vmax), intent(in) :: beta
    ! zmin, zmax - bounds for z dimension
    double precision, intent(in) :: zmin, zmax
    ! tol - error tolerance which determines when to stop SOR iteration
    double precision, optional :: tol
    ! maxiter - maximum number of SOR iterations
    integer, optional :: maxiter
    ! omega - SOR weighting term. Should be between 1 and 2
    double precision, optional :: omega

    ! Step size in z, theta, and phi dimensions
    double precision dz, dtheta, dphi

    ! OUTPUTS:
    ! rad will be updated

    ! BODY:

    ! xx, yy, phi, phi' variables and indices
    double precision xx, yy, phi, phi_p
    integer ii, jj, ll, lp

    ! Iteration counter
    integer iter

    ! Flag to distinguish between odd checkerboard loop and even
    ! 1 for odd, 2 for even
    integer oddeven

    ! CD2 cartesian derivatives
    double precision drdx, drdy

    ! Radial derivative
    double precision drdr

    ! Attenuation coefficient - sum of absorption & scattering
    double precision, dimension(imax,jmax) :: cc

    ! Radiative source term
    double precision source

    ! Gauss-Seidel term
    double precision gs_term

    ! Temporary radiance value for checking error
    double precision tmp_rad

    ! Error (difference between iterations)
    double precision err

    ! Angular variable bounds
    double precision, parameter :: thetamin = 0.D0
    double precision, parameter :: phimin = 0.D0
    double precision, parameter :: thetamax = 2*pi
    double precision, parameter :: phimax = pi

    ! Calculate step size
    dz = (zmax - zmin) / imax
    dtheta = (thetamax - thetamin) / lmax
    dphi = (phimax - phimin) / mmax

    ! Calculate attenuation coefficient
    cc = aa + bb

    ! Assign optional parameters
    if(.not. present(tol)) then
        ! 1.0D2 = 100
        tol = 1.0D2
    end if
    if(.not. present(maxiter)) then
        maxiter = 1000
    end if
    if(.not. present(omega)) then
        ! 1.75D0 = 1.75
        omega = 1.75D0
    end if

    ! Loop through iterations
    do iter = 1, maxiter

        ! Reset error
        err = 0

        ! Checkerboard - perform SOR on odd cells first, then even
        ! oddeven = 1: include upper left
        ! oddeven = 2: Do not include upper left
        do oddeven = 1, 2

            !!! Parallelize this loop !!!
            ! Loop through jj
            ! Start at column 2 (depth layer 2) to accomodate BC
            do jj = 2, jmax

                    write(*,'(A,I3,I3,I3)') 'LOC ', iter, ii, jj

                    ! Loop through phi values
                    ! Do not parallelize
                    do ll = 1, lmax

                        ! Calculate phi
                        phi = (ll-1) * dphi

                        ! Calculate derivatives
                        ! Periodic in x direction
                        drdx = (rad(mod(ii+1,imax),jj,ll) &
                                - rad(mod(ii-1,imax),jj,ll)) / (2*dx)

                        ! y derivative
                        ! Interior
                        if(jj .lt. jmax) then
                            ! CD2 for interior
                            drdy = (rad(ii,jj+1,ll) - rad(ii,jj-1,ll)) / (2*dy)
                        ! Bottom (jj = jmax; greatest depth)
                        else
                            ! BD2 for downwelling radiance
                            if(ll .lt. lmax/2) then
                                drdy = (3*rad(ii,jj,ll) - 4*rad(ii,jj-1,ll) &
                                      + rad(ii,jj-2,ll)) / (2 * dy)
                            ! No upwelling radiance on bottom (phi >= phimax/2; ll >= lmax/2)
                            else
                                drdy = 0
                            end if
                        end if

                        ! Radial derivative in the direction of ray
                        drdr = drdx * cos(phi) + drdy * sin(phi)

                        ! Calculate source term
                        ! Set to zero, then integrate
                        source = 0
                        ! SOURCE ZEROED OUT
                        do lp = 1, lmax
                            ! Exclude scattering from the current direction
                            if(lp .eq. ll) then
                                cycle
                            end if

                            ! Integrate (left endpoint rule)
                            source = source + beta(diff_ind(ll,lp,lmax)) &
                                    * rad(ii,jj,lp) * dphi
                        end do

                        ! Calculate Gauss-Seidel term
                        gs_term = 1/cc(ii,jj) * (bb(ii,jj)*source - drdr)

                        ! Calculate SOR correction and update
                        tmp_rad = (1-omega) * rad(ii,jj,ll) &
                                        +  omega  * gs_term

                        ! Accumulate error
                        err = err + abs(tmp_rad - rad(ii,jj,ll))

                        ! Update radiance
                        rad(ii,jj,ll) = tmp_rad

                        ! Print radiance
                        write(*,'(A,I4,I4,I4,A,E13.4E3)') &
                            '(it,ii,jj,ll): rad = ', ii, jj, ll, &
                            ': ', rad(ii,jj,ll)

                    end do
                end do
            end do
        end do

        ! Check whether tolerance has been met
        write(*,'(A,I3,A,E10.3)') 'After iteration ', iter, &
                ', err = ', err
        write(*,*)
        if(err < tol) then
            exit
        end if
    end do

end subroutine


end module
