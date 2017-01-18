! File Name: rte2d.f90
! Description: Subprograms specific to 2D RTE
! Created: Thu Jan 05, 2017 | 06:30pm EST
! Last Modified: Wed Jan 18, 2017 | 05:38pm EST

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

module rte2d

use rte_core

contains

! Calculate index of beta array which corresponds to beta(phi(ll),phi(lp))
! For use in VSF integration to calculate scattering
function diff_ind(ll,lp,lmax)
    implicit none

    ! INPUTS:
    ! ll, lp - two angular indices
    ! lmax - total number of angles
    ! lmax should be even. Not sure what happens otherwise
    integer, intent(in) :: ll, lp, lmax

    ! OUTPUT:
    ! diff_ind - index of difference
    integer diff_ind

    ! Calculate index
    diff_ind = mod(ll - lp, lmax/2 + 1)
end function

! Checkerboard successive over-relaxation for the 2D RTE
! Periodic boundary conditions in the x direction
subroutine sor(rad, aa, bb, beta, imax, jmax, lmax, &
                xmin, xmax, ymin, ymax, tol, maxiter, omega)
    use utils
    implicit none

    ! INPUTS:
    ! rad - radiance array to be updated.
    ! should be already set with some initial guess
    double precision, dimension(imax,jmax,lmax) :: rad
    ! imax, jmax, lmax - array size in x, y, phi dimensions
    integer, intent(in) :: imax, jmax, lmax
    ! aa - absorption coefficient over space
    ! bb - scattering coefficient over space
    double precision, dimension(imax,jmax), intent(in) :: aa, bb
    ! beta - normalized volume scattering function evenly spaced array
    ! beta_1 = vsf(dphi); vsf(0) not meaningful
    double precision, dimension(lmax-1), intent(in) :: beta
    ! xmin, xmax, ymin, ymax - bounds for x and y dimensions respectively
    double precision, intent(in) :: xmin, xmax, ymin, ymax
    ! tol - error tolerance which determines when to stop iterating
    double precision, optional :: tol
    ! maxiter - maximum number of SOR iterations
    integer, optional :: maxiter
    ! omega - SOR weighting term. Should be between 1 and 2
    double precision, optional :: omega

    ! Step size in each dimension
    double precision dx, dy, dphi

    ! OUTPUTS:
    ! rad will be updated

    ! BODY:

    ! xx, yy, phi, phi' variables and indices
    double precision xx, yy, phi, phi_p
    integer ii, jj, ll, lp

    ! Array size in each dimension (x, y, phi respectively)
    !integer imax, jmax, lmax

    ! Radiance array from previous iteration for tolerance comparison
    double precision, dimension(imax,jmax,lmax) :: prev_rad

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

    ! Calculate step size
    dx = (xmax - xmin) / imax
    dy = (ymax - ymin) / jmax
    dphi = 2 * pi / lmax

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

    ! Copy radiance array to prev_rad
    prev_rad = rad

    ! Loop through iterations
    do iter = 1, maxiter

        ! Checkerboard - perform SOR on odd cells first, then even
        do oddeven = 1, 2

            ! Loop through odd/even jj values
            !!! Parallelize this loop !!!
            !! Add 1 to oddeven to start SOR at 2nd depth layer to accomodate BC
            do jj = oddeven+1, jmax, 2

                ! Calculate yy
                !yy = (jj-1) * dy

                ! Loop through even ii values if jj is even
                ! or odd ii values if jj is odd
                !!! Parallelize this loop !!!
                do ii = 2-mod(jj,2), imax, 2
                    write(*,'(I3,I3)') iter, jj

                    ! Loop through phi values
                    ! Do not parallelize
                    do ll = 1, lmax

                        ! Calculate phi
                        phi = (ll-1) * dphi

                        ! Calculate xx
                        !xx = (ii-1) * dx

                        ! Calculate derivatives
                        ! Periodic in x direction
                        drdx = (rad(mod(ii+1,imax),jj,ll) &
                                - rad(mod(ii-1,imax),jj,ll)) / (2*dx)

                        ! No upwelling radiance on bottom (highest j index)
                        if((jj .lt. jmax) .and. (ll .gt. lmax/2)) then
                            drdx = (rad(ii,jj+1,ll) - rad(ii,jj-1,ll)) / (2*dy)
                        else
                            drdx = (0 - rad(ii,jj-1,ll)) / (2*dy)
                        end if

                        ! Radial derivative in the direction of ray
                        drdr = drdx * cos(phi) + drdy * sin(phi)

                        ! Calculate source term
                        ! Set to zero, then integrate
                        source = 0
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
                        rad(ii,jj,ll) = (1-omega) * rad(ii,jj,ll) &
                                        +  omega  * gs_term
                    end do
                end do
            end do
        end do

        ! Check whether tolerance has been met
        write(*,'(A,I3,A,E10.3)') 'After iteration ', iter, &
                ', err = ', sum(rad-prev_rad)
        write(*,*)
        if(sum(rad-prev_rad) < tol) then
            exit
        end if

        ! Update radiance
        prev_rad = rad
    end do

end subroutine


end module
