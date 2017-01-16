! File Name: test_rte2d.f90
! Description: Test RTE 2D w/ checkerboard SOR, periodic x
! Created: Tue Jan 10, 2017 | 02:54pm EST
! Last Modified: Sun Jan 15, 2017 | 06:00pm EST

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

! Test RTE2D
program test_rte2d
    use rte2d

    ! Box dimensions
    integer, parameter :: imax = 100
    integer, parameter :: jmax = 100
    integer, parameter :: lmax = 100

    ! Surface boundary condition
    double precision, dimension(imax) :: surf_bc

    ! Absorption coefficient
    double precision, parameter :: aa = 1
    ! Scattering coefficient
    double precision, parameter :: bb = 1
    ! Volume Scattering function array (normalized)
    double precision beta

    ! Radiance array
    double precision, dimension(imax,jmax,lmax) :: rad
    ! Irradiance array
    double precision, dimension(imax,jmax) :: irrad

    ! SOR Parameters
    ! Tolerance
    double precision, parameter :: tol = 1.D2
    ! Maximum # of iterations
    integer, parameter :: maxiter = 100
    ! Overrelaxation parameter
    double precision, parameter :: omega = 1.75D0


    ! Read boundary conditions from file
    ! Each row is the radiance value at a corresponding angle
    ! Angle is measured from the x-axis, which points right
    ! towards the downward-pointing y-axis.
    ! Refer to fig. 2 in summary paper
    surf_bc = reshape(read_array(trim(getbasedir())//'/data/surf_bc/surf_bc_50.txt',imax,1),1)

    ! Apply boundary conditions (constant over space)
    do ii = 1, imax
        do ll = 1, int(lmax/2)
            rad(ii,1,ll) = surf_bc(ii,jj)
        end do
    end do

    write(*,*) 'Initial guess'
    ! Initial guess of 0 everywhere
    do ll = 1, lmax
        do jj = 1, jmax
            do ii = 1, imax
                rad(ii,jj,ll) = 0
            end do
        end do
    end do


    ! Perform SOR
    write(*,*) 'SOR'
    call sor(rad, aa, bb, beta, imax, hmax, lmax, tol, maxiter, omega)

    ! Calculate irradiance at each point in space
    write(*,*) 'Irradiance'
    do jj = 1, jmax
        do ii = 1, imax
            ! Integrate radiance over all angles w/ left endpoint rule
            irrad(ii,jj) = lep_rule(rad(ii,jj,:), dphi, lmax)
        end do
    end do

    write(*,*) trim(getbasedir())//'/results/irrad.txt'
    ! Write irradiance to file
    write_array(irrad, imax, jmax, trim(getbasedir())//'/results/irrad.txt', '10.4E')
    write(*,*) 'Done RTE!'

end program


