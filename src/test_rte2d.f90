! File Name: test_rte2d.f90
! Description: Test RTE 2D w/ checkerboard SOR, periodic x
! Created: Tue Jan 10, 2017 | 02:54pm EST
! Last Modified: Sun Jan 29, 2017 | 05:10pm EST

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
    implicit none

    ! Spatial/angular counters
    integer ii, jj, ll

    ! Grid size (prescribed)
    integer, parameter :: imax = 20
    integer, parameter :: jmax = 20
    integer, parameter :: lmax = 20

    ! Grid extent (prescribed)
    double precision, parameter :: xmin = 0, xmax = 1
    double precision, parameter :: ymin = 0, ymax = 1
    double precision, parameter :: phimin = 0, phimax = 2*pi

    ! Grid mesh size (calculated)
    double precision dx, dy, dphi

    ! Surface boundary condition
    double precision, dimension(lmax/2) :: surf_bc

    ! aa - Absorption coefficient
    ! bb - Scattering coefficient
    double precision, dimension(imax,jmax) :: aa, bb

    ! Volume scattering function array (normalized)
    double precision, dimension(lmax-1) :: beta, norm_beta

    ! Absorption/scattering coefficient values (constant over space)
    double precision, parameter :: a_val = 1
    double precision b_val

    ! Volume scattering function data file location, format & # of rows
    character(len=256) vsf_file
    character(len=256) vsf_fmt
    integer, parameter :: vsf_data_rows = 55

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

    ! BODY:

    ! Calculate mesh sizes
    dx = (xmax - xmin) / (imax - 1)
    dy = (ymax - ymin) / (jmax - 1)
    dphi = (phimax - phimin) / lmax

    ! VSF file location & data format
    vsf_file = trim(getbasedir()) // '/data/vsf/nuc_vsf.txt'
    vsf_fmt = 'E13.4'

    ! Generate VSF array
    beta = calc_vsf_arr(trim(vsf_file), vsf_data_rows, lmax, trim(vsf_fmt), 1)
    norm_beta = beta
    ! Determine scattering coefficient & normalize array
    b_val = normalize(norm_beta,dphi,lmax-1)

    ! Fill coefficient arrays
    do ii = 1, imax
        do jj = 1, jmax
            aa(ii,jj) = a_val
            bb(ii,jj) = b_val
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

    ! Read boundary conditions from file
    ! Each row is the radiance value at a corresponding angle
    ! Angle is measured from the x-axis, which points right
    ! towards the downward-pointing y-axis.
    ! Refer to fig. 2 in summary paper

    ! read_array generates imax/2 x 1 array. Reshape this into
    ! an array of rank 1 and length imax.
    surf_bc = reshape(read_array(trim(getbasedir())//'/data/surf_bc/surf_bc_50.txt','E13.4',lmax/2,1),(/lmax/2/))

    ! Apply surface boundary conditions (constant over space)
    do ii = 1, imax
        do ll = 1, lmax/2
            rad(ii,1,ll) = surf_bc(ll)
        end do
    end do

    ! Read beta

    ! Perform SOR
    write(*,*) 'SOR'
    call sor(rad, aa, bb, beta, imax, jmax, lmax, &
            xmin, xmax, ymin, ymax, tol, maxiter, omega)

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
    call write_array(irrad, imax, jmax, trim(getbasedir())//'/results/irrad.txt', 'E13.4E3')
    write(*,*) 'Done RTE!'

end program


