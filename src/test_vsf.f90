! File Name: test_vsf.f90
! Description: Test calc_vsf_array
! Created: Wed Jan 04, 2017 | 06:47pm EST
! Last Modified: Tue Jan 10, 2017 | 11:51am EST
! Author: Oliver Evans <oge1@zips.uakron.edu>

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

program test_vsf
    use rte_core
    implicit none

    ! VSF file to read
    character(len=256) vsf_file
    ! File length
    !integer, parameter :: data_rows = 55
    integer, parameter :: data_rows = 48
    ! Desired number of points to interpolate
    integer, parameter :: lmax = 10
    ! Desired step size
    !double precision, parameter :: dphi = 1.D-1
    ! Fata format
    character(len=256) fmtstr
    ! Output VSF array
    double precision, dimension(lmax) :: beta
    ! angular array
    double precision, dimension(lmax) :: phi
    ! array with both angle and vsf
    double precision, dimension(data_rows,2) :: out_arr

    write(*,*) 'Begin!'

    ! File to read
    !vsf_file = trim(getbasedir()) // '/data/vsf/nuc_vsf.txt'
    vsf_file = trim(getbasedir()) // '/data/sin4x_12rows_10.2e.txt'
    ! Data format
    !fmtstr = 'E13.4'
    !fmtstr = 'F13.7'
    fmtstr = 'E10.2'

    ! Calculate phi
    write(*,*) 'bnd2arr'
    phi = bnd2arr(0.D0, 2*pi, lmax)

    ! Calculate VSF
    write(*,*) 'calc_vsf_arr'
    out_arr = read_array(trim(vsf_file),'E10.2',48,2)
    call print_array(out_arr,48,2)
    beta = calc_vsf_arr(vsf_file, data_rows, lmax, trim(fmtstr), 0)
    call print_array(beta,48,2)

    ! Generate output array
    !out_arr(:,1) = phi
    !out_arr(:,2) = beta

    ! Print result
    !call print_array(out_arr, lmax, 2)

end program
