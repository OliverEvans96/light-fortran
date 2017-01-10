! File Name: rte_core.f90
! Description: Subprograms applicable to RTE in general
! Created: Fri Jan 06, 2017 | 10:54am EST
! Last Modified: Tue Jan 10, 2017 | 02:06pm EST

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

module rte_core

use utils

contains

! Calculate evenly spaced VSF array (beta)
function calc_vsf_arr(vsf_file,data_rows,lmax,fmtstr_in,skiplines_in)
    !use utils
    implicit none

    ! INPUTS:
    ! vsf_file - path to file containing 2D array of experimental VSF data
    character(len=*), intent(in) :: vsf_file
    ! data_rows - size of VSF data
    ! lmax - size of phi array
    integer, intent(in) :: data_rows, lmax
    ! fmtstr - optional - data format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), optional :: fmtstr_in
    character(len=20) :: fmtstr
    ! skiplines - optional - number of lines to skip from data header
    integer, optional :: skiplines_in
    integer :: skiplines

    ! OUTPUT:
    ! evenly spaced VSF array of size lmax-1
    double precision, dimension(lmax-1) :: calc_vsf_arr

    ! BODY:

    ! phi step size
    double precision dphi
    ! Array from data file
    double precision, dimension(data_rows,2) :: vsf_data
    ! Angle at which to evaluate (interpolate) VSF
    double precision phi
    ! Phi index
    integer ll

    ! Default format string
    if(present(fmtstr_in)) then
        fmtstr = fmtstr_in
    else
        fmtstr = 'E10.2'
    end if

    ! Skip no lines if not specified
    if(present(skiplines_in)) then
        skiplines = skiplines_in
    else
        skiplines = 0
    end if

    ! Calculate phi step size
    dphi = pi / lmax

    ! Read data file
    vsf_data = read_array(trim(vsf_file),trim(fmtstr),data_rows,2,skiplines)

    ! Loop through phi values
    do ll = 1, lmax - 1

        ! Calculate angle
        phi = ll * dphi

        ! Interpolate data
        calc_vsf_arr(ll) = interp(phi, vsf_data(:,1), vsf_data(:,2), data_rows)
    end do

end function

end module
