! File Name: test_interp.f90
! Description: Test interp.f90
! Created: Wed Jan 04, 2017 | 06:47pm EST
! Last Modified: Thu Jan 19, 2017 | 04:12pm EST
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

program test_interp
    use utils
    implicit none

    ! File to read
    character(len=256) filename
    ! Input array dimensions
    integer, parameter :: nn=101, mm=2
    ! Number of interpolation points
    integer, parameter :: nout = 150
    ! Input array variable
    double precision, dimension(nn,mm) :: arr
    ! Output array variable
    double precision, dimension(nout,mm) :: out_arr
    ! Interpolation step size
    double precision dx
    ! Input data bounds
    double precision xmin, xmax
    ! Index variable
    integer ii

    ! File to read
    filename = trim(getbasedir()) // '/data/test/exp_101rows_10.3e.txt'

    ! Read data file
    !write(*,*) 'Reading array from "', filename, '"'
    arr = read_array(trim(filename),'E10.3',nn,mm)

    ! Print array
    write(*,*) 'Printing array'
    call print_array(arr,nn,mm,'E13.4')

    ! Calculate interpolation points *INCLUDE ENDPOINTS*
    xmin = arr(1,1)
    xmax = arr(nn,1)
    dx = (xmax - xmin) / (nout - 1)
    do ii = 1, nout
        out_arr(ii,1) = xmin + (ii-1) * dx
    end do

    ! Interpolate
    write(*,*) 'Testing interpolation'
    do ii = 1, nout
        out_arr(ii,2) = interp(out_arr(ii,1),arr(:,1),arr(:,2),nn)
    end do

    ! Write to file
    write(*,*) 'Test write to file'
    call write_array(out_arr,nout,mm, &
            trim(getbasedir())//'/results/test_interp.txt','E13.4')
    write(*,*) 'Done!'

end program
