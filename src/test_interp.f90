! File Name: test_interp.f90
! Description: Test interp.f90
! Created: Wed Jan 04, 2017 | 06:47pm EST
! Last Modified: Sun Jan 15, 2017 | 06:00pm EST
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
    ! Array dimensions
    integer, parameter :: nn=48, mm=2
    ! Array variable
    double precision, dimension(nn,mm) :: arr

    ! File to read
    filename = trim(getbasedir()) // '/data/test/sin4x_12rows_10.2e.txt'

    ! Read data file
    !write(*,*) 'Reading array from "', filename, '"'
    arr = read_array(trim(filename),'E10.2',nn,mm)

    ! Print array
    write(*,*) 'Printing array'
    call print_array(arr,nn,mm,'E10.2')

    ! Interpolate
    write(*,*) 'Testing interpolation'
    write(*,'(A,F5.2)') 'f(1.3) = ', interp(1.3D0,arr(:,1),arr(:,2),12)

end program
