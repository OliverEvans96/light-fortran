! File Name: rte_utils.f90
! Description: General utilities which might be useful in other settings
! Created: Wed Jan 04, 2017 | 06:24pm EST
! Last Modified: Fri Jan 06, 2017 | 10:55am EST

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

! General utilities which might be useful in other settings
module rte_utils
contains

! Interpolate single point from 1D data
function interp(x0,xx,yy,nn)
    implicit none

    ! INPUTS:
    ! x0 - x value at which to interpolate
    double precision, intent(in) :: x0
    ! xx - ordered x values at which y data is sampled
    ! yy - corresponding y values to interpolate
    double precision, dimension (nn), intent(in) :: xx,yy    
    ! nn - length of data
    integer, intent(in) :: nn

    ! OUTPUT:
    ! interp - interpolated y value
    double precision interp

    ! BODY:

    ! Index of lower-adjacent data (xx(i) < x0 < xx(i+1))
    integer ii
    ! Slope of liine between (xx(ii),yy(ii)) and (xx(ii+1),yy(ii+1))
    double precision mm

    ! Determine ii
    do ii = 1, nn
        if (xx(ii) > x0) then
            ! We've now gone one index too far.
            exit
        end if
    end do

    ! Subtract since we went one index too far
    ii = ii - 1

    ! Calculate slope
    mm = (yy(ii+1) - yy(ii)) / (xx(ii+1) - xx(ii))

    ! Return interpolated value
    interp = yy(ii) + mm * (x0 - xx(ii))
end function

! Read 2D array from file
function read_array(filename,fmtstr,nn,mm)
    implicit none

    ! INPUTS:
    ! filename - path to file to be read
    ! fmtstr - input format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), intent(in) :: filename, fmtstr
    ! nn - Number of data rows in file
    ! mm - number of data columns in file
    integer, intent(in) :: nn, mm

    ! OUTPUT:
    double precision, dimension(nn,mm) :: read_array
    
    ! BODY:

    ! Row counter
    integer ii
    ! File unit number
    integer, parameter :: un = 10
    ! Final format to use
    character(len=256) finfmt

    ! Generate final format string
    write(finfmt,'(A,I4,A,A)') '(', mm, fmtstr, ')'

    ! Print message
    write(*,*) 'Reading data from "', trim(filename), '"'

    ! Open file
    open(unit=un, file=trim(filename), status='old', form='formatted')

    ! Loop through lines
    do ii = 1, nn
        ! Read one row at a time
        read(unit=un, fmt=finfmt) read_array(ii,:)
    end do

    ! Close file
    close(unit=un)

end function

! Print 2D array to stdout
subroutine print_array(arr,nn,mm,fmtstr)
    implicit none

    ! INPUTS:
    ! arr - array to print
    double precision, dimension (nn,mm), intent(in) :: arr
    ! nn - number of data rows in file
    ! nn - number of data columns in file
    integer, intent(in) :: nn, mm
    ! fmtstr - output format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), optional :: fmtstr
    
    ! NO OUTPUTS

    ! BODY
    
    ! Row counter
    integer ii
    ! Final format to use
    character(len=256) finfmt

    ! Determine string format
    if(.not. present(fmtstr)) then
        fmtstr = 'E10.2'
    end if

    ! Generate final format string
    write(finfmt,'(A,I4,A,A)') '(', mm, fmtstr, ')'

    ! Loop through rows
    do ii = 1, nn
        ! Print one row at a time
        write(*,finfmt) arr(ii,:) 
    end do

    ! Print blank line after
    write(*,*) ' '

end subroutine

end module
