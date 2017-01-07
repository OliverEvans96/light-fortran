! File Name: utils.f90
! Description: General utilities which might be useful in other settings
! Created: Wed Jan 04, 2017 | 06:24pm EST
! Last Modified: Fri Jan 06, 2017 | 07:22pm EST

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
module utils

! Constants
double precision, parameter :: pi = 4.D0 * datan(1.D0)

contains

! Determine array size from min, max and step
! If alignment is off, array will overstep the maximum

! Determine array size from min, max and step
! If alignment is off, array will overstep the maximum
function bnd2max(xmin,xmax,dx)
    ! INPUTS:
    ! xmin - minimum x value in array
    ! xmax - maximum x value in array (inclusive)
    ! dx - step size
    double precision, intent(in) :: xmin, xmax, dx
    
    ! OUTPUT:
    ! step2max - maximum index of array
    integer step2max

    ! Calculate array size
    step2max = int(ceiling((xmax-xmin)/dx))
end function

! Create array from bounds and number of elements
function bnd2arr(xmin,xmax,imax)
    ! INPUTS:
    ! xmin - minimum x value in array
    ! xmax - maximum x value in array (inclusive)
    double precision, intent(in) :: xmin, xmax
    ! imax - number of elements in array
    integer imax
    
    ! OUTPUT:
    ! bnd2arr - array to generate
    double precision, dimension(imax) :: bnd2arr

    ! Counter
    integer ii
    ! Number of elements

    ! Calculate array size
    !imax = int(ceiling((xmax-xmin)/dx))

    ! Generate array
    do ii = 1, imax
        bnd2arr(ii) = xmin + ii * dx
    end do

end function


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
function read_array(filename,fmtstr,nn,mm,skiplines)
    implicit none

    ! INPUTS:
    ! filename - path to file to be read
    ! fmtstr - input format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), intent(in) :: filename, fmtstr
    ! nn - Number of data rows in file
    ! mm - number of data columns in file
    integer, intent(in) :: nn, mm
    ! skiplines - optional - number of lines to skip from header
    integer, optional :: skiplines

    ! OUTPUT:
    double precision, dimension(nn,mm) :: read_array
    
    ! BODY:

    ! Row,column counters
    integer ii, jj
    ! File unit number
    integer, parameter :: un = 10
    ! Final format to use
    character(len=256) finfmt
    ! Format string for debugging
    character(len=256) dbgfmt

    ! Generate final format string
    write(finfmt,'(A,I1,A,A)') '(', mm, fmtstr, ')'

    ! Generate debug format string
    write(dbgfmt,'(A,I1,A,A)') '(I3, A, ', mm, fmtstr, ', A)'

    ! Print message
    write(*,*) 'Reading data from "', trim(filename), '"'
    write(*,*) 'Using format "', trim(finfmt), '"'

    ! Open file
    open(unit=un, file=trim(filename), status='old', form='formatted')

    ! Skip lines if desired
    if(present(skiplines)) then
        do ii = 1, skiplines
            ! Read without variable ignores the line
            read(un,*)
        end do
    end if

    ! Loop through lines
    do ii = 1, nn - skiplines
        ! Read one row at a time
        !write(*,*) 'Read format: "', trim(finfmt), '"'
        !read(unit=un, fmt=trim(finfmt)) read_array(ii,:)
        read(unit=un, fmt=trim(finfmt)) (read_array(ii,jj),jj=1,mm)
        !write(*,*) 'Write format: "', trim(dbgfmt), '"'
        write(*,fmt=trim(dbgfmt)) ii + skiplines, ': "', read_array(ii,:), '"'
        !write(*,*) 'Done'
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
