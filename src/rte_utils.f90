! File Name: rte_utils.f90
! Description: Linearly interpolate 1D data at specified point
! Created: Wed Jan 04, 2017 | 06:24pm EST
! Last Modified: Thu Jan 05, 2017 | 10:08am EST
! Author: Oliver Evans <oge1@zips.uakron.edu>

module rte_utils
    contains
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
        double precision, intent(out) :: interp

        ! BODY:

        ! Index of lower-adjacent data (xx(i) < x0 < xx(i+1))
        integer ii
        ! Slope of liine between (xx(ii),yy(ii)) and (xx(ii+1),yy(ii+1))
        double precision mm

        ! Determine ii
        do ii = 1, nn
            if (xx(ii) > x0) then
                exit
            end if
        end do

        ! Calculate slope
        mm = (yy(ii+1) - yy(ii)) / (xx(ii+1) - xx(ii))

        ! Return interpolated value
        interp = yy(ii) + mm * (x0 - xx(ii))
    end function interp

    function read_array(filename,fmtstr,nn,mm)
        ! INPUTS:
        ! filename - path to file to be read
        ! fmtstr - data format
        character(len=50), intent(in) :: filename, fmtstr
        ! nn - Number of data rows in file
        ! mm - number of data columns in file
        integer, intent(in) :: nn, mm

        ! OUTPUT:
        double precision, dimension(nn,mm), intent(out) :: read_array
        
        ! BODY:

        ! Row counter
        integer ii

        ! File unit number
        character, parameter :: un = 10

        ! Open file
        open(unit=un, file=filename, status='old', form='formatted')

        ! Loop through lines
        do ii = 1, nn
            ! Read one row at a time
            read(unit=un, fmt=fmtstr) read_array(ii,:)
        end do
    end function read_array
end module
