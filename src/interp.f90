! File Name: interp.f90
! Description: Linearly interpolate 1D data at specified point
! Created: Wed Jan 04, 2017 | 06:24pm EST
! Last Modified: Wed Jan 04, 2017 | 06:46pm EST
! Author: Oliver Evans <oge1@zips.uakron.edu>

! INPUTS:
! x0 - x value at which to interpolate
! xx - ordered x values at which y data is sampled
! yy - corresponding y values to interpolate
! nn - length of data
! OUTPUTS:
! interp - interpolated y value
double precision function interp(x0,xx,yy,nn)
    implicit none
    ! Inputs
    double precision, dimension (n) :: xx,yy    
    double precision x0

    ! Index of lower-adjacent data (xx(i) < x0 < xx(i+1))
    integer ii
    ! Slope of liine between (xx(ii),yy(ii)) and (xx(ii+1),yy(ii+1))
    double precision mm

    ! Determine ii
    do ii = 1, nn
        if xx(ii) > x0
        exit
    end do

    ! Calculate slope
    mm = (yy(ii+1) - yy(i)) / (xx(i+1) - xx(i))

    ! Return interpolated value
    interp = yy(ii) + mm * (x0 - xx(ii))
end function interp

