! File Name: kelp2d.f90
! Description: Create RTE input data from kelp superindividuals
! Created: Mon Jan 16, 2017 | 03:26pm EST
! Last Modified: Mon Jan 16, 2017 | 04:48pm EST

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

module kelp
contains

! Convert superindividual areas & populations to grid of occupancy probability
function si2prob(area,ind,maxdepth,maxsi,nlayers,imax,jmax)
    ! INPUTS:
    ! area - 2D array of frond areas of superindividuals
    ! 1st dimension is depth layer (kk) , 2nd is individuals in depth layer (nn)
    double precision, dimension(nlayers, maxsi) :: area
    ! ind - Number of individuals represented by corresponding superindividual area
    ! Same dimensions as area
    integer, dimension(nlayers, maxsi) :: ind
    ! maxdepth - depth of bottom of simulation
    double precision, intent(in) :: maxdepth
    ! nlayers - number of depth layers to divide simulation into
    integer, intent(in) :: nlayers
    ! imax - number of output bins in x dimension
    ! jmax - number of output bins in y dimension
    integer, intent(in) :: imax, jmax

    ! OUTPUT:
    ! si2prob - spatial array of probability of kelp occupancy over x and y grid points
    double precision, dimension(imax,jmax) si2prob

    ! BODY:

    ! First, determine 

