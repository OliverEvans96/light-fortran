! File Name: pprte.f90
! Description: Solve plane-parallel raditive transfer equation
! numerically given absorption and scattering coefficients as numerical 
! functions of space
! Created: Wed Jan 04, 2017 | 06:07pm EST
! Last Modified:
! Author: Oliver Evans <oge1@zips.uakron.edu>

! INPUTS:
! betafile - path to two column data file describing VSF
!   1st column: angle (rad), 2nd column: normalized transmission
subroutine pprte(betafile,

