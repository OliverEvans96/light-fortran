! File Name: test_vsf.f90
! Description: Test calc_vsf_array
! Created: Wed Jan 04, 2017 | 06:47pm EST
! Last Modified: Fri Jan 06, 2017 | 07:04pm EST
! Author: Oliver Evans <oge1@zips.uakron.edu>

program test_vsf
    use rte_core
    implicit none

    ! VSF file to read
    character(len=256) vsf_file
    ! File length
    integer, parameter :: data_rows = 55
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
    double precision, dimension(lmax,2) :: out_arr

    write(*,*) 'Begin!'

    ! File to read
    vsf_file = '/home/oliver/academic/research/kelp/fortran/data/vsf/nuc_vsf.txt'
    ! Data format
    fmtstr = 'E13.4'

    ! Calculate phi
    write(*,*) 'bnd2arr'
    phi = bnd2arr(0.D0, 2*pi, lmax)

    ! Calculate VSF
    write(*,*) 'calc_vsf_arr'
    beta = calc_vsf_arr(trim(vsf_file), data_rows, lmax, trim(fmtstr), 1)

    ! Generate output array
    out_arr(:,1) = phi
    out_arr(:,2) = beta

    ! Print result
    call print_array(out_arr, lmax, 2)

end program
