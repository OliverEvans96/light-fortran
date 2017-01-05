! File Name: test_interp.f90
! Description: Test interp.f90
! Created: Wed Jan 04, 2017 | 06:47pm EST
! Last Modified: Thu Jan 05, 2017 | 02:56pm EST
! Author: Oliver Evans <oge1@zips.uakron.edu>

program test_interp
    use rte_utils

    ! File to read
    character(len=256) filename
    ! Array dimentsions
    integer, parameter :: nn=12, mm=2
    ! Array variable
    double precision, dimension(nn,mm) :: arr

    ! File to read
    filename = '/home/oliver/academic/research/kelp/fortran/data/sin4x_12rows_10.2e.txt'

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
