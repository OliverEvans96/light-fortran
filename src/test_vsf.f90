! File Name: test_vsf.f90
! Description: Test calc_vsf_array
! Created: Wed Jan 04, 2017 | 06:47pm EST
! Last Modified: Tue Jan 10, 2017 | 08:06pm EST
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

program test_vsf
    use rte_core
    implicit none

    ! VSF file to read
    character(len=256) vsf_file
    ! File length
    !integer, parameter :: data_rows = 55
    integer, parameter :: data_rows = 55
    ! Desired number of points to interpolate
    integer, parameter :: lmax = 20
    ! Step size
    double precision, parameter :: dphi = pi / lmax
    ! Fata format
    character(len=256) fmtstr
    ! Output VSF array
    double precision, dimension(lmax-1) :: beta, norm_beta
    ! angular array
    double precision, dimension(lmax) :: phi
    ! array with both angle and vsf
    double precision, dimension(lmax-1,2) :: out_arr, norm_arr
    double precision integral
    ! Row counter
    integer ii

    write(*,*) 'Begin!'

    ! File to read
    vsf_file = trim(getbasedir()) // '/data/vsf/nuc_vsf.txt'
    ! Data format
    fmtstr = 'E13.4'

    ! Calculate VSF
    write(*,*) 'calc_vsf_arr'
    !out_arr = read_array(trim(vsf_file),'E10.2',48,2)
    !call print_array(out_arr,48,2)

    beta = calc_vsf_arr(vsf_file, data_rows, lmax, trim(fmtstr), 1)
    norm_beta = beta
    integral = normalize(norm_beta,dphi,lmax-1)
    write(*,*) 'Beta'
    call print_array(beta,lmax-1,1)

    ! Generate phi array
    phi = bnd2arr(0.D0,pi,lmax)

    ! Combine phi and vsf
    do ii = 1, lmax-1
        out_arr(ii,1) = phi(ii+1)
        norm_arr(ii,1) = phi(ii+1)
        out_arr(ii,2) = beta(ii)
        norm_arr(ii,2) = norm_beta(ii)
    end do

    ! Print result
    write(*,*) 'Final result'
    call print_array(out_arr, lmax-1, 2)
    call write_array(out_arr, lmax-1, 2,trim(getbasedir())//'/results/vsf_interp.txt','E13.4')
    write(*,*) 'Normalized result'
    call print_array(norm_arr, lmax-1, 2)
    call write_array(norm_arr, lmax-1, 2,trim(getbasedir())//'/results/vsf_norm.txt','E13.4')
    write(*,'(A,E13.4)') 'Integral = ', integral
    write(*,*) 'Written to file'

end program
