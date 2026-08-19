! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! ML Super Resolution Model

MODULE super_resolution_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SUPER_RESOLUTION_MOD'

CONTAINS

SUBROUTINE super_resolution(temp_in, qv_in, temp_sr, qv_sr, um_grid, ml_grid)

USE interpolation_mod, ONLY: linterpolation
USE sr_ennuf_mod, ONLY: sr_ennuf
USE log_mod,                   ONLY: log_event, log_scratch_space, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO
IMPLICIT NONE

INTEGER, PARAMETER :: nc = 2
INTEGER, PARAMETER :: nzi = 29
INTEGER, PARAMETER :: nzf = 128

INTEGER :: k

REAL(kind=real_umphys) :: stats(4)

REAL(kind=real_umphys), INTENT(IN)    :: temp_in(nzi), qv_in(nzi), um_grid(nzi), ml_grid(nzf)
REAL(kind=real_umphys), INTENT(INOUT) :: temp_sr(nzf), qv_sr(nzf)

REAL(kind=real_umphys) :: temp(nzi), qv(nzi)
REAL(kind=real_umphys) :: temp_nn(nzf), qv_nn(nzf)
REAL(kind=real_umphys) :: x_nn(nc,nzf), y_nn(nc,nzf)
REAL(KIND=real_umphys) :: temp_kgi(70)
REAL(KIND=real_umphys) :: q_kgi(70)

INTEGER :: i

!Stats for normalization 

DATA stats(:)/ &
271.25368453,  15.13244937, 0.00368781, 0.00401154/ 

DATA temp_kgi/ &
 294.25,  293.75,  293.12,  292.62,  292.0, &
 291.25,  290.5,  289.75,  288.75,  288.0, &
 287.5,  287.75,  288.12,  288.12,  287.0, &
 285.62,  284.38,  283.0,  281.5,  280.0, &
 278.62,  277.25,  275.88,  274.5,  273.0, &
 271.25,  269.38,  267.12,  264.88,  263.25, &
 260.62,  257.38,  254.5,  252.12,  249.0, &
 245.12,  241.62,  237.75,  233.88,  229.5, &
 225.12,  220.38,  215.38,  210.12,  204.5, &
 198.12,  195.25,  196.38,  193.62,  193.38, &
 198.0,  201.25,  209.0,  209.88,  210.62, &
 213.12,  217.88,  219.25,  223.62,  228.88, &
 233.38,  239.5,  254.38,  263.0,  260.75, &
 256.25,  257.0,  237.5,  216.62, 192.62/
DATA q_kgi/ &
 0.014196,  0.014127,  0.014056,  0.013968,  0.013858, &
 0.013714,  0.013516,  0.013231,  0.012872,  0.012342, &
 0.010987,  0.0088758,  0.0073789,  0.0066231,  0.0066869, &
 0.0067885,  0.0063642,  0.0060629,  0.0056865,  0.0056368, &
 0.0056474,  0.0052389,  0.0044751,  0.003363,  0.0022261, &
 0.0013707,  0.00086451,  0.00063127,  0.00041568,  0.00013107, &
 0.00012183,  0.00012839,  0.00012106,  0.00010854,  0.00010592, &
 0.00010592,  0.00010246,  9.1314e-05,  7.7009e-05,  6.6757e-05, &
 5.4717e-05,  3.9339e-05,  3.2783e-05,  2.7895e-05,  1.8477e-05, &
 9.2983e-06,  4.8876e-06,  3.159e-06,  2.7418e-06,  2.7418e-06, &
 3.0398e-06,  3.0398e-06,  3.3379e-06,  3.4571e-06,  3.7551e-06, &
 3.9935e-06,  3.9339e-06,  3.9935e-06,  3.9339e-06,  3.8743e-06, &
 3.8743e-06,  3.8743e-06,  3.7551e-06,  3.6359e-06,  3.5763e-06, &
 3.0398e-06,  2.6226e-06,  2.563e-06,  2.5034e-06, 2.5034e-06/

temp(:) = temp_kgi(:29)
qv(:) = q_kgi(:29)

!Normalization

x_nn    = 0.0
y_nn    = 0.0
temp_nn = 0.0
qv_nn   = 0.0
temp_sr = 0.0
qv_sr   = 0.0
DO k = 1, nzi
    call log_event( 'Testing KGI: k,temp, q', LOG_LEVEL_INFO )
    write(log_scratch_space,'(I3,F10.4,G12.6)') k, temp(k), qv(k)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
END DO

CALL linterpolation(temp, temp_nn, um_grid, ml_grid)
CALL linterpolation(qv, qv_nn, um_grid, ml_grid)

DO k = 1, 10
    call log_event( 'Testing KGI after interp: k_sr,temp, q', LOG_LEVEL_INFO )
    write(log_scratch_space,'(I3,F10.4,G12.6)') k, temp_nn(k), qv_nn(k)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
END DO

! All looks good until here

 x_nn(1,:) = (temp_nn(:) - stats(1) ) / stats(2) 
 x_nn(2,:) = (qv_nn(:)   - stats(3) ) / stats(4) 

DO k = 1, 10
    call log_event( 'Testing KGI after norm: k_sr,temp, q', LOG_LEVEL_INFO )
    write(log_scratch_space,'(I3,F10.4,G12.6)') k, x_nn(1,k), x_nn(2,k)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
END DO

CALL sr_ennuf(x_nn, y_nn)

DO k = 1, 10
    call log_event( 'Testing KGI after sr (normed): k_sr,temp, q', LOG_LEVEL_INFO )
    write(log_scratch_space,'(I3,F10.4,G12.6)') k, y_nn(1,k), y_nn(2,k)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
END DO

temp_sr(:) = y_nn(1,:) *  stats(2) + stats(1)
qv_sr(:) = y_nn(2,:) * stats(4) + stats(3)

DO k = 1, 10
    call log_event( 'Testing KGI after super res: k_sr,temp_sr, q_sr', LOG_LEVEL_INFO )
    write(log_scratch_space,'(I3,F10.4,G12.6)') k, temp_sr(k), qv_sr(k)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
END DO

RETURN

END SUBROUTINE super_resolution

END MODULE super_resolution_mod
