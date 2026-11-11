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

SUBROUTINE super_resolution(temp, qv, temp_sr, qv_sr, um_grid, ml_grid)

USE interpolation_mod, ONLY: linterpolation
USE sr_ennuf_mod, ONLY: sr_ennuf

IMPLICIT NONE

INTEGER, PARAMETER :: nc = 2
INTEGER, PARAMETER :: nzi = 29
INTEGER, PARAMETER :: nzf = 128


REAL(kind=real_umphys) :: stats(4)

REAL(kind=real_umphys), INTENT(IN)    :: temp(nzi), qv(nzi), um_grid(nzi), ml_grid(nzf)
REAL(kind=real_umphys), INTENT(INOUT) :: temp_sr(nzf), qv_sr(nzf)

REAL(kind=real_umphys) :: temp_nn(nzf), qv_nn(nzf)
REAL(kind=real_umphys) :: x_nn(nc,nzf), y_nn(nc,nzf)
    
INTEGER :: i

!Stats for normalization 

DATA stats(:)/ &
271.25368453,  15.13244937, 0.00368781, 0.00401154/ 

!Normalization

x_nn    = 0.0
y_nn    = 0.0
temp_nn = 0.0
qv_nn   = 0.0
temp_sr = 0.0
qv_sr   = 0.0

CALL linterpolation(temp, temp_nn, um_grid, ml_grid)
CALL linterpolation(qv, qv_nn, um_grid, ml_grid)

 x_nn(1,:) = (temp_nn(:) - stats(1) ) / stats(2) 
 x_nn(2,:) = (qv_nn(:)   - stats(3) ) / stats(4) 

CALL sr_ennuf(x_nn, y_nn)

temp_sr(:) = y_nn(1,:) *  stats(2) + stats(1)
qv_sr(:) = y_nn(2,:) * stats(4) + stats(3)

RETURN

END SUBROUTINE super_resolution

END MODULE super_resolution_mod
