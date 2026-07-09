! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Interpolation methods to use with ML Super Resolution.

MODULE interpolation_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INTERPOLATION_MOD'

CONTAINS

SUBROUTINE linterpolation(y_um, y_ml, x_um, x_ml)

IMPLICIT NONE

INTEGER, PARAMETER :: um_size=29
INTEGER, PARAMETER :: ml_size=128

REAL(kind=real_umphys), INTENT(IN) :: x_um(um_size), y_um(um_size), x_ml(ml_size)
REAL(kind=real_umphys), INTENT(INOUT) :: y_ml(ml_size)
REAL(kind=real_umphys) :: m, c
    
INTEGER :: i, idx, up, dn

DO i=1, ml_size    
   IF (x_ml(i) <= x_um(1)) THEN
        y_ml(i) = y_um(1)
   ELSE IF (x_ml(i) >= x_um(um_size)) THEN
        y_ml(i) = y_um(um_size)
   ELSE
        idx = MINLOC(ABS(x_um - x_ml(i)), DIM=1)
        IF (x_um(idx)>x_ml(i)) THEN
            up = idx
            dn = idx-1
        ELSE
            up = idx+1
            dn = idx
        END IF
        m = (y_um(up) - y_um(dn)) / (x_um(up) - x_um(dn))
        c = y_um(dn) - m * x_um(dn)
        y_ml(i) = x_ml(i) * m + c
    END IF
END DO    

RETURN

END SUBROUTINE linterpolation


SUBROUTINE linterpolation_b(y_ml, y_um, x_ml, x_um)

IMPLICIT NONE

INTEGER, PARAMETER :: um_size=29
INTEGER, PARAMETER :: ml_size=128

REAL(kind=real_umphys), INTENT(IN) :: x_um(um_size), y_ml(ml_size), x_ml(ml_size)
REAL(kind=real_umphys), INTENT(INOUT) :: y_um(um_size)
REAL(kind=real_umphys) :: m, c
    
INTEGER :: i, idx, up, dn

DO i=1, um_size    
   IF (x_um(i) <= x_ml(1)) THEN
        y_um(i) = y_ml(1)
   ELSE IF (x_um(i) >= x_ml(ml_size)) THEN
        y_um(i) = y_ml(ml_size)
   ELSE
        idx = MINLOC(ABS(x_ml - x_um(i)), DIM=1)
        IF (x_ml(idx)>x_um(i)) THEN
            up = idx
            dn = idx-1
        ELSE
            up = idx+1
            dn = idx
        END IF
        m = (y_ml(up) - y_ml(dn)) / (x_ml(up) - x_ml(dn))
        c = y_ml(dn) - m * x_ml(dn)
        y_um(i) = x_um(i) * m + c
    END IF
END DO    

RETURN

END SUBROUTINE linterpolation_b

SUBROUTINE max_interpolation(y_ml, y_um, x_ml, x_um)
  
IMPLICIT NONE

INTEGER, PARAMETER :: um_size=29
INTEGER, PARAMETER :: ml_size=128

REAL(kind=real_umphys), INTENT(IN) :: x_um(um_size), y_ml(ml_size), x_ml(ml_size)
REAL(kind=real_umphys), INTENT(INOUT) :: y_um(um_size)
REAL(kind=real_umphys) :: dn, up
    
INTEGER :: i, j, n, idx


DO i=1, um_size
   n=0
   IF (i==1) THEN
      dn = 0.0
      up = ( x_um(i) + x_um(i+1) ) / 2
   ELSEIF (i==um_size) THEN
      dn = ( x_um(i) + x_um(i-1) ) / 2
      up = x_ml(ml_size)
   ELSE
      dn = ( x_um(i) + x_um(i-1) ) / 2
      up = ( x_um(i) + x_um(i+1) ) / 2
   END IF
   
   DO j=1, ml_size
      IF ( x_ml(j) < up ) THEN
         IF ( x_ml(j) > dn ) THEN
            n = n+1
            IF ( n==1 ) THEN
               idx=j
            ELSE
               IF ( y_ml(j) > y_ml(idx) ) THEN
                  idx=j 
               END IF
            END IF
         END IF
      END IF
   END DO
   y_um(i) = y_ml(idx)     
END DO

END SUBROUTINE max_interpolation

END MODULE interpolation_mod
