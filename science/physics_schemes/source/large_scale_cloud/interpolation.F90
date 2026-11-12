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

! -----------------------------------------------------------------------------------
!
! conservative_interpolation: layer-average (conservative) remapping from the
! high-resolution ML grid (ml_size=128 levels, cell-centres x_ml) down to the
! coarser UM grid (um_size=29 levels, cell-centres x_um).
!
! For each coarse cell i the layer bounds are the midpoints between adjacent
! coarse-level centres (with the domain floor/top as the outermost bounds).
! Every fine cell j whose centre lies inside that layer contributes its value
! weighted by the fraction of the fine cell that overlaps the coarse layer.
! The fine-cell bounds are likewise defined as midpoints between adjacent fine
! centres.  This guarantees exact conservation of the column integral when both
! grids span the same total depth.
!
!  Arguments
!    y_ml  (IN)    : field on the fine (ML) grid  [ml_size]
!    y_um  (OUT)   : field on the coarse (UM) grid [um_size]
!    x_ml  (IN)    : fine-grid cell-centre heights [ml_size], ascending
!    x_um  (IN)    : coarse-grid cell-centre heights [um_size], ascending
!
SUBROUTINE conservative_interpolation(y_ml, y_um, x_ml, x_um)

IMPLICIT NONE

INTEGER, PARAMETER :: um_size = 29
INTEGER, PARAMETER :: ml_size = 128

REAL(kind=real_umphys), INTENT(IN)  :: x_ml(ml_size), y_ml(ml_size)
REAL(kind=real_umphys), INTENT(IN)  :: x_um(um_size)
REAL(kind=real_umphys), INTENT(OUT) :: y_um(um_size)

! Cell-edge arrays (size = n+1 for n cell centres)
REAL(kind=real_umphys) :: edge_ml(ml_size+1)
REAL(kind=real_umphys) :: edge_um(um_size+1)

REAL(kind=real_umphys) :: lo, hi, overlap, total_weight
INTEGER :: i, j

! ------------------------------------------------------------------
! Build fine-grid (ML) cell edges as midpoints; extend at boundaries
! ------------------------------------------------------------------
edge_ml(1) = x_ml(1) - 0.5_real_umphys * (x_ml(2) - x_ml(1))
DO j = 2, ml_size
    edge_ml(j) = 0.5_real_umphys * (x_ml(j-1) + x_ml(j))
END DO
edge_ml(ml_size+1) = x_ml(ml_size) + 0.5_real_umphys * (x_ml(ml_size) - x_ml(ml_size-1))

! Clamp lower edge to zero if it went negative
IF (edge_ml(1) < 0.0_real_umphys) edge_ml(1) = 0.0_real_umphys

! ------------------------------------------------------------------
! Build coarse-grid (UM) cell edges as midpoints; extend at boundaries
! ------------------------------------------------------------------
edge_um(1) = x_um(1) - 0.5_real_umphys * (x_um(2) - x_um(1))
DO i = 2, um_size
    edge_um(i) = 0.5_real_umphys * (x_um(i-1) + x_um(i))
END DO
edge_um(um_size+1) = x_um(um_size) + 0.5_real_umphys * (x_um(um_size) - x_um(um_size-1))

IF (edge_um(1) < 0.0_real_umphys) edge_um(1) = 0.0_real_umphys

! ------------------------------------------------------------------
! Conservative remapping: weighted average over overlapping fine cells
! ------------------------------------------------------------------
DO i = 1, um_size
    total_weight = 0.0_real_umphys
    y_um(i)      = 0.0_real_umphys

    DO j = 1, ml_size
        ! Overlap between fine cell j and coarse cell i
        lo = MAX(edge_ml(j),   edge_um(i))
        hi = MIN(edge_ml(j+1), edge_um(i+1))

        IF (hi > lo) THEN
            overlap      = hi - lo
            y_um(i)      = y_um(i)      + y_ml(j) * overlap
            total_weight = total_weight + overlap
        END IF
    END DO

    IF (total_weight > 0.0_real_umphys) THEN
        y_um(i) = y_um(i) / total_weight
    END IF
END DO

RETURN

END SUBROUTINE conservative_interpolation

! -----------------------------------------------------------------------------------
!
! mean_preserving_adjust: correct a fine (ML) profile so that its layer means
! (as defined by conservative_interpolation) equal the coarse (UM) values,
! keeping the fine-scale structure.  Additive correction for general fields,
! multiplicative for positive-definite fields (falls back to additive where the
! fine layer mean is ~0).  Fine cells straddling two coarse layers get an
! overlap-weighted correction, so a few iterations are used.
!
SUBROUTINE mean_preserving_adjust(y_ml, y_um, x_ml, x_um, l_positive)

IMPLICIT NONE

INTEGER, PARAMETER :: um_size = 29
INTEGER, PARAMETER :: ml_size = 128
INTEGER, PARAMETER :: n_iter  = 3
REAL(kind=real_umphys), PARAMETER :: tiny_val = 1.0e-12_real_umphys

REAL(kind=real_umphys), INTENT(IN OUT) :: y_ml(ml_size)
REAL(kind=real_umphys), INTENT(IN)     :: y_um(um_size)
REAL(kind=real_umphys), INTENT(IN)     :: x_ml(ml_size), x_um(um_size)
LOGICAL,                INTENT(IN)     :: l_positive

REAL(kind=real_umphys) :: edge_ml(ml_size+1), edge_um(um_size+1)
REAL(kind=real_umphys) :: overlap(ml_size, um_size), dz(ml_size)
REAL(kind=real_umphys) :: y_mean(um_size), c_mult(um_size), c_add(um_size)
REAL(kind=real_umphys) :: f_mult, f_add, covered
INTEGER :: i, j, iter

! Same cell edges as conservative_interpolation
edge_ml(1) = x_ml(1) - 0.5_real_umphys * (x_ml(2) - x_ml(1))
DO j = 2, ml_size
  edge_ml(j) = 0.5_real_umphys * (x_ml(j-1) + x_ml(j))
END DO
edge_ml(ml_size+1) = x_ml(ml_size) + 0.5_real_umphys * (x_ml(ml_size) - x_ml(ml_size-1))
IF (edge_ml(1) < 0.0_real_umphys) edge_ml(1) = 0.0_real_umphys

edge_um(1) = x_um(1) - 0.5_real_umphys * (x_um(2) - x_um(1))
DO i = 2, um_size
  edge_um(i) = 0.5_real_umphys * (x_um(i-1) + x_um(i))
END DO
edge_um(um_size+1) = x_um(um_size) + 0.5_real_umphys * (x_um(um_size) - x_um(um_size-1))
IF (edge_um(1) < 0.0_real_umphys) edge_um(1) = 0.0_real_umphys

DO j = 1, ml_size
  dz(j) = edge_ml(j+1) - edge_ml(j)
  DO i = 1, um_size
    overlap(j,i) = MAX(0.0_real_umphys,                                        &
                       MIN(edge_ml(j+1), edge_um(i+1)) - MAX(edge_ml(j), edge_um(i)))
  END DO
END DO

DO iter = 1, n_iter
  CALL conservative_interpolation(y_ml, y_mean, x_ml, x_um)

  DO i = 1, um_size
    IF (l_positive .AND. y_mean(i) > tiny_val) THEN
      c_mult(i) = y_um(i) / y_mean(i)
      c_add(i)  = 0.0_real_umphys
    ELSE
      c_mult(i) = 1.0_real_umphys
      c_add(i)  = y_um(i) - y_mean(i)
    END IF
  END DO

  DO j = 1, ml_size
    covered = SUM(overlap(j,:))
    IF (covered <= 0.0_real_umphys .OR. dz(j) <= 0.0_real_umphys) CYCLE
    ! Fine cells partly above the UM column keep factor 1 / offset 0 there
    f_mult = (SUM(overlap(j,:) * c_mult) + (dz(j) - covered)) / dz(j)
    f_add  =  SUM(overlap(j,:) * c_add) / dz(j)
    y_ml(j) = y_ml(j) * f_mult + f_add
    IF (l_positive) y_ml(j) = MAX(y_ml(j), 0.0_real_umphys)
  END DO
END DO

END SUBROUTINE mean_preserving_adjust


END MODULE interpolation_mod
