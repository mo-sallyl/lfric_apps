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
INTEGER, PARAMETER :: n_um = 29
INTEGER, PARAMETER :: n_ml = 128

CONTAINS

SUBROUTINE linterpolation(y_um, y_ml, x_um, x_ml, l_extrapolate)

IMPLICIT NONE

REAL(kind=real_umphys), INTENT(IN) :: x_um(n_um), y_um(n_um), x_ml(n_ml)
REAL(kind=real_umphys), INTENT(INOUT) :: y_ml(n_ml)
LOGICAL, INTENT(IN), OPTIONAL :: l_extrapolate
REAL(kind=real_umphys) :: m, c
LOGICAL :: extrapolate
INTEGER :: i, idx, up, dn

extrapolate = .FALSE.
IF (PRESENT(l_extrapolate)) extrapolate = l_extrapolate

DO i=1, n_ml
   IF (x_ml(i) <= x_um(1)) THEN
      IF (extrapolate) THEN
         dn = 1
         up = 2
      ELSE
         y_ml(i) = y_um(1)
         CYCLE
      END IF
   ELSE IF (x_ml(i) >= x_um(n_um)) THEN
      IF (extrapolate) THEN
         dn = n_um-1
         up = n_um
      ELSE
         y_ml(i) = y_um(n_um)
         CYCLE
      END IF
   ELSE
        idx = MINLOC(ABS(x_um - x_ml(i)), DIM=1)
        IF (x_um(idx)>x_ml(i)) THEN
            up = idx
            dn = idx-1
        ELSE
            up = idx+1
            dn = idx
        END IF
   END IF
   m = (y_um(up) - y_um(dn)) / (x_um(up) - x_um(dn))
   c = y_um(dn) - m * x_um(dn)
   y_ml(i) = x_ml(i) * m + c
END DO    

RETURN

END SUBROUTINE linterpolation


SUBROUTINE linterpolation_b(y_ml, y_um, x_ml, x_um)

IMPLICIT NONE

REAL(kind=real_umphys), INTENT(IN) :: x_um(n_um), y_ml(n_ml), x_ml(n_ml)
REAL(kind=real_umphys), INTENT(INOUT) :: y_um(n_um)
REAL(kind=real_umphys) :: m, c
    
INTEGER :: i, idx, up, dn

DO i=1, n_um
   IF (x_um(i) <= x_ml(1)) THEN
        y_um(i) = y_ml(1)
   ELSE IF (x_um(i) >= x_ml(n_ml)) THEN
        y_um(i) = y_ml(n_ml)
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

REAL(kind=real_umphys), INTENT(IN) :: x_um(n_um), y_ml(n_ml), x_ml(n_ml)
REAL(kind=real_umphys), INTENT(INOUT) :: y_um(n_um)
REAL(kind=real_umphys) :: dn, up
    
INTEGER :: i, j, n, idx


DO i=1, n_um
   n=0
   IF (i==1) THEN
      dn = 0.0
      up = ( x_um(i) + x_um(i+1) ) / 2
   ELSEIF (i==n_um) THEN
      dn = ( x_um(i) + x_um(i-1) ) / 2
      up = x_ml(n_ml)
   ELSE
      dn = ( x_um(i) + x_um(i-1) ) / 2
      up = ( x_um(i) + x_um(i+1) ) / 2
   END IF
   
   DO j=1, n_ml
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
   ! No fine level inside this layer: use the nearest one
   IF (n == 0) idx = MINLOC(ABS(x_ml - x_um(i)), DIM=1)
   y_um(i) = y_ml(idx)     
END DO

END SUBROUTINE max_interpolation

! -----------------------------------------------------------------------------------
!
! conservative_interpolation: overlap-weighted remapping from ML cell centres
! x_ml to UM cell centres x_um. The coordinates must be strictly ascending.
!
! Cell edges are estimated from centre midpoints, with end edges extrapolated
! from the nearest spacing. The overlap measure is the units of the coordinate:
! e.g. height gives a geometric-thickness mean, while -pressure gives a
! pressure-thickness (mass-weighted) mean. Callers can request a per-layer
! coverage result and must not use partial-overlap means as full-layer means.
!
!  Arguments
!    y_ml  (IN)    : field on the fine (ML) grid [n_ml]
!    y_um  (OUT)   : field on the coarse (UM) grid [n_um]
!    x_ml  (IN)    : fine-grid cell-centre coordinate [n_ml], ascending
!    x_um  (IN)    : coarse-grid cell-centre coordinate [n_um], ascending
!    l_clamp_lower : clamp negative lower edges to zero (for height above
!                    surface; disable for transformed pressure coordinates)
!    full_coverage : (OUT, optional) true when the fine grid covers each
!                    complete coarse-grid cell
!
SUBROUTINE conservative_interpolation(y_ml, y_um, x_ml, x_um, l_clamp_lower, &
                                      full_coverage)

IMPLICIT NONE

REAL(kind=real_umphys), INTENT(IN)  :: x_ml(n_ml), y_ml(n_ml)
REAL(kind=real_umphys), INTENT(IN)  :: x_um(n_um)
REAL(kind=real_umphys), INTENT(OUT) :: y_um(n_um)
LOGICAL, INTENT(IN), OPTIONAL :: l_clamp_lower
LOGICAL, INTENT(OUT), OPTIONAL :: full_coverage(n_um)

! Cell-edge arrays (size = n+1 for n cell centres)
REAL(kind=real_umphys) :: edge_ml(n_ml+1)
REAL(kind=real_umphys) :: edge_um(n_um+1)

REAL(kind=real_umphys) :: lo, hi, overlap, total_weight, target_width
LOGICAL :: clamp_lower
INTEGER :: i, j

! ------------------------------------------------------------------
! Build fine-grid (ML) cell edges as midpoints; extend at boundaries
! ------------------------------------------------------------------
clamp_lower = .TRUE.
IF (PRESENT(l_clamp_lower)) clamp_lower = l_clamp_lower

edge_ml(1) = x_ml(1) - 0.5_real_umphys * (x_ml(2) - x_ml(1))
DO j = 2, n_ml
    edge_ml(j) = 0.5_real_umphys * (x_ml(j-1) + x_ml(j))
END DO
edge_ml(n_ml+1) = x_ml(n_ml) + 0.5_real_umphys * (x_ml(n_ml) - x_ml(n_ml-1))

! Clamp lower edge to zero if it went negative
IF (clamp_lower .AND. edge_ml(1) < 0.0_real_umphys) edge_ml(1) = 0.0_real_umphys

! ------------------------------------------------------------------
! Build coarse-grid (UM) cell edges as midpoints; extend at boundaries
! ------------------------------------------------------------------
edge_um(1) = x_um(1) - 0.5_real_umphys * (x_um(2) - x_um(1))
DO i = 2, n_um
    edge_um(i) = 0.5_real_umphys * (x_um(i-1) + x_um(i))
END DO
edge_um(n_um+1) = x_um(n_um) + 0.5_real_umphys * (x_um(n_um) - x_um(n_um-1))

IF (clamp_lower .AND. edge_um(1) < 0.0_real_umphys) edge_um(1) = 0.0_real_umphys

! ------------------------------------------------------------------
! Conservative remapping: weighted average over overlapping fine cells
! ------------------------------------------------------------------
DO i = 1, n_um
    total_weight = 0.0_real_umphys
    y_um(i)      = 0.0_real_umphys
    target_width = edge_um(i+1) - edge_um(i)

    DO j = 1, n_ml
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
    IF (PRESENT(full_coverage)) THEN
        full_coverage(i) = total_weight >= target_width -                    &
            1.0e-10_real_umphys * MAX(1.0_real_umphys, target_width)
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
SUBROUTINE mean_preserving_adjust(y_ml, y_um, x_ml, x_um, l_positive, l_clamp_lower)

IMPLICIT NONE

INTEGER, PARAMETER :: n_iter  = 3
REAL(kind=real_umphys), PARAMETER :: tiny_val = 1.0e-12_real_umphys

REAL(kind=real_umphys), INTENT(IN OUT) :: y_ml(n_ml)
REAL(kind=real_umphys), INTENT(IN)     :: y_um(n_um)
REAL(kind=real_umphys), INTENT(IN)     :: x_ml(n_ml), x_um(n_um)
LOGICAL,                INTENT(IN)     :: l_positive
LOGICAL,                INTENT(IN), OPTIONAL :: l_clamp_lower

REAL(kind=real_umphys) :: edge_ml(n_ml+1), edge_um(n_um+1)
REAL(kind=real_umphys) :: overlap(n_ml, n_um), dx(n_ml)
REAL(kind=real_umphys) :: y_mean(n_um), c_mult(n_um), c_add(n_um)
REAL(kind=real_umphys) :: f_mult, f_add, covered
LOGICAL :: clamp_lower
INTEGER :: i, j, iter

! Same cell edges as conservative_interpolation
clamp_lower = .TRUE.
IF (PRESENT(l_clamp_lower)) clamp_lower = l_clamp_lower

edge_ml(1) = x_ml(1) - 0.5_real_umphys * (x_ml(2) - x_ml(1))
DO j = 2, n_ml
  edge_ml(j) = 0.5_real_umphys * (x_ml(j-1) + x_ml(j))
END DO
edge_ml(n_ml+1) = x_ml(n_ml) + 0.5_real_umphys * (x_ml(n_ml) - x_ml(n_ml-1))
IF (clamp_lower .AND. edge_ml(1) < 0.0_real_umphys) edge_ml(1) = 0.0_real_umphys

edge_um(1) = x_um(1) - 0.5_real_umphys * (x_um(2) - x_um(1))
DO i = 2, n_um
  edge_um(i) = 0.5_real_umphys * (x_um(i-1) + x_um(i))
END DO
edge_um(n_um+1) = x_um(n_um) + 0.5_real_umphys * (x_um(n_um) - x_um(n_um-1))
IF (clamp_lower .AND. edge_um(1) < 0.0_real_umphys) edge_um(1) = 0.0_real_umphys

DO j = 1, n_ml
  dx(j) = edge_ml(j+1) - edge_ml(j)
  DO i = 1, n_um
    overlap(j,i) = MAX(0.0_real_umphys,                                        &
                       MIN(edge_ml(j+1), edge_um(i+1)) - MAX(edge_ml(j), edge_um(i)))
  END DO
END DO

DO iter = 1, n_iter
  CALL conservative_interpolation(y_ml, y_mean, x_ml, x_um, clamp_lower)

  DO i = 1, n_um
    IF (l_positive .AND. y_mean(i) > tiny_val) THEN
      c_mult(i) = y_um(i) / y_mean(i)
      c_add(i)  = 0.0_real_umphys
    ELSE
      c_mult(i) = 1.0_real_umphys
      c_add(i)  = y_um(i) - y_mean(i)
    END IF
  END DO

  DO j = 1, n_ml
    covered = SUM(overlap(j,:))
    IF (covered <= 0.0_real_umphys .OR. dx(j) <= 0.0_real_umphys) CYCLE
    ! Fine cells partly above the UM column keep factor 1 / offset 0 there
    f_mult = (SUM(overlap(j,:) * c_mult) + (dx(j) - covered)) / dx(j)
    f_add  =  SUM(overlap(j,:) * c_add) / dx(j)
    y_ml(j) = y_ml(j) * f_mult + f_add
    IF (l_positive) y_ml(j) = MAX(y_ml(j), 0.0_real_umphys)
  END DO
END DO

END SUBROUTINE mean_preserving_adjust


END MODULE interpolation_mod
