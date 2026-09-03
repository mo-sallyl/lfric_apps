! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Interpolation methods to use with ML Super Resolution.

MODULE interpolation_mod

use um_types, only: real_umphys

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
! pchip_interpolation: Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
! from the coarse UM grid (um_size=29) to the fine ML grid (ml_size=128).
!
! The method:
!   1. Estimates a derivative d_k at every knot using a weighted harmonic mean of
!      adjacent chord slopes (Fritsch-Carlson, 1980).  Where adjacent slopes change
!      sign the derivative is set to zero, guaranteeing local monotonicity.
!   2. Evaluates the cubic Hermite polynomial H(t) on each sub-interval:
!        H(t) = y_k*(2t^3-3t^2+1) + d_k*h*(t^3-2t^2+t)
!             + y_{k+1}*(-2t^3+3t^2) + d_{k+1}*h*(t^3-t^2)
!      where t = (x - x_k)/h, h = x_{k+1}-x_k.
!
!  Arguments (coarse -> fine, mirrors linterpolation)
!    y_um  (IN)  : field on coarse (UM) grid  [um_size]
!    y_ml  (OUT) : field on fine  (ML) grid   [ml_size]
!    x_um  (IN)  : coarse grid heights [um_size], strictly ascending
!    x_ml  (IN)  : fine   grid heights [ml_size], strictly ascending
!
SUBROUTINE pchip_interpolation(y_um, y_ml, x_um, x_ml)

IMPLICIT NONE

INTEGER, PARAMETER :: um_size = 29
INTEGER, PARAMETER :: ml_size = 128

REAL(kind=real_umphys), INTENT(IN)  :: x_um(um_size), y_um(um_size)
REAL(kind=real_umphys), INTENT(IN)  :: x_ml(ml_size)
REAL(kind=real_umphys), INTENT(OUT) :: y_ml(ml_size)

REAL(kind=real_umphys) :: h(um_size-1)      ! knot intervals
REAL(kind=real_umphys) :: delta(um_size-1)  ! chord slopes
REAL(kind=real_umphys) :: dk(um_size)       ! PCHIP derivatives at each knot

REAL(kind=real_umphys) :: w1, w2, hk, t, phi0, phi1, psi0, psi1, denom
REAL(kind=real_umphys), PARAMETER :: eps = 1.0e-12_real_umphys
INTEGER :: i, j, lo, hi, mid

! ----------------------------------------------------------------
! Step 1: chord slopes and knot spacings
! ----------------------------------------------------------------
DO i = 1, um_size-1
    h(i)     = x_um(i+1) - x_um(i)
    IF (h(i) > eps) THEN
        delta(i) = (y_um(i+1) - y_um(i)) / h(i)
    ELSE
        delta(i) = 0.0_real_umphys
    END IF
END DO

! ----------------------------------------------------------------
! Step 2: PCHIP derivatives (Fritsch-Carlson weighted harmonic mean)
! ----------------------------------------------------------------

! --- Interior knots ---
DO i = 2, um_size-1
    IF (delta(i-1) * delta(i) <= 0.0_real_umphys) THEN
        ! Local extremum or flat: zero derivative enforces monotonicity
        dk(i) = 0.0_real_umphys
    ELSE
        ! Weighted harmonic mean of adjacent chord slopes
        w1    = 2.0_real_umphys * h(i)   + h(i-1)
        w2    = 2.0_real_umphys * h(i-1) + h(i)
        denom = w1 / delta(i-1) + w2 / delta(i)
        IF (ABS(denom) > eps) THEN
            dk(i) = (w1 + w2) / denom
        ELSE
            dk(i) = 0.0_real_umphys
        END IF
    END IF
END DO

! --- Left endpoint (one-sided three-point formula, then clamp) ---
IF (h(1) + h(2) > eps) THEN
    dk(1) = ((2.0_real_umphys * h(1) + h(2)) * delta(1) - h(1) * delta(2)) &
            / (h(1) + h(2))
ELSE
    dk(1) = 0.0_real_umphys
END IF
IF (dk(1) * delta(1) < 0.0_real_umphys) THEN
    dk(1) = 0.0_real_umphys
ELSE IF (delta(1) * delta(2) < 0.0_real_umphys .AND. &
         ABS(dk(1)) > 3.0_real_umphys * ABS(delta(1))) THEN
    dk(1) = 3.0_real_umphys * delta(1)
END IF

! --- Right endpoint ---
IF (h(um_size-1) + h(um_size-2) > eps) THEN
    dk(um_size) = ((2.0_real_umphys * h(um_size-1) + h(um_size-2)) * delta(um_size-1) &
                   - h(um_size-1) * delta(um_size-2)) &
                  / (h(um_size-1) + h(um_size-2))
ELSE
    dk(um_size) = 0.0_real_umphys
END IF
IF (dk(um_size) * delta(um_size-1) < 0.0_real_umphys) THEN
    dk(um_size) = 0.0_real_umphys
ELSE IF (delta(um_size-1) * delta(um_size-2) < 0.0_real_umphys .AND. &
         ABS(dk(um_size)) > 3.0_real_umphys * ABS(delta(um_size-1))) THEN
    dk(um_size) = 3.0_real_umphys * delta(um_size-1)
END IF

! ----------------------------------------------------------------
! Step 3: evaluate cubic Hermite polynomial at each ML point
! ----------------------------------------------------------------
DO j = 1, ml_size

    ! Flat extrapolation outside the data range
    IF (x_ml(j) <= x_um(1)) THEN
        y_ml(j) = y_um(1)
        CYCLE
    ELSE IF (x_ml(j) >= x_um(um_size)) THEN
        y_ml(j) = y_um(um_size)
        CYCLE
    END IF

    ! Binary search for the bracketing knot interval
    lo = 1
    hi = um_size
    DO WHILE (hi - lo > 1)
        mid = (lo + hi) / 2
        IF (x_um(mid) <= x_ml(j)) THEN
            lo = mid
        ELSE
            hi = mid
        END IF
    END DO

    ! Normalised local coordinate t in [0,1]
    hk  = h(lo)
    IF (hk <= eps) THEN
        y_ml(j) = y_um(lo)
        CYCLE
    END IF
    t   = (x_ml(j) - x_um(lo)) / hk

    ! Cubic Hermite basis functions
    phi0 =  2.0_real_umphys*t**3 - 3.0_real_umphys*t**2 + 1.0_real_umphys
    phi1 = -2.0_real_umphys*t**3 + 3.0_real_umphys*t**2
    psi0 = (t**3 - 2.0_real_umphys*t**2 + t) * hk
    psi1 = (t**3 -                 t**2    ) * hk

    y_ml(j) = phi0 * y_um(lo) + phi1 * y_um(lo+1) &
            + psi0 * dk(lo)   + psi1 * dk(lo+1)

END DO

RETURN

END SUBROUTINE pchip_interpolation

! -----------------------------------------------------------------------------------
!
! pchip_interpolation_b: Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
! from the fine ML grid (ml_size=128) back to the coarse UM grid (um_size=29).
!
! This is the fine->coarse companion to pchip_interpolation and mirrors the
! argument order of linterpolation_b.
!
SUBROUTINE pchip_interpolation_b(y_ml, y_um, x_ml, x_um)

IMPLICIT NONE

INTEGER, PARAMETER :: um_size = 29
INTEGER, PARAMETER :: ml_size = 128

REAL(kind=real_umphys), INTENT(IN)  :: x_um(um_size)
REAL(kind=real_umphys), INTENT(IN)  :: x_ml(ml_size), y_ml(ml_size)
REAL(kind=real_umphys), INTENT(OUT) :: y_um(um_size)

REAL(kind=real_umphys) :: h(ml_size-1)
REAL(kind=real_umphys) :: delta(ml_size-1)
REAL(kind=real_umphys) :: dk(ml_size)

REAL(kind=real_umphys) :: w1, w2, hk, t, phi0, phi1, psi0, psi1, denom
REAL(kind=real_umphys), PARAMETER :: eps = 1.0e-12_real_umphys
INTEGER :: i, j, lo, hi, mid

DO i = 1, ml_size-1
    h(i)     = x_ml(i+1) - x_ml(i)
    IF (h(i) > eps) THEN
        delta(i) = (y_ml(i+1) - y_ml(i)) / h(i)
    ELSE
        delta(i) = 0.0_real_umphys
    END IF
END DO

DO i = 2, ml_size-1
    IF (delta(i-1) * delta(i) <= 0.0_real_umphys) THEN
        dk(i) = 0.0_real_umphys
    ELSE
        w1    = 2.0_real_umphys * h(i)   + h(i-1)
        w2    = 2.0_real_umphys * h(i-1) + h(i)
        denom = w1 / delta(i-1) + w2 / delta(i)
        IF (ABS(denom) > eps) THEN
            dk(i) = (w1 + w2) / denom
        ELSE
            dk(i) = 0.0_real_umphys
        END IF
    END IF
END DO

IF (h(1) + h(2) > eps) THEN
    dk(1) = ((2.0_real_umphys * h(1) + h(2)) * delta(1) - h(1) * delta(2)) &
            / (h(1) + h(2))
ELSE
    dk(1) = 0.0_real_umphys
END IF
IF (dk(1) * delta(1) < 0.0_real_umphys) THEN
    dk(1) = 0.0_real_umphys
ELSE IF (delta(1) * delta(2) < 0.0_real_umphys .AND. &
         ABS(dk(1)) > 3.0_real_umphys * ABS(delta(1))) THEN
    dk(1) = 3.0_real_umphys * delta(1)
END IF

IF (h(ml_size-1) + h(ml_size-2) > eps) THEN
    dk(ml_size) = ((2.0_real_umphys * h(ml_size-1) + h(ml_size-2)) * delta(ml_size-1) &
                   - h(ml_size-1) * delta(ml_size-2)) &
                  / (h(ml_size-1) + h(ml_size-2))
ELSE
    dk(ml_size) = 0.0_real_umphys
END IF
IF (dk(ml_size) * delta(ml_size-1) < 0.0_real_umphys) THEN
    dk(ml_size) = 0.0_real_umphys
ELSE IF (delta(ml_size-1) * delta(ml_size-2) < 0.0_real_umphys .AND. &
         ABS(dk(ml_size)) > 3.0_real_umphys * ABS(delta(ml_size-1))) THEN
    dk(ml_size) = 3.0_real_umphys * delta(ml_size-1)
END IF

DO j = 1, um_size

    IF (x_um(j) <= x_ml(1)) THEN
        y_um(j) = y_ml(1)
        CYCLE
    ELSE IF (x_um(j) >= x_ml(ml_size)) THEN
        y_um(j) = y_ml(ml_size)
        CYCLE
    END IF

    lo = 1
    hi = ml_size
    DO WHILE (hi - lo > 1)
        mid = (lo + hi) / 2
        IF (x_ml(mid) <= x_um(j)) THEN
            lo = mid
        ELSE
            hi = mid
        END IF
    END DO

    hk  = h(lo)
    IF (hk <= eps) THEN
        y_um(j) = y_ml(lo)
        CYCLE
    END IF
    t   = (x_um(j) - x_ml(lo)) / hk

    phi0 =  2.0_real_umphys*t**3 - 3.0_real_umphys*t**2 + 1.0_real_umphys
    phi1 = -2.0_real_umphys*t**3 + 3.0_real_umphys*t**2
    psi0 = (t**3 - 2.0_real_umphys*t**2 + t) * hk
    psi1 = (t**3 -                 t**2    ) * hk

    y_um(j) = phi0 * y_ml(lo) + phi1 * y_ml(lo+1) &
            + psi0 * dk(lo)   + psi1 * dk(lo+1)

END DO

RETURN

END SUBROUTINE pchip_interpolation_b

! -----------------------------------------------------------------------------------
!
! mass_weighted_interpolation: mass-weighted (density-weighted) coarsening from the
! high-resolution ML grid (ml_size=128 levels, cell-centres x_ml) down to the
! coarser UM grid (um_size=29 levels, cell-centres x_um).
!
! Identical to conservative_interpolation in how cell edges and overlap fractions
! are computed, but each fine cell's contribution is additionally weighted by the
! hydrostatic air density approximated as
!
!     rho(z) = exp( -z / H_scale )
!
! with H_scale = 8500 m (standard isothermal-atmosphere scale height).
! The combined weight for fine cell j overlapping coarse cell i is therefore
!
!     w_ij = overlap_ij * exp( -x_ml(j) / H_scale )
!
! and the result is
!
!     y_um(i) = SUM_j [ y_ml(j) * w_ij ] / SUM_j [ w_ij ]
!
! This preserves mass-weighted column means in the same way that the standard
! pressure-layer approach does in NWP post-processing.
!
!  Arguments  (fine -> coarse, mirrors conservative_interpolation)
!    y_ml  (IN)    : field on the fine (ML) grid  [ml_size]
!    y_um  (OUT)   : field on the coarse (UM) grid [um_size]
!    x_ml  (IN)    : fine-grid cell-centre heights [ml_size], ascending
!    x_um  (IN)    : coarse-grid cell-centre heights [um_size], ascending
!
SUBROUTINE mass_weighted_interpolation(y_ml, y_um, x_ml, x_um)

IMPLICIT NONE

INTEGER, PARAMETER :: um_size  = 29
INTEGER, PARAMETER :: ml_size  = 128

! Standard isothermal-atmosphere scale height (m)
REAL(kind=real_umphys), PARAMETER :: H_scale = 8500.0_real_umphys

REAL(kind=real_umphys), INTENT(IN)  :: x_ml(ml_size), y_ml(ml_size)
REAL(kind=real_umphys), INTENT(IN)  :: x_um(um_size)
REAL(kind=real_umphys), INTENT(OUT) :: y_um(um_size)

REAL(kind=real_umphys) :: edge_ml(ml_size+1)
REAL(kind=real_umphys) :: edge_um(um_size+1)

REAL(kind=real_umphys) :: lo, hi, overlap, rho_j, total_weight
INTEGER :: i, j

! ------------------------------------------------------------------
! Build fine-grid (ML) cell edges as midpoints; extend at boundaries
! ------------------------------------------------------------------
edge_ml(1) = x_ml(1) - 0.5_real_umphys * (x_ml(2) - x_ml(1))
DO j = 2, ml_size
    edge_ml(j) = 0.5_real_umphys * (x_ml(j-1) + x_ml(j))
END DO
edge_ml(ml_size+1) = x_ml(ml_size) + 0.5_real_umphys * (x_ml(ml_size) - x_ml(ml_size-1))

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
! Mass-weighted remapping
! ------------------------------------------------------------------
DO i = 1, um_size
    total_weight = 0.0_real_umphys
    y_um(i)      = 0.0_real_umphys

    DO j = 1, ml_size
        lo = MAX(edge_ml(j),   edge_um(i))
        hi = MIN(edge_ml(j+1), edge_um(i+1))

        IF (hi > lo) THEN
            overlap      = hi - lo
            rho_j        = EXP(-x_ml(j) / H_scale)
            y_um(i)      = y_um(i)      + y_ml(j) * overlap * rho_j
            total_weight = total_weight + overlap * rho_j
        END IF
    END DO

    IF (total_weight > 0.0_real_umphys) THEN
        y_um(i) = y_um(i) / total_weight
    END IF
END DO

RETURN

END SUBROUTINE mass_weighted_interpolation

END MODULE interpolation_mod
