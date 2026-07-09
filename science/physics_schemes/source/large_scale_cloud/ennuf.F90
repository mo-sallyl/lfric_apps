! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE ennuf_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ENNUF_MOD'

CONTAINS


SUBROUTINE conv_1d(  &
! in/out
    x_in, y_out, &
! dimensions of the arrays
    channels_in, &
    channels_out, &
    length_in, &
    length_out, &
! weights and biases
    size_kernel, &
    weights, &
    biases, &
! padding, stride and dilation
    pad_mode, &
    padding, &
    stride, &
    dilation)

IMPLICIT NONE

INTEGER, INTENT(in) :: &
    channels_in &
            !channels of input
    ,channels_out &
            !channels of output
    ,length_in &
            !length of input
    ,length_out &
            !length of output
    ,size_kernel
            !size of the kernel
    
REAL(kind = real_umphys), INTENT(IN) :: &
    x_in(channels_in, length_in)
REAL(kind = real_umphys), INTENT(OUT) :: &
    y_out(channels_out, length_out)
    
REAL(kind = real_umphys), INTENT(IN) :: &
    weights(channels_out, channels_in, size_kernel)
REAL(kind = real_umphys), INTENT(IN) :: &
    biases(channels_out)   
    
INTEGER , INTENT(IN) :: &
    padding &
    ,stride &
    ,dilation 

INTEGER :: &
    length_inter
     
CHARACTER (LEN=7), INTENT(IN) :: &
    pad_mode 
    
REAL(kind = real_umphys), ALLOCATABLE :: &
    inter(:,:)
    
INTEGER :: h,i,j,k,s

IF( length_out /= INT( 1 + ( length_in + 2 * padding - dilation * (size_kernel - 1) - 1) / stride ) ) THEN
    PRINT*, "ERROR: "
    PRINT*, "The dimensions of the output array do not correpond to the expected"
    PRINT*, "Check the values of padding, stride and dilation"
    PRINT*, "Ensure that length_out = 1 + ( length_in + 2 * padding - dilation * (size_kernel - 1) - 1) / stride"
END IF
        

length_inter = length_in + 2 * padding    
ALLOCATE(inter(channels_in, length_inter))
    
SELECT CASE (pad_mode)
    CASE ("none   ")
        inter = x_in
    CASE ("zeros  ")
        inter = 0.0
        DO i=1, channels_in
            inter(i,padding+1:-padding-1) = x_in(i,:)
        END DO
    CASE ("reflect")
        DO i=1, channels_in
            inter(i,padding+1:-padding-1) = x_in(i,:)
            inter(i,:padding) = x_in(i,padding+1:2:-1)
            inter(i,length_inter-padding+1:) = x_in(i,length_in-1:length_in-padding-1:-1)
        END DO
END SELECT
    
DO h=1, channels_out
    DO i=1, channels_in
        s = 0
        DO j=1, length_out
            DO k=1, size_kernel
                y_out(h,j) = y_out(h,j) + weights(h,i,k) * inter(i, s+k)
            END DO
        s = s + stride
        END DO
    END DO
    y_out(h, :) = y_out(h, :) + biases(h)
END DO

IF ( ALLOCATED( inter         ) ) DEALLOCATE ( inter         )
     
END SUBROUTINE conv_1d

!-----------------------------------------------------------------------------------------

SUBROUTINE activation_function(                                                & 
! in/out
        x_in                                                                   &
        , y_out                                                                & 
! dimensions of data
        , channels                                                             & 
       , length                                                                & 
! activation choice
      , activation)

IMPLICIT NONE

REAL(KIND=real_umphys), PARAMETER :: alpha=0.2  !slope for leaky relu (make it optional)

INTEGER, INTENT(IN):: channels,length  

CHARACTER (LEN=10), INTENT(IN) :: activation  !choice of activation function

REAL(KIND=real_umphys), INTENT(IN)::                                           &    
  x_in(channels,length)                                                             
REAL(KIND=real_umphys), INTENT(OUT)::                                          &      
  y_out(channels, length)    
  
INTEGER ::                                                                     &
  i,j      ! loop counter

y_out = 0.0
! Main Calculations


SELECT CASE (activation)

    CASE ("relu      ")
        DO i=1, channels
            DO j=1, length
                y_out(i,j) = max(0.0, x_in(i,j))
            END DO
        END DO
    CASE ("leakyrelu ")
        DO i=1, channels
            DO j=1, length
                y_out(i,j) = max(alpha*x_in(i,j), x_in(i,j))
            END DO
        END DO
    CASE ("sigmoid   ")
        DO i=1, channels
            DO j=1, length
                y_out(i,j) = 1.0 / (1.0 + exp(-x_in(i,j)))	
            END DO
        END DO
    CASE ("tanh      ")
        DO i=1, channels
            DO j=1, length
                y_out(i,j) = tanh(x_in(i,j))
            END DO
        END DO
    CASE ("softmax   ")
        DO i=1, channels
            DO j=1, length
                y_out(i,j) = exp(x_in(i,j)) / sum(x_in)
            END DO
        END DO   
    CASE default
        PRINT*,'ERROR'
        PRINT*,'You have asked for an activation function that is not available.'
        PRINT*,'Please check your spelling or add it as an option.'
        PRINT*,'Currently available functions are:' 
        PRINT*,'ReLU, Leaky ReLU, Sigmoid, Tanh and Softmax'

END SELECT


RETURN
END SUBROUTINE activation_function
    
!--------------------------------------------------------------------------------------------------

SUBROUTINE pixel_shuffle( &
! arrays for input and output 
    x_in, &
    y_out, &
! dimensions of data
    channels_in, &
    channels_out, &
    length_in, &
    length_out, &
! upscale factor 
    up_factor)
    
IMPLICIT NONE

INTEGER, INTENT(IN) :: &
    channels_in &
    ,channels_out &
    ,length_in &
    ,length_out &
    ,up_factor
    
REAL(kind=real_umphys), INTENT(IN) :: &
    x_in(channels_in, length_in)
REAL(kind=real_umphys), INTENT(OUT) :: &
    y_out(channels_out, length_out)
        
REAL(kind=real_umphys) :: &
    inter1(up_factor, channels_out, length_in) &
    ,inter2(channels_out, length_in, up_factor) 
    
INTEGER :: i,j,k
    
    ! -------------------------------------------------
    ! Check for consistency in array dimensions  
    ! -------------------------------------------------
    
IF( channels_out /= INT(channels_in / up_factor) ) THEN
    PRINT*, "ERROR: "
    PRINT*, "The dimensions of the output do not correspond to expected"
    PRINT*, "Ensure that channels_out = channels_in / upscale_factor"
ELSE IF ( length_out /= INT(length_in * up_factor) ) THEN
    PRINT*, "ERROR: "
    PRINT*, "The dimensions of the output do not correspond to expected"
    PRINT*, "Ensure that length_out = length_in * upscale_factor"    
END IF
        
j=1
DO i=1, up_factor
    inter1(i,:,:) = x_in(j:j+channels_out,:)
    j = j + channels_out
END DO
        
DO i=1, channels_out
    DO j=1, up_factor
        inter2(i,:,j) = inter1(j,i,:)
    END DO
END DO
        
DO i=1, channels_out
    j=1
    DO k=1, length_in
        y_out(i,j:j+up_factor) = inter2(i,k,:)
        j = j + up_factor
    END DO
END DO

END SUBROUTINE pixel_shuffle

! -----------------------------------------------------------------------------------------------------

SUBROUTINE skip_connection( &
! in/out
    x1_in, &
    x2_in, &
    y_out, &
! dimensions of data
    channels_in, &
    channels_out, &
    length)
    
IMPLICIT NONE

INTEGER, INTENT(IN) :: &
    channels_in &
    ,channels_out &
    ,length
    
REAL(kind=real_umphys), INTENT(IN) :: &
    x1_in(channels_in, length) &
    ,x2_in(channels_in, length) 
REAL(kind=real_umphys), INTENT(OUT) :: &
    y_out(channels_out, length)
    
INTEGER :: i

IF( channels_out /= INT(channels_in * 2) ) THEN
    PRINT*, "ERROR: "
    PRINT*, "The dimensions of the output do not correspond to expected"
    PRINT*, "Ensure that channels_out = channels_in * 2"
END IF

! Main calculations

DO i=1, channels_in
    y_out(i,:) = x1_in(i,:)
    y_out(i + channels_in,:) = x2_in(i,:)
END DO
    
END SUBROUTINE skip_connection

! -----------------------------------------------------------------------------------

SUBROUTINE addition( &
! in/out
    x1_in, &
    x2_in, &
    y_out, &
! dimensions of data
    channels, &
    length)
    
IMPLICIT NONE

INTEGER, INTENT(IN) :: &
    channels           &
    ,length
    
REAL(kind=real_umphys), INTENT(IN) :: &
    x1_in(channels, length) &
    ,x2_in(channels, length)
    
REAL(kind=real_umphys), INTENT(OUT) :: &
    y_out(channels, length)
    
INTEGER :: i
    
DO i=1, channels
    y_out(i,:) = x1_in(i,:) + x2_in(i,:)
END DO
    
END SUBROUTINE addition

END MODULE ennuf_mod
