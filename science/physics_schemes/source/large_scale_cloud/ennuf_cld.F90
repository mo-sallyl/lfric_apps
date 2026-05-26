! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Easy Neural Networks in the Um in Fortran: Cloud Scheme

MODULE ennuf_cld_mod

USE um_types, ONLY: real_umphys

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ENNUF_CLD_MOD'
CONTAINS

SUBROUTINE ennuf_cld(p_theta_levels, temp, qv, qcl, qcf, bcf, cfl, cff,        &
                     topography, sigma_h, landfrac, horiz_scale,               &
                     test_kgi, test_kgo )

USE planet_constants_mod,      ONLY: lcrcp, lsrcp
USE atm_fields_bounds_mod,     ONLY: pdims, tdims
USE crmml_ennuf_mod,           ONLY: crmml_ennuf
USE yomhook,                   ONLY: lhook, dr_hook
USE parkind1,                  ONLY: jprb, jpim

IMPLICIT NONE

LOGICAL, INTENT(IN) :: test_kgi
LOGICAL, INTENT(IN) :: test_kgo

REAL(KIND=real_umphys), INTENT(IN) ::                                          &
 p_theta_levels(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end,                                     &
                pdims%k_start:pdims%k_end)
!    Pressure at all points (Pa)

REAL(KIND=real_umphys), INTENT(IN OUT) ::                                      &
 temp(          tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end),                                    &
!    Temperature (K)
 qv(            tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end),                                    &
!    Vapour content (kg water per kg air)
 qcl(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end),                                    &
!    Liquid content (kg water per kg air)
 qcf(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end)
!    Ice water content (kg water per kg air)

REAL(KIND=real_umphys), INTENT(IN OUT) ::                                      &
   bcf(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!    Total cloud fraction (no units)
   cfl(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end),                                  &
!    Liquid cloud fraction (no units)
   cff(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:tdims%k_end)
!    Ice cloud fraction (no units)

  REAL(KIND=real_umphys), INTENT(IN) ::                                          &
   topography(tdims%i_start:tdims%i_end,                                         &
        tdims%j_start:tdims%j_end),                                        &
   sigma_h(tdims%i_start:tdims%i_end,                                            &
     tdims%j_start:tdims%j_end),                                           &
   landfrac(tdims%i_start:tdims%i_end,                                           &
      tdims%j_start:tdims%j_end),                                          &
   horiz_scale(tdims%i_start:tdims%i_end,                                        &
         tdims%j_start:tdims%j_end)

INTEGER, PARAMETER       :: nn_nz     = 70

! Local arrays for intermediate calculations
REAL(KIND=real_umphys) ::                                                      &    
!  qt(             tdims%i_start:tdims%i_end,                                   &
!                  tdims%j_start:tdims%j_end,                                   &
!                              1:nn_nz),                                        &
!  t_li(           tdims%i_start:tdims%i_end,                                   &
!                  tdims%j_start:tdims%j_end,                                   &
!                              1:nn_nz),                                        &
  cfl_minoverlap( tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:nn_nz),                                        &
  cff_minoverlap( tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:nn_nz),                                        &
  cfl_maxoverlap( tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:nn_nz),                                        &
  cff_maxoverlap( tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,                                   &
                              1:nn_nz)

INTEGER :: i,j,k

REAL(KIND=real_umphys) :: t_normed(1:nn_nz)
REAL(KIND=real_umphys) :: q_normed(1:nn_nz)
REAL(KIND=real_umphys) :: p_normed(1:nn_nz)
REAL(KIND=real_umphys) :: bcf_out(1:nn_nz)
REAL(KIND=real_umphys) :: lwc_out(1:nn_nz)
REAL(KIND=real_umphys) :: iwc_out(1:nn_nz)
REAL(KIND=real_umphys) :: cloud_out(1:3*nn_nz)

REAL(KIND=real_umphys) :: temp_1d(1:nn_nz)
REAL(KIND=real_umphys) :: q_1d(1:nn_nz)
REAL(KIND=real_umphys) :: pressure_1d(1:nn_nz)

REAL(KIND=real_umphys) :: temp_kgi(1:nn_nz)
REAL(KIND=real_umphys) :: q_kgi(1:nn_nz)
REAL(KIND=real_umphys) :: pressure_kgi(1:nn_nz)

REAL(KIND=real_umphys), PARAMETER :: max_temp    = 320.0
REAL(KIND=real_umphys), PARAMETER :: min_temp    = 140.0
REAL(KIND=real_umphys), PARAMETER :: max_qv      = 0.025
REAL(KIND=real_umphys), PARAMETER :: max_pres    = 106000.0
REAL(KIND=real_umphys), PARAMETER :: max_orog    = 4000.0
REAL(KIND=real_umphys), PARAMETER :: max_sd_orog = 1200.0
REAL(KIND=real_umphys), PARAMETER :: eps         = 1.0e-3
REAL(KIND=real_umphys), PARAMETER :: cloud_tol   = 0.001
REAL(KIND=real_umphys), PARAMETER :: max_fraction= 0.99
REAL(KIND=real_umphys), PARAMETER :: my_cond_limit= 1.0e-9
REAL(KIND=real_umphys), PARAMETER :: epsilon     =0.622
INTEGER, PARAMETER       :: n_nodes_0 = 215

REAL(KIND=real_umphys), DIMENSION(n_nodes_0) :: tqp_in
REAL(KIND=4), DIMENSION(1,n_nodes_0) :: nn_input
REAL(KIND=4), DIMENSION(1,3*nn_nz) :: nn_output
REAL(KIND=real_umphys), DIMENSION(n_nodes_0) :: norm_avg
REAL(KIND=real_umphys), DIMENSION(n_nodes_0) :: norm_std
!REAL(KIND=real_umphys), DIMENSION(n_nodes_0) :: x_input
REAL(KIND=real_umphys), PARAMETER            :: overlap = 1.0
REAL(KIND=real_umphys)                       :: tmp
REAL(KIND=real_umphys)                       :: num, den
REAL(KIND=real_umphys)                       :: esatwat, esatice
REAL(KIND=real_umphys)                       :: esat,flag,qsat,rht,cmax

DATA norm_avg/ &
 0.7614816,  0.7606579,  0.7594498,  0.758102,  0.7565116, &
 0.7547233,  0.7526109,  0.7501312,  0.747299,  0.7442085, &
 0.7408941,  0.7373815,  0.7336683,  0.7297046,  0.7253591, &
 0.7205253,  0.7151337,  0.7091362,  0.7025111,  0.6952057, &
 0.6872428,  0.6785923,  0.6692305,  0.6591362,  0.6483177, &
 0.6367705,  0.6245682,  0.6115647,  0.5977788,  0.5832471, &
 0.5680537,  0.5523201,  0.5363379,  0.5204732,  0.5052811, &
 0.4916645,  0.4801852,  0.4707024,  0.4626225,  0.4555017, &
 0.4486077,  0.4413044,  0.4336585,  0.425329,  0.4166661, &
 0.4073882,  0.3983682,  0.391034,  0.3865004,  0.3869121, &
 0.393008,  0.4009321,  0.4083993,  0.4157611,  0.4236864, &
 0.4331081,  0.4454411,  0.4601295,  0.478983,  0.5038472, &
 0.5341998,  0.5731899,  0.6232171,  0.6596172,  0.6632298, &
 0.6453263,  0.6134138,  0.5659575,  0.4504302,  0.2756139, &
 0.2742856,  0.2722686,  0.2706556,  0.2689677,  0.2668896, &
 0.2641276,  0.2602384,  0.2548095,  0.247833,  0.2394832, &
 0.2300841,  0.2197695,  0.2086304,  0.19701,  0.1853512, &
 0.1740018,  0.1629188,  0.1519095,  0.1410297,  0.1302883, &
 0.1197892,  0.1095871,  0.0997841,  0.0905079,  0.0817694, &
 0.0734371,  0.0653391,  0.057343,  0.0497797,  0.0427473, &
 0.0362448,  0.0303516,  0.0250791,  0.0203946,  0.0162843, &
 0.0127151,  0.0097084,  0.0072671,  0.0053378,  0.003861, &
 0.0027609,  0.0019561,  0.00138,  0.0009683,  0.0006712, &
 0.000456,  0.0003075,  0.0002138,  0.0001581,  0.00013, &
 0.0001319,  0.0001368,  0.0001407,  0.0001437,  0.0001456, &
 0.0001471,  0.0001481,  0.0001485,  0.0001482,  0.0001466, &
 0.0001445,  0.0001425,  0.0001404,  0.000137,  0.0001308, &
 0.00012,  0.0001107,  0.0001067,  0.0001007,  0.0001, &
 0.9134005,  0.9097782,  0.9047289,  0.8982714,  0.8904339, &
 0.8812503,  0.8707578,  0.8589991,  0.8460214,  0.831876, &
 0.816619,  0.8003125,  0.7830201,  0.76481,  0.7457525, &
 0.725919,  0.7053807,  0.6842111,  0.6624833,  0.6402699, &
 0.6176459,  0.5946846,  0.5714591,  0.5480434,  0.5245084, &
 0.5009269,  0.4773694,  0.4539031,  0.4305927,  0.4075033, &
 0.3846975,  0.3622378,  0.3401874,  0.3186093,  0.2975693, &
 0.2771351,  0.2573777,  0.2383579,  0.2201171,  0.2026773, &
 0.1860416,  0.1701966,  0.1551201,  0.1407862,  0.1271668, &
 0.1142348,  0.1019686,  0.0903625,  0.0794267,  0.0691935, &
 0.0597045,  0.0509698,  0.0429706,  0.0356925,  0.0291352, &
 0.0233073,  0.0182197,  0.0138746,  0.010259,  0.0073426, &
 0.0050714,  0.0033707,  0.0021519,  0.0013136,  0.0007573, &
 0.0004045,  0.0001955,  8.31e-05,  2.95e-05,  8e-06, &
 0.0785298,  0.0322872,  0.3652388,  0.3435108, 0.9424711/
DATA norm_std/ &
 0.0943203,  0.0931822,  0.0919129,  0.0902416,  0.0884476, &
 0.0865239,  0.0846346,  0.0829287,  0.0815012,  0.080401, &
 0.0795941,  0.0790311,  0.0787034,  0.0785102,  0.0783255, &
 0.0781079,  0.0778279,  0.0775562,  0.0773066,  0.0771211, &
 0.0769778,  0.0768746,  0.0767924,  0.0767912,  0.0769425, &
 0.0773531,  0.078051,  0.078841,  0.0796792,  0.0804267, &
 0.0809845,  0.0811737,  0.080691,  0.0793225,  0.0769029, &
 0.0731081,  0.0677457,  0.0614121,  0.0550089,  0.0488834, &
 0.0440105,  0.0416528,  0.0420784,  0.0448282,  0.0502694, &
 0.0571786,  0.0647698,  0.0721742,  0.0772613,  0.0762677, &
 0.0705187,  0.0662425,  0.0640104,  0.0629902,  0.0630606, &
 0.0637174,  0.0641381,  0.0645398,  0.06376,  0.0615994, &
 0.0601811,  0.0612139,  0.0729507,  0.0829573,  0.0913297, &
 0.083133,  0.0570613,  0.0569694,  0.0427269,  0.0449908, &
 0.247596,  0.2452462,  0.2433998,  0.2415298,  0.2393378, &
 0.2364204,  0.232256,  0.2262954,  0.218681,  0.2099915, &
 0.2010373,  0.1922211,  0.1836305,  0.1752317,  0.1669266, &
 0.1586284,  0.1502049,  0.141669,  0.133133,  0.1246795, &
 0.1165071,  0.1087748,  0.1015813,  0.0951033,  0.0893213, &
 0.0839884,  0.0784867,  0.0719455,  0.0651045,  0.0582879, &
 0.0516446,  0.0452545,  0.0391975,  0.0334502,  0.0280935, &
 0.0232458,  0.0190183,  0.015496,  0.0126663,  0.0104419, &
 0.0086894,  0.0072902,  0.0061638,  0.0052036,  0.0043158, &
 0.0034517,  0.002583,  0.0017952,  0.0010202,  0.0001523, &
 2.65e-05,  2.14e-05,  1.77e-05,  1.53e-05,  1.3e-05, &
 1.06e-05,  8.5e-06,  7.2e-06,  8e-06,  1.09e-05, &
 1.39e-05,  1.49e-05,  1.51e-05,  1.57e-05,  1.64e-05, &
 1.69e-05,  1.47e-05,  1.24e-05,  2e-06,  1.2e-06, &
 0.0713495,  0.0709042,  0.0702829,  0.0694901,  0.0685323, &
 0.067418,  0.0661565,  0.0647579,  0.0632332,  0.0615946, &
 0.0598555,  0.0580304,  0.0561344,  0.0541834,  0.052193, &
 0.0501783,  0.0481539,  0.0461343,  0.0441335,  0.0421649, &
 0.0402419,  0.0383766,  0.0365805,  0.0348643,  0.0332379, &
 0.0317113,  0.0302931,  0.0289876,  0.027795,  0.0267118, &
 0.0257299,  0.0248359,  0.0240097,  0.0232249,  0.02245, &
 0.02165,  0.0207894,  0.019842,  0.0187967,  0.0176554, &
 0.016429,  0.0151363,  0.0137999,  0.0124428,  0.0110859, &
 0.0097482,  0.0084549,  0.0072411,  0.0061523,  0.0052413, &
 0.0045302,  0.0039846,  0.0035447,  0.0031604,  0.0028004, &
 0.0024484,  0.002098,  0.001751,  0.0014149,  0.0010995, &
 0.0008183,  0.0005826,  0.0003979,  0.0002636,  0.0001702, &
 0.0001034,  5.46e-05,  2.29e-05,  7.4e-06,  2.1e-06, &
 0.1745441,  0.0819576,  0.475383,  0.2883232, 0.014763/
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
DATA pressure_kgi/ &
 91002.0,  90691.0,  90255.0,  89696.0,  89017.0, &
 88218.0,  87302.0,  86272.0,  85131.0,  83881.0, &
 82529.0,  81080.0,  79542.0,  77921.0,  76218.0, &
 74438.0,  72582.0,  70659.0,  68674.0,  66632.0, &
 64542.0,  62409.0,  60242.0,  58046.0,  55828.0, &
 53593.0,  51349.0,  49099.0,  46850.0,  44611.0, &
 42387.0,  40178.0,  37991.0,  35836.0,  33716.0, &
 31631.0,  29586.0,  27585.0,  25632.0,  23728.0, &
 21874.0,  20073.0,  18324.0,  16627.0,  14984.0, &
 13392.0,  11865.0,  10432.0,  9090.5,  7835.8, &
 6691.1,  5660.8,  4741.5,  3921.5,  3184.0, &
 2531.9,  1969.6,  1493.1,  1098.0,  782.25, &
 537.62,  355.0,  225.88,  138.62,  80.125, &
 42.625,  20.75,  8.875,  3.125, 0.875/

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ENNUF_CLD'

!- End of Header

! ==Main Block==--------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tqp_in(:)=0.0
cloud_out(:)=0.0

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(i,j,k,temp_1d,q_1d,pressure_1d,tqp_in,nn_input,nn_output,        &
!$OMP         cloud_out,                                                        &
!$OMP         bcf_out,lwc_out,iwc_out,tmp,num,den,esatwat,esatice,esat,flag,qsat,rht,cmax)&
!$OMP SHARED(tdims,test_kgo,test_kgi,                                          &
!$OMP        temp_kgi,q_kgi,pressure_kgi,temp,lcrcp,qcl,lsrcp,qcf,qv,          &
!$OMP        p_theta_levels,topography,sigma_h,landfrac,horiz_scale,           &
!$OMP        bcf,cfl,cff,norm_std,norm_avg,                                    &
!$OMP        cfl_maxoverlap,cfl_minoverlap,                                    &
!$OMP        cff_maxoverlap,cff_minoverlap)

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end

      IF (test_kgo) THEN
        ! Use the hard-wired test inputs defined earlier in subroutine.
        DO k = 1, nn_nz
          temp_1d(k)     = temp_kgi(k)
          q_1d(k)        = q_kgi(k)
          pressure_1d(k) = pressure_kgi(k)
        END DO
        tqp_in( 210 + 1 ) = 0.0
        tqp_in( 210 + 2 ) = 0.0
        tqp_in( 210 + 3 ) = 0.0
        tqp_in( 210 + 4 ) = 1.0
        tqp_in( 210 + 5 ) = 1.0
      ELSE
        ! Use the actual model variables as inputs.
        DO k = 1, nn_nz
          temp_1d(k)     = temp(i,j,k) - (lcrcp * qcl(i,j,k)) - (lsrcp * qcf(i,j,k))
          q_1d(k)        = qv(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
          pressure_1d(k) = p_theta_levels(i,j,k)
        END DO
        tqp_in( 210 + 1 ) = topography(i,j) / max_orog
        tqp_in( 210 + 2 ) = sigma_h(i,j) / max_sd_orog
        tqp_in( 210 + 3 ) = landfrac(i,j)
        tqp_in( 210 + 4 ) = horiz_scale(i,j) / ( 144.0 * 1.0e3 )
        tqp_in( 210 + 5 ) = (p_theta_levels(i,j,1)+(1.0*9.81*topography(i,j)))/ max_pres
      END IF ! test_kgo

      DO k = 1, nn_nz
        ! Set everything to zero for good measure
        !qv(i,j,k) =0.0
        !qcl(i,j,k)=0.0
        !qcf(i,j,k)=0.0
        !bcf(i,j,k)=0.0
        !cfl(i,j,k)=0.0
        !cff(i,j,k)=0.0

        ! First renormalisation and pasting of data into the single input vector.
        tqp_in(       k ) =  (temp_1d(k)-min_temp) / (max_temp-min_temp)
        tqp_in(  70 + k ) =  q_1d(k) / max_qv
        tqp_in( 140 + k ) =  pressure_1d(k) / max_pres
      END DO

      ! Do 2nd set of rescaling 
      DO k = 1, 210
        IF (norm_std(k) > 1.0e-4) THEN
          tqp_in(k) = ( tqp_in(k) - norm_avg(k) ) / norm_std(k)
        ELSE
          tqp_in(k) = tqp_in(k) - norm_avg(k)
        END IF
        ! Cap to +/- 3 standard deviations
        tqp_in(k)=min(max(tqp_in(k),-3.0),3.0)
      END DO

      ! For the pressure input cap to +/- 1 standard deviations
      DO k = 141, 210
        tqp_in(k)=min(max(tqp_in(k),-1.0),1.0)
      END DO

      ! Do 2nd set of rescaling for orog, std_orog, landfrac and dx
      tqp_in( 210 + 1 ) = (tqp_in( 210 + 1 )-0.5)*2.0
      tqp_in( 210 + 2 ) = (tqp_in( 210 + 2 )-0.5)*2.0
      tqp_in( 210 + 3 ) = (tqp_in( 210 + 3 )-0.5)*2.0
      tqp_in( 210 + 4 ) = (tqp_in( 210 + 4 )-0.5)*2.0

      nn_input(1,:) = REAL(tqp_in(:), KIND=4)
      CALL crmml_ennuf(nn_input, nn_output)
      cloud_out(:) = REAL(nn_output(1,:), KIND=real_umphys)

      DO k = 1, nn_nz
        ! Extract the different fields from the 1d neural network output.
        ! Store them in 1d stand-alone vectors first (partially to allow KGO checking).
        bcf_out(k)=cloud_out(       k )
        lwc_out(k)=cloud_out(  70 + k )
        iwc_out(k)=cloud_out( 140 + k )

        ! Set the cloud fraction
        bcf(i,j,k)=bcf_out(k)

        !Eqn 17 and 18 from Huang
        !https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
        num=exp(34.494-(4924.99/((temp_1d(k)-273.16)+237.1)))
        den=((temp_1d(k)-273.16)+105.0)**1.57
        esatwat=num/den
        num=exp(43.494-(6545.8/((temp_1d(k)-273.16)+278)))
        den=((temp_1d(k)-273.16)+868.0)**2.0
        esatice=num/den

        IF (temp_1d(k)<273.16) THEN
          flag=1.0
        ELSE
          flag=0.0
        END IF
        esat=(flag*esatice)+((1.0-flag)*esatwat)
        qsat=epsilon*esat/pressure_1d(k)
        rht=q_1d(k)/qsat

        ! Fairly important sanity check
        !cmax=(2.5*rht)-1.5
        cmax=(2.5*rht)-1.625
        cmax=min(max(cmax,0.0),1.0)
        cmax=cmax**2.0

        bcf(i,j,k)=min(cmax,bcf(i,j,k))
        bcf(i,j,k)=min(max(bcf(i,j,k),0.0),1.0)

        ! lwc_out and iwc_out are actually in-cloud values multiplied by 1000.
        ! So multiply by predicted cloud fraction to get grid-box means and / 1000.
        qcl(i,j,k) = lwc_out(k)*bcf(i,j,k)/1000.0
        qcf(i,j,k) = iwc_out(k)*bcf(i,j,k)/1000.0

        ! Retuning, to account for training data being for LAMs in January.
        ! Rescaling obtained by comparing 80 LAMs-worth of data in global model in Jan
        ! compared whole globe for annual data
        !        bcf(i,j,k) = bcf(i,j,k) * 1.29
        !        qcl(i,j,k) = qcl(i,j,k) * 0.91
        !        qcf(i,j,k) = qcf(i,j,k) * 1.00
        ! Multiplying above can push BCF to be greater than 1.0
        bcf(i,j,k)=min(max(bcf(i,j,k),0.0),1.0)

        ! Pragmatic tuning

        ! Henn et al (2024) remove any cloud fraction less than 6%
        ! https://doi.org/10.1029/2023MS003949
        ! So try 0.03, 0.06, 0.09
        IF (bcf(i,j,k)<0.06) THEN
          bcf(i,j,k) = 0.0
          qcl(i,j,k) = 0.0
          qcf(i,j,k) = 0.0 
        END IF

        ! Thresholding at 0.06, has the effect of introducing a step and reducing 0.06 
        ! to half way between 0.0 and 0.06, i.e. 0.03. A similar reduction can be obtained by raising to 1.25
        !bcf(i,j,k)=bcf(i,j,k)**1.25
        !bcf(i,j,k)=bcf(i,j,k)**1.5
        !bcf(i,j,k)=bcf(i,j,k)**2.0

        !IF (rht<0.5) THEN
        !  bcf(i,j,k)=0.0
        !  qcl(i,j,k)=0.0
        !  qcf(i,j,k)=0.0
        !END IF


        ! Check whether sum of condensates is more than total available humidity
        tmp = ( qcl(i,j,k) + qcf(i,j,k) ) / q_1d(k)
        IF (tmp > max_fraction) THEN
          ! Rescale, keeping same ratio of condensate as liquid and ice.
          qcl(i,j,k) = max_fraction * qcl(i,j,k) / tmp
          qcf(i,j,k) = max_fraction * qcf(i,j,k) / tmp
        END IF

        ! Rescaling back to physical units of kg/kg.
        !qcl(i,j,k) = qcl(i,j,k) * q_1d(k)
        !qcf(i,j,k) = qcf(i,j,k) * q_1d(k)

        ! Ensure that if there is some condensate, there is some cloud, 
        ! and enough to prevent in-cloud condensate being more than 5.0e-3 kg/kg,
        ! but not more than 1.0
        !bcf(i,j,k)=min(max(bcf_out(k),(qcl(i,j,k)+qcf(i,j,k))/1.0e-3),1.0)

        ! Use the bulk cloud fraction as the liquid and ice cloud fraction. 
        ! Note: this implies maximum overlap between the liquid and ice regions.
        cfl_maxoverlap(i,j,k) = bcf(i,j,k)
        cff_maxoverlap(i,j,k) = bcf(i,j,k)

        ! Ensure no cloud if no condensate and
        ! set the liquid and ice cloud fractions to be fractions of the bulk cloud fraction
        ! weighted by the fraction of condensate that is liquid and ice respectively.
        ! So not allowing any mixed-phase cloud fraction, hence minimum overlap.

        ! Do this for liquid first.
        IF (qcl(i,j,k) <= my_cond_limit) THEN
          qcl(i,j,k)=0.0
          cfl_maxoverlap(i,j,k)=0.0
          cfl_minoverlap(i,j,k)=0.0
        ELSE
          bcf(i,j,k)=max(bcf(i,j,k),0.001)
          cfl_minoverlap(i,j,k) = bcf(i,j,k) * qcl(i,j,k) / (qcl(i,j,k)+qcf(i,j,k))
        END IF

        ! Then ice.
        IF (qcf(i,j,k) <= my_cond_limit) THEN
          qcf(i,j,k)=0.0
          cff_maxoverlap(i,j,k)=0.0
          cff_minoverlap(i,j,k)=0.0
        ELSE
          bcf(i,j,k)=max(bcf(i,j,k),0.001)
          cff_minoverlap(i,j,k) = bcf(i,j,k) * qcf(i,j,k) / (qcl(i,j,k)+qcf(i,j,k))
        END IF

        ! Blend the 2 types of overlap to get something in between.
        cfl(i,j,k)=(overlap*cfl_maxoverlap(i,j,k))+((1.0-overlap)*cfl_minoverlap(i,j,k))
        cff(i,j,k)=(overlap*cff_maxoverlap(i,j,k))+((1.0-overlap)*cff_minoverlap(i,j,k))

        ! Take the liquid and ice water contents off the total humidity
        qv(i,j,k)=q_1d(k) -qcl(i,j,k) -qcf(i,j,k)
        ! Account for phase changes
        temp(i,j,k)= temp_1d(k) + (lcrcp * qcl(i,j,k)) + (lsrcp * qcf(i,j,k))

        IF (test_kgo) THEN
          ! Use the outputs straight out of the NN 
          ! without any rescaling or consistency checking.
          bcf(i,j,k) = bcf_out(k)
          qcl(i,j,k) = lwc_out(k)
          qcf(i,j,k) = iwc_out(k)
        END IF

        IF (test_kgi) THEN
          bcf(i,j,k) = temp_kgi(k)
          qcl(i,j,k) = q_kgi(k)
          qcf(i,j,k) = pressure_kgi(k)
        END IF ! test_kgi

      END DO !k

    END DO ! i
  END DO ! j

!$OMP END DO

!$OMP END PARALLEL

! End of the subroutine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ennuf_cld
END MODULE ennuf_cld_mod
