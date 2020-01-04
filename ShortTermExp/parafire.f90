!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                      PARAMETER FILE FOR
!      MODELLING MEDITERRANEAN FOREST FIREA UNDER INCREASED ARIDITY    !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! THIS CODE HAS BEEN DEVELOPED WITHIN THE PAPER:
! BAUDENA ET AL 2019 NEW PHYTOLOGIST doi: 10.1111/nph.16252
! https://doi.org/10.1111/nph.16252
!
! Please refer to the paper for the underlying assumptions and equations.
! For more information on this file, contact its author:
! Dr. Mara Baudena m.baudena@uu.nl
! mainfile: mainfire_shortterm.f90
!++++++++++++++++

module parafire

implicit none

! COLONIZATION AND MORTALITY RATES
real :: c1  = 0.047       ! y^-1
real :: m1  = 1/400.      ! y^-1 0.0025
real :: c2  = 0.053        ! y^-1
real :: m2  = 1/125.      ! y^-1 0.008
real :: c3  = .045         ! y^-1
real :: m3  = 1/50.       ! y^-1 0.02
real :: c4  = .067         ! y^-1
real :: m4  = 1/25.       ! y^-1 0.04
real :: c5  = .11         ! y^-1
real :: m5  = 1/15.       ! y^-1 0.066667
real :: c6  = .22         ! y^-1
real :: m6  = 1/40.       ! y^-1 0.025


double precision  :: b0ic(6,9) =reshape((/0.01, 0.01, 0.01, 0.01, 0.01, 0.01, & ! low all
0.89, 0.01, 0.01, 0.01, 0.01, 0.01, &  !  high Q
0.01, 0.89, 0.01,0.01, 0.01, 0.01, &  !  high P
0.01, 0.1, 0.59, 0.1, 0.1, 0.1,     & ! high R+other shrubs
0.01, 0.1, 0.1, 0.59, 0.1, 0.1,    &  ! high U+other shrubs
0.01,0.1,0.1,0.1, 0.59, 0.1,       &  ! high C+other shrubs
0.01, 0.1, 0.1, 0.1, 0.1,0.59, &  !   high B+shrubs
0.15, 0.15, 0.15, 0.15, 0.15, 0.15, & ! successional mixture
0.01, 0.01, 0.15, 0.15, 0.15, 0.15/),(/6,9/)) ! shrubs +B

real  :: Sbank0(6) = (/0., 0.,50., 50.,50.,0. /)
real  :: SbankF(6) = 0.
double precision  :: seedpropb=50. !Proportion of seed to cover

! FILE NAMES
character (len=12)::  file1='bTime.dat'
character (len=12)::  filep='prova.dat'
character (len=9) ::  fileff='fire1.dat'
character (len=12)::  filealpha='alpha.dat'
character (len=14)::  filel='bLoopdi.dat'
character (len=15)::  filelfdb='bLoopFdbdi.dat'


! OTHER PARAMS
integer, parameter :: maxT   = 100              ! y,  Number of years FROM 1ST FIRE when to stop integration (LONGER THAN integration time)
integer,parameter :: dimfv =maxT*365                    ! DIMENSION OF FIRE VECTOR (ONE A DAY)
! ALSO: d,  Number of years for integration
real,parameter   :: h   = 1               ! d, delta t of integration LEAVE IT TO 1 DAY,DO
real    :: hy=h/365                      ! y, delta t of integration
!real    :: Tvect(4)=(/5, 20,50,100/)! y, periods of fire occurence
real    :: Tvect(4)=(/10,20,30,40/)! y, periods of fire occurence
integer :: ifire  = 2               ! 1 to 4, index of Tvect
!integer :: Tmoist = 20              ! d, periods of colonization perturbation

!real :: l(6)=(/.002, 0.1, .02, 0.05, .143, .4/) ! y^-1 flammability frequencies (before Ramon) -> I INVENTED THESE NUMBERS IT SEEMS! SO NO USE
!real :: l(6)=(/.002, 0.05, .0666667, 0.1, .1, .1/) ! y^-1 flammability frequencies (advice Ramon)
real :: l0(6)=(/.002, 0.05, .0666667, 0.1, .1, .1/) ! y^-1 flammability frequencies (from Ramon's, increased for experiments)
!real :: l(6)
!real :: l(6)=(/.002, 0.04, .0666667, 0.1, .1, .2/) ! y^-1 flammability frequencies (advice J&V)
real :: l(6)=(/.00, 0., .0, 0.0, .0, .0/) ! y^-1 zero flammability frequencies -> no fire
real :: flammincr(5)=(/.0, .2, .5, 2., 4./)       ! flammability increase because of drought

real    :: eps    = .0001             ! y^-1, minimum fire frequency when feedback
real    :: fracR(6)  = (/0.,0.,0.,0.,0.,0.4/)             ! fraction of resprouter surviving a fire
!real    :: M      = 100             ! value of the max seed production 1-2y after fire
real    :: aseed = 1               ! y^-1; param (a)dimensionalising the Seeds
!real    :: amax(6)= (/0., 20., .001, .001,.001 ,0./)      ! value of the max seed production (FOR PINES:1-2y after fire)
real    :: amp    = 0.1             ! tanh width for seed germination
real    :: amp2   = 4.             ! tanh width for pine maturity
!real    :: beta   = 25             ! y, shrub seeds half life
real    :: pmat   = 12              ! y, pine maturity
real    :: pseed  = 2               ! y, end of pine seed germination stimulation after fire
real    :: gp(6)  = (/0., .13, .09, .08, .03, 0./) ! %, seed germination percent after fire
real    :: sgp(6) = (/0., .31, .36, .36, .26, 0./) ! %, seedling survival
real    :: e(6) = (/0., 1., 0.5, 0.1, .1, 0./) ! % 1/{# seed that would guarantee dominance, pine=1}
real    :: ms(6) = (/0., .0, .1, .01, .02, 0./)    ! mean seed decay rate (y^-1) from expert estimation of seed durability
real    :: ss(6) = (/0., 0., .45, .8, .75, 0./)! (/0., 0., .3, .8, .75, 0./)      ! %, seed survival to fire (NB THIS CAN CHANGE  SEE DOC -> VARY THIS PARAM)
real    :: sp(6) = (/0., 1., 48., 18., 242., 0./)  ! seed production in proportion to pine
!real    :: C(6) =  (/.0, .01, .01, .01, .01, .01/)       ! conversion rate from seeds to (%space per y-1)
real    :: C = .014                  ! FROM OPTIMIZATION
!real    :: Cvec(4) = (/.01, .1, 1., 10./)
real    :: minfirerettime = 2       ! y
integer :: idumout ! seed for fire stoch

! FIRE VECTORS
real,dimension(dimfv)    :: firev ! DAYS FROM LAST FIRE EVENT: 0=FIRE, #=DAYS AFTER FIRES
real,dimension(dimfv)    :: deltaf ! # DAYS BETWEEN LAST TWO FIRE EVENTS
integer,dimension(dimfv) :: ifv    !  0=NO FIRE, 1=FIRE


! CONSTANTS
DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589d0

end module parafire
