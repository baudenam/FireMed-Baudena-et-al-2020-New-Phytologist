!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! parameters for "mediterrenanean fires"                                            !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

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


! INITIAL CONDITIONS
double precision  :: b0ic(6,6) =reshape((/0.0039, 0.01, 0.01, 0.01, 0.01, 0.02, & ! low all
0.89, 0.01, 0.01,0.01, 0.01, 0.02, &  !  high Q
0.0039, 0.1, 0.59, 0.1, 0.1, 0.1,     & ! high R
0.0039, 0.2, 0.1, 0.59, 0.1, 0.1,    &  ! high U
0.0039,0.2,0.1,0.1, 0.59, 0.1,       &  ! high C
0.0039, 0.2, 0.1, 0.1, 0.1,0.59 /),(/6,6/)) !   high B


! SEED BANK INITIAL CONDITIONS
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
integer, parameter :: maxT   = 10000              ! y,  Number of years FROM 1ST FIRE when to stop integration (LONGER THAN integration time)
integer,parameter :: dimfv =maxT*365    ! DIMENSION OF FIRE VECTOR (ONE A DAY)
                                        ! ALSO: d,  Number of years for integration
real,parameter   :: h   = 1             ! d, delta t of integration LEAVE IT TO 1 DAY,DO
real    :: hy=h/365                     ! y, delta t of integration
real    :: Tvect(4)=(/10,20,30,40/)     ! y, periods of fire occurence
integer :: ifire  = 2                   ! 1 to 4, index of Tvect

! flammabilities (y^-1 zero flammability frequencies -> no fire)
real :: l(6)=(/.00, 0., .0, 0.0, .0, .0/)           ! IT GETS RECALCULATED IN THE LOOP SCRIPTS
real :: l0(6)=(/.002, 0.05, .0666667, 0.1, .1, .1/) ! y^-1, BASELINE flammability frequencies ( increased for experiments)
real :: flammincr(5)=(/.0, .2, .5, 2., 4./)       ! flammability increase because of drought

real    :: eps    = .0001             ! y^-1, minimum fire frequency when feedback
real    :: fracR(6)  = (/0.,0.,0.,0.,0.,0.4/)             ! fraction of resprouter surviving a fire

! SEEDBANK AND GERMINATION:
real    :: aseed = 1                ! y^-1; param (a)dimensionalising the Seeds
real    :: amp    = 0.1             ! tanh width for seed germination
real    :: amp2   = 4.              ! tanh width for pine maturity
real    :: pmat   = 12              ! y, pine maturity
real    :: pseed  = 2               ! y, end of pine seed germination stimulation after fire
real    :: gp(6)  = (/0., .13, .09, .08, .03, 0./)  ! %, seed germination percent after fire
real    :: sgp(6) = (/0., .31, .36, .36, .26, 0./)  ! %, seedling survival
real    :: e(6) = (/0., 1., 0.5, 0.1, .1, 0./)      ! % 1/{# seed that would guarantee dominance, pine=1}
real    :: ms(6) = (/0., .0, .1, .01, .02, 0./)     ! mean seed decay rate (y^-1) from expert estimation of seed durability
real    :: ss(6) = (/0., 0., .45, .8, .75, 0./)     ! %, seed survival to fire
real    :: sp(6) = (/0., 1., 48., 18., 242., 0./)   ! seed production in proportion to pine
real    :: C = .014                  ! ! conversion rate from seeds to (%space per y-1), obtained FROM OPTIMIZATION

real    :: minfirerettime = 2       ! y
integer :: idumout ! seed for fire stoch

! FIRE VECTORS (FOR NO FDB RUNS)
real,dimension(dimfv)    :: firev ! DAYS FROM LAST FIRE EVENT: 0=FIRE, #=DAYS AFTER FIRES
real,dimension(dimfv)    :: deltaf ! # DAYS BETWEEN LAST TWO FIRE EVENTS
integer,dimension(dimfv) :: ifv    !  0=NO FIRE, 1=FIRE


! CONSTANTS
DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589d0

end module parafire
