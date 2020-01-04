!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                      PROGRAM TO MODEL FIRES
!                IN THE MEDITERRANEAN FORESTS UNDER INCREASED ARIDITY  !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! THIS CODE HAS BEEN DEVELOPED WITHIN THE PAPER:
! BAUDENA ET AL 2019 NEW PHYTOLOGIST doi: 10.1111/nph.16252
! https://doi.org/10.1111/nph.16252
!
! Please refer to the paper for the underlying assumptions and equations.
! For more information on this routine, contact its author:
! Dr. Mara Baudena m.baudena@uu.nl
!++++++++++++++++

!  THE MODEL INCLUDES 6 VEGETATION TYPES (IN ORDER OF COMPETITION)
! 1- QUERCUS: RESPROUTER
! 2- P. HALEPENSIS: SEEDER
! 3- ROSEMARY     |
! 4- ULEX         | SHRUB SEEDERS
! 5- CISTUS       |
! 6- BRACHIPODIUM: GRASS RESPROUTER

! THIS VERSION IS USED TO SIMULATE LONG TERM ECOSYSTEM DYNAMICS (10,000 YR)
! program is looping over params (to obtain runs for Fig. 3)

program mainfire_longterm
use parafire

implicit none

double precision,dimension(6) ::  b,db,f,df,fout,bav,bavf,alpha, Sbank
double precision  :: numd

integer           ::  i,NN,ii

! FOR FDB FILE:
real            ::  firevf=100000,deltaff=100000
integer         ::  iifire=0, idummy!=627
integer         ::  num, numok
real            ::  dummy, ran3

! for loop file:
integer             :: iflamm,initcond, ifrac, id, icol
character (len=1)   :: iflstr,initstr, icolstr
character (len=2)   :: idummystr,ifracstr
character (len=20)  :: filename

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! CODE STARTS


NN=dimfv   ! NUMBER OF INTEGRATION TIMESTEP (SAME AS dimfv BECAUSE h=1)

open (22,file=filelfdb)

do id=10,10             ! stochastic fire loop
do iflamm=1,4           ! increasing flammabilities UP TO x3
do initcond=1,6         ! loop on init cond
do ifrac=60,90,15        ! loop on fracR values =ifrac/100 CHANGES R1
do icol=1,4              ! loop on c1 values
    write(idummystr,'(I2)') id
    write(iflstr,'(I1)') iflamm
    write(initstr,'(I1)') initcond
    write(ifracstr,'(I2)') ifrac
    write(icolstr,'(I1)') icol

    ! WRITES FILE NAME FOR SAVING EACH RUN
    filename='bTime'//idummystr//'_'//iflstr//'_'//initstr//'_'//ifracstr//'_'//icolstr//'.dat'
    open (21,file=filename)

    !     params initialization
    firevf=100000
    deltaff=100000
    iifire=0
    fracR(1)=ifrac/100.
    idummy=id
    l=l0*(1+flammincr(iflamm))
    c1=0.047-(icol-1.)*0.012


    !     variables initialization
    b=b0ic(:,initcond)
    bav=0.*b
    bavf=b0ic(:,initcond)*deltaff/100. ! THIS IS B SUMMED UP BETWEEN FIRES, AND AVERAGED AT FIRES (CALCULATED TO BE USED FOR SEEDERS)

    ! intialize alpha (no fire)
    alpha=0.!alp0
    Sbank=Sbank0
    !  IN PROPORTION TO B0 FOR SHRUB SEEDERS
    Sbank(3:5)=seedpropb*b0ic(3:5,initcond)
    ! seed bank initialisation
    call alphacalc(firevf,deltaff,iifire,bavf,b,Sbank,alpha)
          
    !     temporal loop
    do i=1,(NN/2)
          
          ! integration continous dynamics
          f=b 
          call derivs(i,f,df, alpha)
          call rk4(i,f,df,fout,alpha)
          f=fout
          ! stochastic fire dynamics (with fdb)
          dummy=ran3(idummy)

          numd=1/(l(1)*b(1)+l(2)*b(2)+l(3)*b(3)+l(4)*b(4)+l(5)*b(5)+l(6)*b(6)+eps)   ! years (from about 2 to 500 y return time)
          numok=nint(numd)
          num=nint(dummy*365*numok)              ! days
          if (num.eq.numok*365 .and. firevf>minfirerettime*365.) then
            iifire=1                            ! fire!
            deltaff=firevf
            firevf=0.
            bavf=bavf/deltaff*100. ! it is used by alphacalc and then reset to zero there
          else
            iifire=0
            firevf=firevf+h
            bavf=bavf+b/100.
          end if

          ! vegetation burns:
          call fireocc(f,iifire,fout)
          b=fout
          call alphacalc(firevf,deltaff,iifire,bavf,b,Sbank,alpha)

          if (mod(i,100).eq.0 .or. iifire>0) then
            write (21,'(1i15,6f15.4,1i15)') i,b,iifire
          end if
    end do


    ii=0
    do i=(NN/2)+1,NN
        ii=ii+1
        ! integration continous dynamics
        f=b
        call derivs(i,f,df, alpha)
        call rk4(i,f,df,fout,alpha)
        f=fout

        ! stochastic fire dynamics (with fdb)
        dummy=ran3(idummy)
        numd=1/(l(1)*b(1)+l(2)*b(2)+l(3)*b(3)+l(4)*b(4)+l(5)*b(5)+l(6)*b(6)+eps)   ! years (from about 2 to 500 y return time)
        numok=nint(numd)
        num=nint(dummy*365*numok)              ! days
        if (num.eq.numok*365 .and. firevf>minfirerettime*365.) then
            iifire=1                            ! fire!
            deltaff=firevf
            firevf=0.
           bavf=bavf/deltaff*100. ! it is used by alphacalc and then reset to zero there
        else
            iifire=0
            firevf=firevf+h
            bavf=bavf+b/100.
        end if

        ! vegetation burns:
        call fireocc(f,iifire,fout)
        b=fout
        call alphacalc(firevf,deltaff,iifire,bavf,b,Sbank,alpha)
        if (mod(i,100).eq.0 .or. iifire>0) then
            write (21,'(1i15,6f15.4,1i15)') i,b,iifire
        end if

        bav=bav+b/1000
    end do
    bav=bav/ii*1000
    write (* ,'(1i15,6f10.4,5i8,6f10.4)') i,bav, id, iflamm, initcond,ifrac,icol,l
    write (22 ,'(1i15,6f10.4,5i8,6f10.4)') i,bav, id, iflamm, initcond,ifrac,icol,l
end do
end do
end do
end do
end do

end program mainfire_longterm


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        subroutine rk4(tt,y,dydx,yout,g)
!       THIS SUBROUTING INTEGRATE THE DIFFERENTIAL EQUATIONS USING RUNGE-KUPTA 4 INTEGRATION SCHEME
        use parafire
        
        implicit none
        
        real                          :: h6
        double precision,dimension(6) :: dydx,y,yout,k1,k2,k3,k4,y1,y2,y3, g(6)
        integer                       :: tt
        h6=hy/6.

        k1=dydx
        y1=y+0.5*k1*hy
        call derivs(tt,y1,k2,g)
        y2=y+0.5*k2*hy
        call derivs(tt,y2,k3,g)
        y3=y+k3*hy
        call derivs(tt,y3,k4,g)

        yout=y+h6*(k1+2*k2+2*k3+k4)

        return
        end


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        subroutine derivs(tt,v,dv,a)
!       ROUTINE CALCULATING THE MODEL DIFFERENTIAL EQUATIONS
!       v vegetation cover, dv derivatives
        use parafire
        
        implicit none
        integer               :: tt
        double precision,dimension(6)  :: v,dv
        double precision :: a(6), aa(6) ! ALPHA IS ZERO ALWAYS FOR RESPROUTERS (1&6)

        ! NORMALIZES ALPHA (DENOMINATOR IN EQUATION SUPPLEMENTARY MATERIAL) AND MULTIPLIES BY C
        aa=C*a/(a(2)+a(3)+a(4)+a(5)+.0001)*(1-v(1)-v(2)-v(3)-v(4)-v(5)-v(6))
        
        dv(1)=c1*v(1)*(1-v(1))-m1*v(1)
        dv(2)=(c2*v(2))*(1-v(1)-v(2))-m2*v(2)-c1*v(1)*v(2)+aa(2)
        dv(3)=(c3*v(3))*(1-v(1)-v(2)-v(3))-m3*v(3)-c2*v(2)*v(3)-c1*v(1)*v(3)+aa(3)
        dv(4)=(c4*v(4))*(1-v(1)-v(2)-v(3)-v(4))-m4*v(4)-(c2*v(2)+c1*v(1)+c3*v(3))*v(4)+aa(4)
        dv(5)=(c5*v(5))*(1-v(1)-v(2)-v(3)-v(4)-v(5))-m5*v(5)-(c2*v(2)+c1*v(1)+c3*v(3)+c4*v(4))*v(5)+aa(5)
        dv(6)=c6*v(6)*(1-v(1)-v(2)-v(3)-v(4)-v(5)-v(6))-m6*v(6)-(c2*v(2)+c1*v(1)+c3*v(3)+c4*v(4)+c5*v(5))*v(6)
        
        return
        end
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        subroutine alphacalc(tl,dtl,iifire,vavf,v,sb,a)
!       ROUTINE FOR CALCULATING ALPHAs (NB WITHOUT NORMALIZATION AND C, SO IN FACT
!       ROUTINE IS CALCULATING \gamma*S; C and NORMALIZ. ARE IN THE DERIVS SUBROUTINE)
        use parafire
        
        implicit none
        integer :: iifire ! 1-THERE IS A FIRE TODAY, 0 OTHERW.
        real :: tl,dtl, tpine ! tl=t last fire, dtl= t between last two fires
        double precision :: fsa, fsl
        double precision,dimension(6) :: vavf,v,a, sb


        ! PINES
        ! FUNCTION OF THE TIME FROM LAST FIRE tl: GIVES A PEAK FOR TWO YEARS AND THEN ZERO
        fsa=1/2.-1/2.*tanh((tl/365.-pseed+amp)/amp)
        ! FUNCTION THAT DETERMINES WHETHER THERE WERE MATURE PINES BEFORE LAST FIRE  
        fsl=1/2.*(1+tanh((dtl/365.-pmat)/amp2))

!       TO ASSURE THIS TERM DOESN'T GET NEGATIVE: SB2 IS PROPORTIONAL TO THE TIME BETWEEN LAST 2 FIRES -MINUS- THE MATURITY AGE OF PINE
        tpine=(dtl/365-pmat)*fsl
        sb(2)=iifire*(vavf(2)*sp(2)*aseed*tpine*fsl)+sb(2)*fsa*(1-iifire)

!        ! OTHER SEEDERS: THE SEED BANK IS CALCULATED EACH STEP
!        BETWEEN FIRES IT DECREASES ONLY WITH TIME FROM FIRE
!       AT FIRES IT STEEPLY DECRESES IN PROPORTION TO SURVIVAL RATE AND IT IS INCREASED IN PROPORTION TO THE AVERAGE COVER BEFORE FIRE
!
        SbankF(3)=iifire*ss(3)*(sb(3)+vavf(3)*sp(3)*dtl/365)+(1-iifire)*SbankF(3)
        SbankF(4)=iifire*ss(4)*(sb(4)+vavf(4)*sp(4)*dtl/365)+(1-iifire)*SbankF(4)
        SbankF(5)=iifire*ss(5)*(sb(5)+vavf(5)*sp(5)*dtl/365)+(1-iifire)*SbankF(5)

        sb(3)=SbankF(3)*exp(-3*tl/365*ms(3))
        sb(4)=SbankF(4)*exp(-3*tl/365*ms(4))
        sb(5)=SbankF(5)*exp(-3*tl/365*ms(5))

        ! ALPHA IS CALCULATED FROM THE SEED BANKS KEEPING INTO ACCOUNT
        ! GERMINATION RATE AND SEEDLING SURVIVAL
        ! NB C is included in the derivative subroutine
        !NORMALIZATION IS ALSO CALCULATED IN THE derivs subroutine
        
        a(2)=e(2)*sb(2)*gp(2)*sgp(2)
        a(3)=e(3)*sb(3)*gp(3)*sgp(3)
        a(4)=e(4)*sb(4)*gp(4)*sgp(4)
        a(5)=e(5)*sb(5)*gp(5)*sgp(5)

        ! RESET vavf in case of fire
        vavf=vavf*(1-iifire)+0.*iifire

        return
        end

        
        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      subroutine fireocc(b,iifire,bout)
!       THIS ROUTINE RECALCULATE VEGETATION RIGHT AFTER A FIRE (OTHERWISE IT LEAVES IT UNTOUCHED)
      use parafire
      
      implicit none
      
      double precision,dimension(6)       :: b,bout!,bbef
      integer                             :: iifire


      bout=b
      bout(1)=b(1)*(1-iifire+fracR(1)*iifire) ! RESPROUTER
      bout(2)=b(2)*(1-iifire)!+eps*iifire   ! SEEDER eps part is needed only without alphas
      bout(3)=b(3)*(1-iifire)!+eps*iifire   ! SEEDER
      bout(4)=b(4)*(1-iifire)!+eps*iifire   ! SEEDER
      bout(5)=b(5)*(1-iifire)!+eps*iifire   ! SEEDER      
      bout(6)=b(6)*(1-iifire+fracR(6)*iifire) ! RESPROUTER
      
     
      end subroutine
      


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      FUNCTION RAN3(IDUM)
!       RANDOM NUMBER GENERATOR - FROM NUMERICAL RECIPES IN FORTRAN
      SAVE
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
