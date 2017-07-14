












MODULE module_mp_wsm6
!
   IMPLICIT NONE
!
   REAL, PARAMETER, PUBLIC :: dtcldcr     = 120. ! maximum time step for minor loops
   REAL, PARAMETER, PUBLIC :: n0r = 8.e6         ! intercept parameter rain
   REAL, PARAMETER, PUBLIC :: n0g = 4.e6         ! intercept parameter graupel
   REAL, PARAMETER, PUBLIC :: avtr = 841.9       ! a constant for terminal velocity of rain
   REAL, PARAMETER, PUBLIC :: bvtr = 0.8         ! a constant for terminal velocity of rain
   REAL, PARAMETER, PUBLIC :: r0 = .8e-5         ! 8 microm  in contrast to 10 micro m
   REAL, PARAMETER, PUBLIC :: peaut = .55        ! collection efficiency
   REAL, PARAMETER, PUBLIC :: xncr = 3.e8        ! maritime cloud in contrast to 3.e8 in tc80
   REAL, PARAMETER, PUBLIC :: xmyu = 1.718e-5    ! the dynamic viscosity kgm-1s-1
   REAL, PARAMETER, PUBLIC :: avts = 11.72       ! a constant for terminal velocity of snow
   REAL, PARAMETER, PUBLIC :: bvts = .41         ! a constant for terminal velocity of snow
   REAL, PARAMETER, PUBLIC :: avtg = 330.        ! a constant for terminal velocity of graupel
   REAL, PARAMETER, PUBLIC :: bvtg = 0.8         ! a constant for terminal velocity of graupel
   REAL, PARAMETER, PUBLIC :: deng = 500.        ! density of graupel
   REAL, PARAMETER, PUBLIC :: n0smax =  1.e11    ! maximum n0s (t=-90C unlimited)
   REAL, PARAMETER, PUBLIC :: lamdarmax = 8.e4   ! limited maximum value for slope parameter of rain
   REAL, PARAMETER, PUBLIC :: lamdasmax = 1.e5   ! limited maximum value for slope parameter of snow
   REAL, PARAMETER, PUBLIC :: lamdagmax = 6.e4   ! limited maximum value for slope parameter of graupel
   REAL, PARAMETER, PUBLIC :: dicon = 11.9       ! constant for the cloud-ice diamter
   REAL, PARAMETER, PUBLIC :: dimax = 500.e-6    ! limited maximum value for the cloud-ice diamter
   REAL, PARAMETER, PUBLIC :: n0s = 2.e6         ! temperature dependent intercept parameter snow
   REAL, PARAMETER, PUBLIC :: alpha = .12        ! .122 exponen factor for n0s
   REAL, PARAMETER, PUBLIC :: pfrz1 = 100.       ! constant in Biggs freezing
   REAL, PARAMETER, PUBLIC :: pfrz2 = 0.66       ! constant in Biggs freezing
   REAL, PARAMETER, PUBLIC :: qcrmin = 1.e-9     ! minimun values for qr, qs, and qg
   REAL, PARAMETER, PUBLIC :: eacrc = 1.0        ! Snow/cloud-water collection efficiency
   REAL, PARAMETER, PUBLIC :: dens  =  100.0     ! Density of snow
   REAL, PARAMETER, PUBLIC :: qs0   =  6.e-4     ! threshold amount for aggretion to occur
   REAL, SAVE ::                                      &
             qc0, qck1,bvtr1,bvtr2,bvtr3,bvtr4,g1pbr, &
             g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,    &
             bvtr6,g6pbr,                             &
             precr1,precr2,roqimax,bvts1,             &
             bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,     &
             g5pbso2,pvts,pacrs,precs1,precs2,pidn0r, &
             pidn0s,xlv1,pacrc,pi,                    &
             bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,           &
             g3pbg,g4pbg,g5pbgo2,pvtg,pacrg,          &
             precg1,precg2,pidn0g,                    &
             rslopermax,rslopesmax,rslopegmax,        &
             rsloperbmax,rslopesbmax,rslopegbmax,     &
             rsloper2max,rslopes2max,rslopeg2max,     &
             rsloper3max,rslopes3max,rslopeg3max
CONTAINS
!===================================================================
!
!===================================================================
!
!--------------------------------------------------------------------------
      subroutine slope_wsm6(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  REAL, DIMENSION( its:ite , kts:kte,3) ::                                     &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &                                                 
                                                                      rslope2, &                                                 
                                                                      rslope3, &                                                 
                                                                           vt
  REAL, DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  REAL       ::  lamdar, lamdas, lamdag, x, y, z, supcol
  integer :: i, j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
!
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
      do k = kts, kte
!DIR$ VECTOR ALIGNED
        do i = its, ite
          supcol = t0c-t(i,k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
          if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
          vt(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)
          vt(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)
          if(qrs(i,k,1).le.0.0) vt(i,k,1) = 0.0
          if(qrs(i,k,2).le.0.0) vt(i,k,2) = 0.0
          if(qrs(i,k,3).le.0.0) vt(i,k,3) = 0.0
        enddo
      enddo
  END subroutine slope_wsm6
!-----------------------------------------------------------------------------
      subroutine slope_rain(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   & 
                            vt,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               kts,kte
  REAL, DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &      
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL       ::  lamdar, x, y, z, supcol
  integer :: k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
!
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
!DIR$ VECTOR ALIGNED
      do k = kts, kte
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopermax
            rslopeb(k) = rsloperbmax
            rslope2(k) = rsloper2max
            rslope3(k) = rsloper3max
          else
            rslope(k) = 1./lamdar(qrs(k),den(k))
            rslopeb(k) = rslope(k)**bvtr
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvtr*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_rain
!------------------------------------------------------------------------------
      subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               kts,kte
  REAL, DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL       ::  lamdas, x, y, z, supcol
  integer :: k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
!DIR$ VECTOR ALIGNED
      do k = kts, kte
          supcol = t0c-t(k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopesmax
            rslopeb(k) = rslopesbmax
            rslope2(k) = rslopes2max
            rslope3(k) = rslopes3max
          else
            rslope(k) = 1./lamdas(qrs(k),den(k),n0sfac(k))
            rslopeb(k) = rslope(k)**bvts
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvts*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_snow
!----------------------------------------------------------------------------------
      subroutine slope_graup(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,  &
                            vt,kts,kte)
  IMPLICIT NONE
  INTEGER       :: kts,kte
  REAL, DIMENSION( kts:kte) ::                                                 &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &  
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( kts:kte ) ::                                                &
                                                                       n0sfac
  REAL       ::  lamdag, x, y, z, supcol
  integer :: j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
!
!DIR$ ASSUME_ALIGNED qrs:64,den:64,denfac:64,t:64,rslope:64,rslopeb:64,rslope2:64,rslope3:64,vt:64
!DIR$ VECTOR ALIGNED
      do k = kts, kte
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          if(qrs(k).le.qcrmin)then
            rslope(k) = rslopegmax
            rslopeb(k) = rslopegbmax
            rslope2(k) = rslopeg2max
            rslope3(k) = rslopeg3max
          else
            rslope(k) = 1./lamdag(qrs(k),den(k))
            rslopeb(k) = rslope(k)**bvtg
            rslope2(k) = rslope(k)*rslope(k)
            rslope3(k) = rslope2(k)*rslope(k)
          endif
          vt(k) = pvtg*rslopeb(k)*denfac(k)
          if(qrs(k).le.0.0) vt(k) = 0.0
      enddo
  END subroutine slope_graup
!---------------------------------------------------------------------------------
!-------------------------------------------------------------------
      SUBROUTINE nislfv_rain_plm(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,precip,dt,id,iter)
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
      implicit none
      integer  its,ite,kts,kte,id
      real  dt
      real  dzl(its:ite,kts:kte),wwl(its:ite,kts:kte),rql(its:ite,kts:kte),precip(its:ite)
      real  denl(its:ite,kts:kte),denfacl(its:ite,kts:kte),tkl(its:ite,kts:kte)
!
      integer  i,k,n,m,kk,kb,kt,iter
      real  tl,tl2,qql,dql,qqd
      real  th,th2,qqh,dqh
      real  zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real  allold, allnew, zz, dzamin, cflmax, decfl
      real  dz(kts:kte), ww(kts:kte), qq(kts:kte), wd(kts:kte), wa(kts:kte), was(kts:kte)
      real  den(kts:kte), denfac(kts:kte), tk(kts:kte)
      real  wi(kts:kte+1), zi(kts:kte+1), za(kts:kte+1)
      real  qn(kts:kte), qr(kts:kte),tmp(kts:kte),tmp1(kts:kte),tmp2(kts:kte),tmp3(kts:kte)
      real  dza(kts:kte+1), qa(kts:kte+1), qmi(kts:kte+1), qpi(kts:kte+1)
!
!DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,precip:64
      precip(:) = 0.0
!
      i_loop : do i=its,ite
! -----------------------------------
      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)
! skip for no precipitation for all layers
      allold = 0.0
      do k=kts,kte
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif
!
! compute interface values
      zi(kts)=0.0
      do k=kts,kte
        zi(k+1) = zi(k)+dz(k)
      enddo
!
! save departure wind
      wd(:) = ww(:)
      n=1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(kts) = ww(kts)
      wi(kte+1) = ww(kte)
      do k=kts+1,kte
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(kts) = ww(kts)
      wi(kts+1) = 0.5*(ww(kts+1)+ww(kts))
      do k=kts+2,kte-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(kte) = 0.5*(ww(kte)+ww(kte-1))
      wi(kte+1) = ww(kte)
!
! terminate of top of raingroup
      do k=kts+1,kte
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo
!
! diffusivity of wi
      con1 = 0.05
      do k=kte,kts,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
! compute arrival point
      do k=kts,kte+1
        za(k) = zi(k) - wi(k)*dt
      enddo
!
      do k=kts,kte
        dza(k) = za(k+1)-za(k)
      enddo
      dza(kte+1) = zi(kte+1) - za(kte+1)
!
! computer deformation at arrival point
      do k=kts,kte
        qa(k) = qq(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
      enddo
      qa(kte+1) = 0.0
!     call maxmin(kte-kts+1,1,qa,' arrival points ')
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_rain(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,kts,kte)
        if( n.ge.2 ) wa(kts:kte)=0.5*(wa(kts:kte)+was(kts:kte))
        do k=kts,kte
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k),ww(k),wa(k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        was(:) = wa(:)
        n=n+1
        go to 100
      endif
!
! estimate values at arrival cell interface with monotone
      do k=kts+1,kte
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(kts)=qa(kts)
      qmi(kts)=qa(kts)
      qmi(kte+1)=qa(kte+1)
      qpi(kte+1)=qa(kte+1)
!
! interpolation to regular point
      qn = 0.0
      kb=kts
      kt=kts
      intp : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             if( zi(k).ge.za(kte+1) ) then
               exit intp
             else
               find_kb : do kk=kb,kte
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,kte
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
!
       enddo intp
!
! rain out
      sum_precip: do k=kts,kte
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip
!
! replace the new values
!DIR$ VECTOR ALIGNED
      rql(i,:) = qn(:)
!
! ----------------------------------
      enddo i_loop
!
  END SUBROUTINE nislfv_rain_plm
!-------------------------------------------------------------------
      SUBROUTINE nislfv_rain_plm6(its,ite,kts,kte,denl,denfacl,tkl,dzl,wwl,rql,rql2, precip1, precip2,dt,id,iter)
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
      implicit none
      integer  its,ite,kts,kte,id
      real  dt
      real  dzl(its:ite,kts:kte),wwl(its:ite,kts:kte),rql(its:ite,kts:kte),rql2(its:ite,kts:kte),precip(its:ite),precip1(its:ite),precip2(its:ite)
      real  denl(its:ite,kts:kte),denfacl(its:ite,kts:kte),tkl(its:ite,kts:kte)
!
      integer  i,k,n,m,kk,kb,kt,iter,ist
      real  tl,tl2,qql,dql,qqd
      real  th,th2,qqh,dqh
      real  zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real  allold, allnew, zz, dzamin, cflmax, decfl
      real  dz(kts:kte), ww(kts:kte), qq(kts:kte), qq2(kts:kte), wd(kts:kte), wa(kts:kte), wa2(kts:kte), was(kts:kte)
      real  den(kts:kte), denfac(kts:kte), tk(kts:kte)
      real  wi(kts:kte+1), zi(kts:kte+1), za(kts:kte+1)
      real  qn(kts:kte), qr(kts:kte),qr2(kts:kte),tmp(kts:kte),tmp1(kts:kte),tmp2(kts:kte),tmp3(kts:kte)
      real  dza(kts:kte+1), qa(kts:kte+1), qa2(kts:kte+1),qmi(kts:kte+1), qpi(kts:kte+1)
!
!DIR$ ASSUME_ALIGNED denl:64,denfacl:64,tkl:64,dzl:64,wwl:64,rql:64,rql2:64,precip1:64,precip2:64
      precip(:) = 0.0
      precip1(:) = 0.0
      precip2(:) = 0.0
!
      i_loop : do i=its,ite
! -----------------------------------
      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      qq2(:) = rql2(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)
! skip for no precipitation for all layers
      allold = 0.0
      do k=kts,kte
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif
!
! compute interface values
      zi(kts)=0.0
      do k=kts,kte
        zi(k+1) = zi(k)+dz(k)
      enddo
!
! save departure wind
      wd(:) = ww(:)
      n=1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(kts) = ww(kts)
      wi(kte+1) = ww(kte)
      do k=kts+1,kte
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(kts) = ww(kts)
      wi(kts+1) = 0.5*(ww(kts+1)+ww(kts))
      do k=kts+2,kte-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(kte) = 0.5*(ww(kte)+ww(kte-1))
      wi(kte+1) = ww(kte)
!
! terminate of top of raingroup
      do k=kts+1,kte
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo
!
! diffusivity of wi
      con1 = 0.05
      do k=kte,kts,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
! compute arrival point
      do k=kts,kte+1
        za(k) = zi(k) - wi(k)*dt
      enddo
!
      do k=kts,kte
        dza(k) = za(k+1)-za(k)
      enddo
      dza(kte+1) = zi(kte+1) - za(kte+1)
!
! computer deformation at arrival point
      do k=kts,kte
        qa(k) = qq(k)*dz(k)/dza(k)
        qa2(k) = qq2(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
        qr2(k) = qa2(k)/den(k)
      enddo
      qa(kte+1) = 0.0
      qa2(kte+1) = 0.0
!     call maxmin(kte-kts+1,1,qa,' arrival points ')
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,kts,kte)
        call slope_graup(qr2,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa2,kts,kte)
        do k = kts, kte
          tmp(k) = max((qr(k)+qr2(k)), 1.E-15)
          IF ( tmp(k) .gt. 1.e-15 ) THEN
            wa(k) = (wa(k)*qr(k) + wa2(k)*qr2(k))/tmp(k)
          ELSE
            wa(k) = 0.
          ENDIF
        enddo
        if( n.ge.2 ) wa(kts:kte)=0.5*(wa(kts:kte)+was(kts:kte))
        do k=kts,kte
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k), &
!           ww(k),wa(k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        was(:) = wa(:)
        n=n+1
        go to 100
      endif
      ist_loop : do ist = 1, 2
      if (ist.eq.2) then
       qa(:) = qa2(:)
      endif
!
      precip(i) = 0.
!
! estimate values at arrival cell interface with monotone
      do k=kts+1,kte
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(kts)=qa(kts)
      qmi(kts)=qa(kts)
      qmi(kte+1)=qa(kte+1)
      qpi(kte+1)=qa(kte+1)
!
! interpolation to regular point
      qn = 0.0
      kb=kts
      kt=kts
      intp : do k=kts,kte
             kb=max(kb-1,kts)
             kt=max(kt-1,kts)
! find kb and kt
             if( zi(k).ge.za(kte+1) ) then
               exit intp
             else
               find_kb : do kk=kb,kte
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,kte
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
!
       enddo intp
!
! rain out
      sum_precip: do k=kts,kte
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip
!
! replace the new values
      if(ist.eq.1) then
        rql(i,:) = qn(:)
        precip1(i) = precip(i)
      else
        rql2(i,:) = qn(:)
        precip2(i) = precip(i)
      endif
      enddo ist_loop
!
! ----------------------------------
      enddo i_loop
!
  END SUBROUTINE nislfv_rain_plm6
!---------------------------------------------------------------------------------
!+---+-----------------------------------------------------------------+
  ! Read array from unitno and convert from "no-chunk" to "chunk"
  SUBROUTINE readarray2(arr,arrname,unitno,ips,ipe)
    REAL,             INTENT(OUT) :: arr(:,:)
    CHARACTER(LEN=*), INTENT(IN)  :: arrname
    INTEGER,          INTENT(IN)  :: unitno
    INTEGER,          INTENT(IN)  :: ips,ipe
    REAL, ALLOCATABLE :: tmparr(:)
    INTEGER :: i,j,ij,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    jsize = size(arr,2)
    ALLOCATE(tmparr(ips:ipe))
    read(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ! replicate last column if needed
        ipn = min(ij+i+ips-1,ipe)
        arr(i,j) = tmparr(ipn)
      enddo
    enddo
    DEALLOCATE(tmparr)
  END SUBROUTINE readarray2

!+---+-----------------------------------------------------------------+
  ! Read array from unitno and convert from "no-chunk" to "chunk"
  SUBROUTINE readarray3(arr,arrname,unitno,ips,ipe)
    REAL,             INTENT(OUT) :: arr(:,:,:)
    CHARACTER(LEN=*), INTENT(IN)  :: arrname
    INTEGER,          INTENT(IN)  :: unitno
    INTEGER,          INTENT(IN)  :: ips,ipe
    REAL, ALLOCATABLE :: tmparr(:,:)
    INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    jsize = size(arr,3)
    ALLOCATE(tmparr(ips:ipe,ksize))
    read(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ! replicate last column if needed
        ipn = min(ij+i+ips-1,ipe)
        do k = 1,ksize
          arr(i,k,j) = tmparr(ipn,k)
        enddo
      enddo
    enddo
    DEALLOCATE(tmparr)
  END SUBROUTINE readarray3

!+---+-----------------------------------------------------------------+
  ! Read array from unitno and convert from "no-chunk" to "chunk"
  SUBROUTINE readarray4(arr,arrname,unitno,ips,ipe)
    REAL,             INTENT(OUT) :: arr(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN)  :: arrname
    INTEGER,          INTENT(IN)  :: unitno
    INTEGER,          INTENT(IN)  :: ips,ipe
    REAL, ALLOCATABLE :: tmparr(:,:,:)
    INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn,m,msize
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    msize = size(arr,3)
    jsize = size(arr,4)
    ALLOCATE(tmparr(ips:ipe,ksize,msize))
    read(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do m = 1,msize
        do i = 1,CHUNK
          ! replicate last column if needed
          ipn = min(ij+i+ips-1,ipe)
          do k = 1,ksize
            arr(i,k,m,j) = tmparr(ipn,k,m)
          enddo
        enddo
      enddo
    enddo
    DEALLOCATE(tmparr)
  END SUBROUTINE readarray4

!+---+-----------------------------------------------------------------+
  ! Convert array from "chunk" to "no-chunk" and write to unitno.  
  SUBROUTINE writearray2(arr,arrname,unitno,ips,ipe)
    REAL,             INTENT(IN) :: arr(:,:)
    CHARACTER(LEN=*), INTENT(IN) :: arrname
    INTEGER,          INTENT(IN) :: unitno
    INTEGER,          INTENT(IN) :: ips,ipe
    REAL, ALLOCATABLE :: tmparr(:)
    INTEGER :: i,j,ij,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    jsize = size(arr,2)
    ALLOCATE(tmparr(ips:ipe))
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ipn = ij+i+ips-1
        ! skip any replicated columns
        if ((ips<=ipn).and.(ipn<=ipe)) then
          tmparr(ipn) = arr(i,j)
        endif
      enddo
    enddo
    write(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    DEALLOCATE(tmparr)
  END SUBROUTINE writearray2

!+---+-----------------------------------------------------------------+
  ! Convert array from "chunk" to "no-chunk" and write to unitno.  
  SUBROUTINE writearray3(arr,arrname,unitno,ips,ipe)
    REAL,             INTENT(IN) :: arr(:,:,:)
    CHARACTER(LEN=*), INTENT(IN) :: arrname
    INTEGER,          INTENT(IN) :: unitno
    INTEGER,          INTENT(IN) :: ips,ipe
    REAL, ALLOCATABLE :: tmparr(:,:)
    INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    jsize = size(arr,3)
    ALLOCATE(tmparr(ips:ipe,ksize))
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do i = 1,CHUNK
        ipn = ij+i+ips-1
        ! skip any replicated columns
        if ((ips<=ipn).and.(ipn<=ipe)) then
          do k = 1,ksize
            tmparr(ipn,k) = arr(i,k,j)
          enddo
        endif
      enddo
    enddo
    write(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    DEALLOCATE(tmparr)
  END SUBROUTINE writearray3

!+---+-----------------------------------------------------------------+
  ! Convert array from "chunk" to "no-chunk" and write to unitno.  
  SUBROUTINE writearray4(arr,arrname,unitno,ips,ipe)
    REAL,             INTENT(IN) :: arr(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN) :: arrname
    INTEGER,          INTENT(IN) :: unitno
    INTEGER,          INTENT(IN) :: ips,ipe
    REAL, ALLOCATABLE :: tmparr(:,:,:)
    INTEGER :: i,j,k,ij,ksize,jsize,CHUNK,ipn,m,msize
    CHUNK = size(arr,1)
    ksize = size(arr,2)
    msize = size(arr,3)
    jsize = size(arr,4)
    ALLOCATE(tmparr(ips:ipe,ksize,msize))
    do j = 1,jsize
      ij = (j-1)*CHUNK
      do m = 1,msize
        do i = 1,CHUNK
          ipn = ij+i+ips-1
          ! skip any replicated columns
          if ((ips<=ipn).and.(ipn<=ipe)) then
            do k = 1,ksize
              tmparr(ipn,k,m) = arr(i,k,m,j)
            enddo
          endif
        enddo
      enddo
    enddo
    write(unitno) tmparr
    print *,' max ',trim(arrname),' = ',maxval(tmparr),' at ',maxloc(tmparr)
    print *,' min ',trim(arrname),' = ',minval(tmparr),' at ',minloc(tmparr)
    DEALLOCATE(tmparr)
  END SUBROUTINE writearray4

  SUBROUTINE firstTouch(t, q                                      &   
                   ,qci, qrs, den, p, delz                        &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                                                                  )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!  Mimic memory access patterns of wsm62D() while setting physics arrays 
!  to zero.  
!
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL, DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                               t
  REAL, DIMENSION( its:ite , kts:kte, 2 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qci
  REAL, DIMENSION( its:ite , kts:kte, 3 ),                        &
        INTENT(INOUT) ::                                          &
                                                             qrs
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
  REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL,                  &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv
! LOCAL VAR
  INTEGER :: i, k
!
      do k = kts, kte
        do i = its, ite
          t(i,k) = 0.
          q(i,k) = 0.
          qci(i,k,1) = 0.
          qci(i,k,2) = 0.
          qrs(i,k,1) = 0.
          qrs(i,k,2) = 0.
          qrs(i,k,3) = 0.
          den(i,k) = 0.
          p(i,k) = 0.
          delz(i,k) = 0.
        enddo
      enddo
      do i = its, ite
        rain(i) = 0.
        rainncv(i) = 0.
        sr(i) = 0.
      enddo
      if (PRESENT(snow)) then
        do i = its, ite
          snow(i,lat) = 0.
        enddo
      endif
      if (PRESENT(snowncv)) then
        do i = its, ite
          snowncv(i,lat) = 0.
        enddo
      endif
      if (PRESENT(graupel)) then
        do i = its, ite
          graupel(i,lat) = 0.
        enddo
      endif
      if (PRESENT(graupelncv)) then
        do i = its, ite
          graupelncv(i,lat) = 0.
        enddo
      endif
  END SUBROUTINE firstTouch

END MODULE module_mp_wsm6




  PROGRAM wsm6
!-------------------------------------------------------------------
  USE module_mp_wsm6
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
  INTEGER ::   ids, ide,  jds, jde, kds,kde , &
               ims, ime,  jms, jme, kms,kme , &
               its, ite,  jts, jte, kts,kte
  ! "chunk" indices
  INTEGER ::  iids,iide, jjds,jjde,           &
              iims,iime, jjms,jjme,           &
              iits,iite, jjts,jjte
! timers
!
! $Id: gptl.inc.template,v 1.3 2011-03-28 20:55:19 rosinski Exp $
!
! Author: Jim Rosinski
!
! GPTL header file to be included in user code. Values match
! their counterparts in gptl.h. See that file or man pages
! or web-based documenation for descriptions of each value
!
      integer GPTLsync_mpi
      integer GPTLwall
      integer GPTLcpu
      integer GPTLabort_on_error
      integer GPTLoverhead
      integer GPTLdepthlimit
      integer GPTLverbose
      integer GPTLnarrowprint
      integer GPTLpercent
      integer GPTLpersec
      integer GPTLmultiplex
      integer GPTLdopr_preamble
      integer GPTLdopr_threadsort
      integer GPTLdopr_multparent
      integer GPTLdopr_collision
      integer GPTLdopr_memusage
      integer GPTLprint_method
      integer GPTLtablesize
      integer GPTLmaxthreads

      integer GPTL_IPC
      integer GPTL_CI
      integer GPTL_FPC
      integer GPTL_FPI
      integer GPTL_LSTPI
      integer GPTL_DCMRT
      integer GPTL_LSTPDCM
      integer GPTL_L2MRT
      integer GPTL_LSTPL2M
      integer GPTL_L3MRT

      integer GPTLnanotime
      integer GPTLmpiwtime
      integer GPTLclockgettime
      integer GPTLgettimeofday
      integer GPTLpapitime
      integer GPTLplacebo
      integer GPTLread_real_time

      integer GPTLfirst_parent
      integer GPTLlast_parent
      integer GPTLmost_frequent
      integer GPTLfull_tree

      parameter (GPTLsync_mpi       = 0)
      parameter (GPTLwall           = 1)
      parameter (GPTLcpu            = 2)
      parameter (GPTLabort_on_error = 3)
      parameter (GPTLoverhead       = 4)
      parameter (GPTLdepthlimit     = 5)
      parameter (GPTLverbose        = 6)
      parameter (GPTLnarrowprint    = 7)
      parameter (GPTLpercent        = 9)
      parameter (GPTLpersec         = 10)
      parameter (GPTLmultiplex      = 11)
      parameter (GPTLdopr_preamble  = 12)
      parameter (GPTLdopr_threadsort= 13)
      parameter (GPTLdopr_multparent= 14)
      parameter (GPTLdopr_collision = 15)
      parameter (GPTLdopr_memusage  = 27)
      parameter (GPTLprint_method   = 16)
      parameter (GPTLtablesize      = 50)
      parameter (GPTLmaxthreads     = 51)

      parameter (GPTL_IPC           = 17)
      parameter (GPTL_CI            = 18)
      parameter (GPTL_FPC           = 19)
      parameter (GPTL_FPI           = 20)
      parameter (GPTL_LSTPI         = 21)
      parameter (GPTL_DCMRT         = 22)
      parameter (GPTL_LSTPDCM       = 23)
      parameter (GPTL_L2MRT         = 24)
      parameter (GPTL_LSTPL2M       = 25)
      parameter (GPTL_L3MRT         = 26)

      parameter (GPTLgettimeofday   = 1)
      parameter (GPTLnanotime       = 2)
      parameter (GPTLmpiwtime       = 4)
      parameter (GPTLclockgettime   = 5)
      parameter (GPTLpapitime       = 6)
      parameter (GPTLplacebo        = 7)
      parameter (GPTLread_real_time = 3)

      parameter (GPTLfirst_parent   = 1)
      parameter (GPTLlast_parent    = 2)
      parameter (GPTLmost_frequent  = 3)
      parameter (GPTLfull_tree      = 4)

! Externals

      integer gptlsetoption
      integer gptlinitialize
      integer gptlstart
      integer gptlstart_handle
      integer gptlinit_handle
      integer gptlstop
      integer gptlstop_handle
      integer gptlstamp 
      integer gptlpr
      integer gptlpr_file
      integer gptlpr_summary
      integer gptlpr_summary_file
      integer gptlbarrier
      integer gptlreset 
      integer gptlfinalize
      integer gptlget_memusage
      integer gptlprint_memusage
      integer gptlprint_rusage
      integer gptlenable
      integer gptldisable
      integer gptlsetutr
      integer gptlquery
      integer gptlquerycounters
      integer gptlget_wallclock
      integer gptlget_eventvalue
      integer gptlget_nregions
      integer gptlget_regionname
      integer gptl_papilibraryinit
      integer gptlevent_name_to_code
      integer gptlevent_code_to_name
      integer gptlnum_errors
      integer gptlnum_warn
      integer gptlget_count

      external gptlsetoption
      external gptlinitialize
      external gptlstart
      external gptlstart_handle
      external gptlinit_handle
      external gptlstop
      external gptlstop_handle
      external gptlstamp 
      external gptlpr
      external gptlpr_file
      external gptlpr_summary
      external gptlpr_summary_file
      external gptlbarrier
      external gptlreset 
      external gptlfinalize
      external gptlget_memusage
      external gptlprint_memusage
      external gptlprint_rusage
      external gptlenable
      external gptldisable
      external gptlsetutr
      external gptlquery
      external gptlquerycounters
      external gptlget_wallclock
      external gptlget_eventvalue
      external gptlget_nregions
      external gptlget_regionname
      external gptl_papilibraryinit
      external gptlevent_name_to_code
      external gptlevent_code_to_name
      external gptlnum_errors
      external gptlnum_warn
      external gptlget_count
  INTEGER :: ret
  REAL*8  :: totaltime
  INTEGER :: handle = 0
! Contains definition of _CHUNK_

! Contains definition of _NZ_

  INTEGER :: CHUNK, NZ_BUILD
  INTEGER :: num_tiles_C
  INTEGER ::  itimestep
  REAL, ALLOCATABLE :: t(:,:,:)
  REAL, ALLOCATABLE :: qci(:,:,:,:)
  REAL, ALLOCATABLE :: qrs(:,:,:,:)
  REAL, ALLOCATABLE :: q(:,:,:)
  REAL, ALLOCATABLE :: den(:,:,:)
  REAL, ALLOCATABLE :: p(:,:,:)
  REAL, ALLOCATABLE :: delz(:,:,:)
  REAL :: delt, g, rd, rv, t0c, den0, cpd, cpv, ep1, ep2, qmin, &
          XLS, XLV0, XLF0, cliq, cice, psat, denr
  REAL, ALLOCATABLE :: rain(:,:)
  REAL, ALLOCATABLE :: rainncv(:,:)
  REAL, ALLOCATABLE :: sr(:,:)
  REAL, ALLOCATABLE :: snow(:,:)
  REAL, ALLOCATABLE :: snowncv(:,:)
  REAL, ALLOCATABLE :: graupel(:,:)
  REAL, ALLOCATABLE :: graupelncv(:,:)

!+---+-----------------------------------------------------------------+
! LOCAL VAR
  INTEGER ::               i,j,k
  INTEGER :: ios, unitno
  CHARACTER(LEN=64) :: fn  ! file name

!+---+-----------------------------------------------------------------+
! EXTERNAL FUNCTIONS
  integer, external :: omp_get_max_threads
!+---+-----------------------------------------------------------------+

  call gptlprocess_namelist ('GPTLnamelist', 77, ret)
  ret = gptlinitialize ()
  ret = gptlstart('Total')

  ! read constants
  PRINT *,'wsm6init():  read constants'
  fn = 'wsm6_constants.dat'
  unitno=31
  open (unitno,file=trim(fn),form="unformatted",action='read', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open constant file ',trim(fn), &
               ' . stopping'
    stop
  endif
  read(unitno) pi, xlv1
  read(unitno) qc0, qck1
  read(unitno) bvtr1, bvtr2, bvtr3, bvtr4, bvtr6, g1pbr, g3pbr, &
               g4pbr, g6pbr, g5pbro2, pvtr, eacrr, pacrr, &
               precr1, precr2, roqimax
  read(unitno) bvts1, bvts2, bvts3, bvts4, g1pbs, g3pbs, g4pbs,  &
               g5pbso2, pvts, pacrs, precs1, precs2, pidn0r, pidn0s
  read(unitno) pacrc
  read(unitno) bvtg1, bvtg2, bvtg3, bvtg4, g1pbg, g3pbg, g4pbg,  &
               pacrg, g5pbgo2, pvtg, precg1, precg2, pidn0g
  read(unitno) rslopermax, rslopesmax, rslopegmax, rsloperbmax,  &
               rslopesbmax, rslopegbmax, rsloper2max, rslopes2max,  &
               rslopeg2max, rsloper3max, rslopes3max, rslopeg3max
  close(unitno)

  ! read input data
  PRINT *,'wsm62D():  read input state'
  fn = 'wsm6_input.dat'
  unitno=31
  open (unitno,file=trim(fn),form="unformatted",action='read', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open input file ',trim(fn), &
               ' . stopping'
    stop
  endif
  read(unitno) itimestep
  PRINT *,'wsm62D():  itimestep == ',itimestep
  ! read serial versions of indices
  read(unitno) ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte

  ! assert constraints on indices
  ! NOTE that [ikj]d[se] are ignored
  if ((ims/=its).or.(ime/=ite).or.(jms/=jts).or.(jme/=jte)) then
    print *,'ERROR:  index mismatch reading file ',trim(fn), &
            ' . stopping'
    stop
  endif
  if ((ims/=1).or.(jms/=1)) then
    print *,'ERROR:  incorrect start index reading file ',trim(fn), &
            ' . stopping'
    stop
  endif

  ! set default values of "chunk" indices
  iids = ids
  iide = ide
  iims = ims
  iime = ime
  iits = its
  iite = ite
  jjds = jds
  jjde = jde
  jjms = jms
  jjme = jme
  jjts = jts
  jjte = jte

  ! set up optional "i" chunking and optional fixed vertical extent
   CHUNK = 4
   NZ_BUILD = 32
   ! if specified, NZ_BUILD must match namelist setting for nvl
   if (NZ_BUILD/=kte-kts+1) then
     print *, 'ERROR:  Build-time-specified NZ must equal namelist nz, values are: ',NZ_BUILD,kte-kts+1
     call flush(6)
     stop
   endif
   num_tiles_C = (iite-iits+1) / CHUNK
   if (mod((iite-iits+1),CHUNK) > 0) then
     num_tiles_C = num_tiles_C + 1
   endif
   iime = CHUNK
   iite = CHUNK
   jjme = num_tiles_C
   jjte = num_tiles_C
   PRINT *,'ims,ime,iims,iime',ims,ime,iims,iime
   PRINT "('omp_get_max_threads() returned:  ',I9)",omp_get_max_threads()
   PRINT *,'CHUNK =', CHUNK
   PRINT *,'NUMBER OF CHUNKS =', num_tiles_C
   PRINT *,'NZ_BUILD =', NZ_BUILD

  ! allocate arrays
  ALLOCATE(t(iits:iite,kts:kte,jjts:jjte))
  ALLOCATE(qci(iits:iite,kts:kte,2,jjts:jjte))
  ALLOCATE(qrs(iits:iite,kts:kte,3,jjts:jjte))
  ALLOCATE(q(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(den(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(p(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(delz(iims:iime,kms:kme,jjms:jjme))
  ALLOCATE(rain(iims:iime,jjms:jjme))
  ALLOCATE(rainncv(iims:iime,jjms:jjme))
  ALLOCATE(sr(iims:iime,jjms:jjme))
  ALLOCATE(snow(iims:iime,jjms:jjme))
  ALLOCATE(snowncv(iims:iime,jjms:jjme))
  ALLOCATE(graupel(iims:iime,jjms:jjme))
  ALLOCATE(graupelncv(iims:iime,jjms:jjme))

!$OMP PARALLEL DO &
!$OMP PRIVATE ( j ) &
!$OMP SCHEDULE(runtime)
  do j = jjts,jjte
    CALL firstTouch(t(iits,kts,j), q(iims,kms,j)              &
               ,qci(iits,kts,1,j), qrs(iits,kts,1,j)          &
               ,den(iims,kms,j)                               &
               ,p(iims,kms,j), delz(iims,kms,j)               &
               ,j                                             &
               ,rain(iims,j),rainncv(iims,j)                  &
               ,sr(iims,j)                                    &
               ,iids,iide, jjds,jjde, kds,kde                 &
               ,iims,iime, jjms,jjme, kms,kme                 &
               ,iits,iite, jjts,jjte, kts,kte                 &
               ,snow,snowncv                                  &
               ,graupel,graupelncv                            &
                                                              )
  enddo  ! j loop
!$OMP END PARALLEL DO

  ! read remaining input data
  call readarray3(t,'t',unitno,its,ite)
  call readarray4(qci,'qci',unitno,its,ite)
  call readarray4(qrs,'qrs',unitno,its,ite)
  call readarray3(q,'q',unitno,its,ite)
  call readarray3(den,'den',unitno,its,ite)
  call readarray3(p,'p',unitno,its,ite)
  call readarray3(delz,'delz',unitno,its,ite)
  read(unitno) delt
  print *,' delt = ',delt
  read(unitno) g
  print *,' g = ',g
  read(unitno) cpd
  print *,' cpd = ',cpd
  read(unitno) cpv
  print *,' cpv = ',cpv
  read(unitno) t0c
  print *,' t0c = ',t0c
  read(unitno) den0
  print *,' den0 = ',den0
  read(unitno) rd
  print *,' rd = ',rd
  read(unitno) rv
  print *,' rv = ',rv
  read(unitno) ep1
  print *,' ep1 = ',ep1
  read(unitno) ep2
  print *,' ep2 = ',ep2
  read(unitno) qmin
  print *,' qmin = ',qmin
  read(unitno) XLS
  print *,' XLS = ',XLS
  read(unitno) XLV0
  print *,' XLV0 = ',XLV0
  read(unitno) XLF0
  print *,' XLF0 = ',XLF0
  read(unitno) cliq
  print *,' cliq = ',cliq
  read(unitno) cice
  print *,' cice = ',cice
  read(unitno) psat
  print *,' psat = ',psat
  read(unitno) denr
  print *,' denr = ',denr
  call readarray2(rain,'rain',unitno,its,ite)
  call readarray2(rainncv,'rainncv',unitno,its,ite)
  call readarray2(sr,'sr',unitno,its,ite)
  call readarray2(snow,'snow',unitno,its,ite)
  call readarray2(snowncv,'snowncv',unitno,its,ite)
  call readarray2(graupel,'graupel',unitno,its,ite)
  call readarray2(graupelncv,'graupelncv',unitno,its,ite)
  close(unitno)

  ! minimize timer overhead inside OpenMP loop
  ret = gptlinit_handle ('WSM62D', handle)
  ! call WSM6
  ret = gptlstart('WSM62D+OpenMP')
!$OMP PARALLEL DO &
!$OMP PRIVATE ( j,ret ) &
!$OMP SCHEDULE(runtime)
  do j = jjts,jjte
    ret = gptlstart_handle ('WSM62D', handle)
    CALL wsm62D(t(iits,kts,j), q(iims,kms,j)                  &
               ,qci(iits,kts,1,j), qrs(iits,kts,1,j)          &
               ,den(iims,kms,j)                               &
               ,p(iims,kms,j), delz(iims,kms,j)               &
               ,delt,g, cpd, cpv, rd, rv, t0c                 &
               ,ep1, ep2, qmin                                &
               ,XLS, XLV0, XLF0, den0, denr                   &
               ,cliq,cice,psat                                &
               ,j                                             &
               ,rain(iims,j),rainncv(iims,j)                  &
               ,sr(iims,j)                                    &
               ,iids,iide, jjds,jjde, kds,kde                 &
               ,iims,iime, jjms,jjme, kms,kme                 &
               ,iits,iite, jjts,jjte, kts,kte                 &
               ,snow,snowncv                                  &
               ,graupel,graupelncv                            &
                                                              )
    ret = gptlstop_handle ('WSM62D', handle)
  enddo  ! j loop
!$OMP END PARALLEL DO
   ret = gptlstop('WSM62D+OpenMP')

  ! write output data
  PRINT *,'wsm62D():  itimestep == ',itimestep,', write output state'
  fn = 'wsm6_output.dat'
  unitno=31
  open (unitno,file=trim(fn),form="unformatted",action='write', &
        iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR: failed to open output file ',trim(fn), &
               ' . stopping'
    stop
  endif
  call writearray3(t,'t',unitno,its,ite)
  call writearray4(qci,'qci',unitno,its,ite)
  call writearray4(qrs,'qrs',unitno,its,ite)
  call writearray3(q,'q',unitno,its,ite)
  call writearray2(rain,'rain',unitno,its,ite)
  call writearray2(rainncv,'rainncv',unitno,its,ite)
  call writearray2(sr,'sr',unitno,its,ite)
  call writearray2(snow,'snow',unitno,its,ite)
  call writearray2(snowncv,'snowncv',unitno,its,ite)
  call writearray2(graupel,'graupel',unitno,its,ite)
  call writearray2(graupelncv,'graupelncv',unitno,its,ite)
  close(unitno)

  ! deallocate arrays
  DEALLOCATE(t)
  DEALLOCATE(qci)
  DEALLOCATE(qrs)
  DEALLOCATE(q)
  DEALLOCATE(den)
  DEALLOCATE(p)
  DEALLOCATE(delz)
  DEALLOCATE(rain)
  DEALLOCATE(rainncv)
  DEALLOCATE(sr)
  DEALLOCATE(snow)
  DEALLOCATE(snowncv)
  DEALLOCATE(graupel)
  DEALLOCATE(graupelncv)

  ret = gptlstop('Total')
  ret = gptlget_wallclock ('Total', 0, totaltime)  ! The "0" is thread number
  print*,''
  print*,'Total time =' , totaltime

!+---+-----------------------------------------------------------------+

  ! print timing info
  ret = gptlpr (0)
  ret = gptlpr_summary (0)

  END PROGRAM wsm6

