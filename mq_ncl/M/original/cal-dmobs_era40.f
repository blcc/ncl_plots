c..*********************************************************************
c..    Calculate Deltam from observational data and some important
c..    thermodynamic parameters (with new definition of deltam)
c..         (after separation of variables)
c..
c..    Jia-yuh Yu (1993/11)
c..    rewrited by Yin-Min Cho (2005/10/26)
c..*********************************************************************
c..   nx     : points in x-direction
c..   ny     : points in y-direction
c..   nm     : months
c..   nz     : levels
c..   sfc    : the based of the vertical direction [mb]
c..   pbl    : the planetary boundary layer [mb]
c..   toa    : the top of the vertical direction [mb]
c..   mbt    : interval of vertical direction [mb]
c..   z      : geopotental hight [m]
c..   t      : temperature [K]
c..   q      : specific humidity [kg/kg]
c..   p      : pressure [Pascale]
c..   iop    : 0 for fixed cloud-top level at 150mb
c..          : 1 for determined by cloud-top level
c..          : 2 for determined by lifting condensation level (LCL)
c..          : 3 for determined by convective condensation level (CCL)
c..   deltam : gross moist stability [J*Kg**-1]
c..   deltamq: gross moisture stratification [J*Kg**-1]
c..*********************************************************************
      parameter ( nx=144      , ny= 73      , nm=12   , nz=12
     &,           sfc=1000.   , toa=100.    , mbt=1
     &,           nz1=toa/mbt , nz2=sfc/mbt , dp=100.*mbt
     &,           iop=1       , ify=1958    , spv=-99.99 )

      real l , kappa

      parameter ( l=2.44e+6 , cp=1.005e+3 , rv=4.6151e+2 , g=9.81
     &,           kappa=0.286 )

      character fn(5)*100 , year*4   , vname*10 , zname*10
     &,         path1*29  , path2*30 , cas*3

      dimension p(nz2) , dd(nx,ny)     , lev(nz)
     &,         z(nz2) , x1(nx,ny,nz)  , zz1(nz2-nz1+1)
     &,         t(nz2) , x2(nx,ny,nz)  , tt1(nz2-nz1+1)
     &,         q(nz2) , x3(nx,ny,nz)  , pp1(nz2-nz1+1)
     &,         h(nz2) , a(nz2)        , gamma(nz2)
     &,         ps(nx,ny)  , ts(nx,ny)  , ds(nx,ny)
     &,         ah(nx,ny)  , aph(nx,ny) , aplus(nz2)
     &,         top(nx,ny) , deltam(nx,ny)
     &,         bot(nx,ny) , deltamq(nx,ny)
     &,         d1(nx,ny,nz2) , d2(nx,ny,nz2) , d3(nx,ny,nz2)

      data lev /  100,  150,  200,  250,  300,  400,  500,  600,
     &            700,  850,  925, 1000/

c..*********************************************************************
c.. open I/O file and read data
c..*********************************************************************
      print*,'********************************************'
      print*,'please input start year  (ex: 1958)'
      print*,'********************************************'
      read(*,'(i4)') isyear
      print*,'********************************************'
      print*,'please input  end  year  (ex: 2001)'
      print*,'********************************************'
      read(*,'(i4)') ieyear

      path1 = '/data1/ylchen/DATA_old/ERA40/'
      path2 = '/data1/ylchen/DATA_out/DM_CLT/'
      if (iop.eq.0) cas = 'FIX'
      if (iop.eq.1) cas = 'CLT'
      if (iop.eq.2) cas = 'LCL'
      if (iop.eq.3) cas = 'CCL'

      do 1000 ity = isyear , ieyear
      
      write(year,'(i4)') ity

      open (21,file=path2//'ERA40_DM_'//year//'_mon_'//cas//'.bin'
     &,status='unknown',form='unformatted',access='direct',recl=4*nx*ny)
c      open (22,file=path2//'ERA40_DMQ_'//year//'_pnt_'//cas//'.bin'
c     &,status='unknown',form='unformatted',access='direct',recl=4*nx*ny)
c      open (23,file=path2//'ERA40_CLT_'//year//'_pnt_'//cas//'.bin'
c     &,status='unknown',form='unformatted',access='direct',recl=4*nx*ny)
c     open (24,file=path2//'ERA40_GAM_'//year//'_mon_'//cas//'.bin'
c    &,status='unknown',form='unformatted',access='direct',recl=nx*ny)
c     open (25,file=path2//'ERA40_A_'//year//'_mon_'//cas//'.bin'
c    &,status='unknown',form='unformatted',access='direct',recl=nx*ny)
c     open (26,file=path2//'ERA40_AH_'//year//'_mon_'//cas//'.bin'
c    &,status='unknown',form='unformatted',access='direct',recl=nx*ny)
c     open (27,file=path2//'ERA40_AP_'//year//'_mon_'//cas//'.bin'
c    &,status='unknown',form='unformatted',access='direct',recl=nx*ny)
c     open (28,file=path2//'ERA40_APH_'//year//'_mon_'//cas//'.bin'
c    &,status='unknown',form='unformatted',access='direct',recl=nx*ny)
      if (iop.ge.2)
     &   open (29,file=path2//'ERA40_CLB_'//year//'_mon_'//cas//'.bin'
     &,status='unknown',form='unformatted',access='direct',recl=4*nx*ny)

c..*********************************************************************

      do 900 itm = 1 , nm

         icnt = itm + (ity-ify)*nm

      do 150 ivarb = 1 , 3

         call find_var(ivarb,vname,numv)

      do 150 k = 1 , nz

         call find_lev(lev(k),zname,numz)

         fn(1) = path1//'prs/'//vname(1:numv)//'/'//zname(1:numz)//
     &           'hPa/ERA40_'//vname(1:numv)//'_'//zname(1:numz)//
     &           'hPa_1958_2001_pnt.bin'
c        print*,'input file=',fn(1)

         open (11,file=fn(1),status='old',form='unformatted'
     &                      ,access='direct',recl=4*nx*ny)

         read (11,rec=icnt) dd

      do 110 j = 1 , ny
      do 110 i = 1 , nx
         if (ivarb.eq.1) x1(i,j,k)=dd(i,j)
         if (ivarb.eq.2) x2(i,j,k)=dd(i,j)+273.15
         if (ivarb.eq.3) x3(i,j,k)=dd(i,j)*10.0**-3
 110  continue

 150  continue

c***********************************************************************

      if (iop.ge.2) then
         fn(1) = path1//'sfc/2T/ERA40_2T_SFC_1958_2001_pnt.bin'
         fn(2) = path1//'sfc/2D/ERA40_2D_SFC_1958_2001_pnt.bin'
         fn(3) = path1//'sfc/MSL/ERA40_MSL_SFC_1958_2001_pnt.bin'
         open (12,file=fn(1),status='old',form='unformatted'
     &                      ,access='direct',recl=4*nx*ny)
         open (13,file=fn(2),status='old',form='unformatted'
     &                      ,access='direct',recl=4*nx*ny)
         open (14,file=fn(3),status='old',form='unformatted'
     &                      ,access='direct',recl=4*nx*ny)

         read (12,rec=icnt) ts
         read (13,rec=icnt) ds
         read (14,rec=icnt) ps
      endif

c..*********************************************************************
c.. define vertical levels of integration ***
c..*********************************************************************
      print*,'calculate DM & DMQ -  year = ',ity,' month = ',itm

      do 210 k = 1 , nz2
         p(k)=float(k)*100.*float(mbt)
 210  continue

c..*********************************************************************
c.. vertically linear interpolation ***
c..*********************************************************************

      do 800 j = 1 , ny
      do 800 i = 1 , nx

      do 250 level = 1 , nz-1

         if (level.eq.1) then
            k1=100/mbt
            k2=150/mbt
         else if (level.eq.2) then
            k1=150/mbt
            k2=200/mbt
         else if (level.eq.3) then
            k1=200/mbt
            k2=250/mbt
         else if (level.eq.4) then
            k1=250/mbt
            k2=300/mbt
         else if (level.eq.5) then
            k1=300/mbt
            k2=400/mbt
         else if (level.eq.6) then
            k1=400/mbt
            k2=500/mbt
         else if (level.eq.7) then
            k1=500/mbt
            k2=600/mbt
         else if (level.eq.8) then
            k1=600/mbt
            k2=700/mbt
         else if (level.eq.9) then
            k1=700/mbt
            k2=850/mbt
         else if (level.eq.10) then
            k1=850/mbt
            k2=925/mbt
         else if (level.eq.11) then
            k1=925/mbt
            k2=1000/mbt
         endif

      do 220 k = k1 , k2
         z(k)=x1(i,j,level)+(k-k1)*(x1(i,j,level+1)
     &        -x1(i,j,level))/float(k2-k1)
         t(k)=x2(i,j,level)+(k-k1)*(x2(i,j,level+1)
     &        -x2(i,j,level))/float(k2-k1)

      if (k.lt.300/mbt) then
         q(k)=0.
      else
         q(k)=x3(i,j,level)+(k-k1)*(x3(i,j,level+1)
     &        -x3(i,j,level))/float(k2-k1)
c        print*,z(k)
c        print*,t(k)
c        print*,q(k)
      endif
 220  continue

 250  continue

c..*********************************************************************
c  *** calculate gamma profile ***
c..*********************************************************************

      do 300 k = nz1 , nz2
         e0=100.*6.112*exp(17.67*(t(nz2)-273.15)/(t(nz2)-273.15+243.5))
         gamma(k)=0.622*(l**2)*e0/(rv*cp)/(t(k)**2)/p(k)*
     &            exp(l/rv*(1./t(nz2)-1./t(k)))
 300  continue

c..*********************************************************************
c  *** calculate moist-static energy ***
c..*********************************************************************

      do 320 k = nz1 , nz2
         h(k)=cp*t(k)+l*q(k)+g*z(k)
 320  continue

c..*********************************************************************
c  *** estimate (or specify) cloud-top level ***
c..*********************************************************************
c***********************************************************************
c  << data: extend from 1000mb up to 100mb >>
c
c  << using fixed cloud-top level: 150mb >>
c     ktop : cloud-top index
c***********************************************************************

      if (iop.eq.0) ktop=150/mbt

c***********************************************************************
c  << estimate cloud-top level using h_sat=h_1000mb >>
c  Note 1: check cloud-top level from PBL(950mb) up to 100mb
c  Note 2: lowest ct level set at 450mb
c  Note 3: highest ct level determined by basic states
c..*********************************************************************

      if (iop.eq.1) then

      pbl=950.

      do 340 k = int(pbl)/mbt , nz1 , -1

         dif=h(nz2)-h(k)

      if (h(nz2).le.h(nz1)) then
         if (dif.le.0.) then
            ktop=k+1
            go to 350
         endif
         ktop=450/mbt
      else
         ktop=nz1
      endif

 340  continue

 350  if (ktop.ge.450/mbt) ktop=450/mbt

      endif

c***********************************************************************
c  << estimate cloud-based level using LCL or CCL >>
c***********************************************************************

      if (iop.ge.2) then

      pbl=0.

      if ((ts(i,j).eq.spv).or.(ds(i,j).eq.spv).or.(ps(i,j).eq.spv))
     &   go to 700

         do n = 1000 , 100 , -1
            zz1(nz2-n+1)=z(n)
            tt1(nz2-n+1)=t(n)-273.15
            pp1(nz2-n+1)=p(n)*0.01
         enddo

         zz0=-99.99
         pp0=ps(i,j)
         tt0=ts(i,j)
         ds0=ds(i,j)

         call cal_cape(pp0,zz0,tt0,ds0,nz2-nz1+1,pp1,zz1,tt1
     &                ,p_lcl,p_ccl,p_lfc,p_el,cape,cin)
c        print*,p_ccl

         if (iop.eq.2) pbl=p_lcl
         if (iop.eq.3) pbl=p_ccl
         if (pbl.gt.1000.) pbl=1000.
         if ((pbl.eq.0.).or.(pbl.eq.spv)) go to 700
         kbot=int(pbl)/mbt

         do 360 k = int(pbl)/mbt-1 , nz1 , -1

            dif=h(kbot)-h(k)

         if (h(kbot).le.h(nz1)) then
            if (dif.le.0.) then
               ktop=k+1
               go to 370
            endif
            ktop=450/mbt
         else
            ktop=nz1
         endif

 360  continue

 370  if (ktop.ge.450/mbt) ktop=450/mbt

      endif

c..*********************************************************************
c  *** calculate hhat and qhat ***
c..*********************************************************************

      hhat=0.
      qhat=0.
      kk=ktop

      do 380 k = kk , nz2
         hhat=hhat+h(k)/float(nz2-kk+1)
         qhat=qhat+l*q(k)/float(nz2-kk+1)
 380  continue

c..*********************************************************************
c  *** calculate A(p,x,y) profile ***
c..*********************************************************************
c***********************************************************************
c  << Above PBL (cloud base) >>
c***********************************************************************

      npbl=int(sfc-pbl)/mbt

      do 410 k = nz1 , nz2-npbl

      sum=0.

      do 400 kk = k , nz2-npbl
         sum=sum+(1./(1.+gamma(kk)))/p(kk)*dp
 400  continue

      a(k)=(1./(1.+gamma(k)))*exp(-(kappa*sum))

 410  continue

c***********************************************************************
c  << Below PBL (cloud base) >>
c***********************************************************************

      do 420 k = nz2-npbl , nz2
         a(k)=(1./(1.+gamma(nz2-npbl)))*(float(k*mbt)/pbl)**kappa
 420  continue

c..*********************************************************************
c  *** calculate ahat ***
c..*********************************************************************

      ahat=0.
      kk=ktop

      do 440 k = kk , nz2
         ahat=ahat+a(k)/float(nz2-kk+1)
 440  continue

c..*********************************************************************
c  *** calculate aplus(p,x,y) profile ***
c..*********************************************************************

      do 470 k = nz1 , nz2-1

      sum=0.

      do 460 kk = k , nz2-1
         sum=sum+a(kk)/p(kk)*dp
 460  continue

      aplus(k)=sum

 470  continue

c..*********************************************************************
c  *** calculate aplushat ***
c..*********************************************************************

      aplushat=0.
      kk=ktop

      do 480 k = kk , nz2
         aplushat=aplushat+aplus(k)/float(nz2-kk+1)
 480  continue

c..*********************************************************************
c  *** calculate hahat & qahat ***
c..*********************************************************************

      hahat=0.
      qahat=0.
      kk=ktop

      do 500 k = kk , nz2
         hahat=hahat+(h(k)*aplus(k))/float(nz2-kk+1)
         qahat=qahat+(l*q(k)*aplus(k))/float(nz2-kk+1)
 500  continue

c..*********************************************************************
c  *** calculate deltam with fixed/flexible cloud-top level ***
c..*********************************************************************
c***********************************************************************
c  Note: No value for deltam once the cloud-top level exists
c        below 450 mb
c***********************************************************************

      deltam(i,j)=hahat-aplushat*hhat
      deltamq(i,j)=-(qahat-aplushat*qhat)
      ah(i,j)=ahat
      aph(i,j)=aplushat

      if ((ktop.lt.(100/mbt)).or.(ktop.ge.(450/mbt))) then

         deltam(i,j)=spv
         deltamq(i,j)=spv
         top(i,j)=spv
         if (iop.ge.2) bot(i,j)=spv

      else

         top(i,j)=float(ktop*mbt)
         if (iop.ge.2) bot(i,j)=float(kbot*mbt)

      endif

      do 600 k = 1 , nz2
         d1(i,j,k)=a(k)
         d2(i,j,k)=gamma(k)
         d3(i,j,k)=aplus(k)
 600  continue

      go to 800

 700  deltam(i,j)=spv
      deltamq(i,j)=spv
      top(i,j)=spv
      bot(i,j)=spv

 800  continue

      print *,'Delta M(73,37)',deltam(73,37)
c      print *,'Delta Mq(73,37)',deltamq(73,37)
      print *,'k-top=',top(73,37)
      print *,'k-bot=',bot(73,37)

      write (21,rec=itm) ((deltam(i,j),i=1,nx),j=1,ny)
      write (22,rec=itm) ((deltamq(i,j),i=1,nx),j=1,ny)
      write (23,rec=itm) ((top(i,j),i=1,nx),j=1,ny)
c     write (26,rec=itm) ((ah(i,j),i=1,nx),j=1,ny)
c     write (28,rec=itm) ((aph(i,j),i=1,nx),j=1,ny)
      if (iop.ge.2) write (29,rec=itm) ((bot(i,j),i=1,nx),j=1,ny)

      do 880 k = 1 , nz2

      icmt = k + (itm-1)*nz2

c     write (24,rec=icmt) ((d1(i,j,k),i=1,nx),j=1,ny)
c     write (25,rec=icmt) ((d2(i,j,k),i=1,nx),j=1,ny)
c     write (27,rec=icmt) ((d3(i,j,k),i=1,nx),j=1,ny)

 880  continue

 900  continue
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)
      if (iop.ge.2) close(29)

1000  continue
      stop
      end

c***********************************************************************
      subroutine find_var(ivar,name,num)
c***********************************************************************
      character*10 name
      dimension ivnum(3)
      data ivnum/ 1 , 1 , 1 /

      if (ivar.eq.1) then
          name = 'Z'
      elseif (ivar.eq.2) then
          name = 'T'
      elseif (ivar.eq.3) then
          name = 'Q'
      endif
      num=ivnum(ivar)

      return
      end

c***********************************************************************
      subroutine find_lev(ilev,name,num)
c***********************************************************************
      character*10 name

      if (ilev.eq.1000) then
         write (name,'(i4)') ilev
         num = 4
      elseif ((ilev.lt.1000).and.(ilev.ge.100)) then
         write (name,'(i3)') ilev
         num = 3
      elseif (ilev.lt.100) then
         write (name,'(i2)') ilev
         num = 2
      endif

      return
     end
