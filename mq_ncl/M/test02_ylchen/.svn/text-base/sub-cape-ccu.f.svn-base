c***********************************************************************
c subroutine name  : cal_cape                                          *
c                                                                      *
c description      : To compute lifting condensation level(LCL),       *
c                    convective condensation level(CCL),level of free  *
c                    convection(LFC), equilibrium level(EL),           *
c                    convective available potential energy(CAPE) and   * 
c                    negative energy(CIN).                             *
c                                                                      *
c I/O parameters   :                                                   *
c  input:  name       type        description                          *
c          p0         real        the initial pressure of air parcel   *
c                                 in unit of hPa.                      *
c          zz0        real        the initial height of air parcel     *
c                                 in unit of meters.                   *
c          t0         real        the initial temperature of air parcel*
c                                 in unit of degree C.                 *
c          td0        real        the initial dew-point temperature of *
c                                 air parcel in unit of degree C.      *
c          nk         integer     the number of sounding profile.      *
c          p(nk)      real array  the pressure of sounding profile in  *
c                                 unit of hPa.                         *
c                                 p(1) > p(2) > p(3) > ... > p(nk)     *
c          z(nk)      real array  the height of sounding profile in    *
c                                 unit of meters.                      *
c                                 z(1) < z(2) < z(3) < ... < z(nk)     *
c          t(nk)      real array  the temperature of sounding profile  *
c                                 unit of degree C.                    *
c                                 t(1) >= t(2) >= t(3) >= ... >= t(nk) *
c  output: name       type        description                          *
c          p_lcl      real        the pressure of LCL in unit of hPa.  *
c          p_ccl      real        the pressure of CCL in unit of hPa.  *
c          p_lfc      real        the pressure of LFC in unit of hPa.  *
c          p_el       real        the pressure of EL in unit of hPa.   * 
c          cape       real        the positive parcel energy of CAPE   *
c                                 in unit of (m**2)(s**2).             *
c          cin        real        the negative parcel energy of CIN    *
c                                 in unit of (m**2)(s**2).             *
c other parameters :                                                   *
c          name       type        description                          *
c          g          real        gravity at sea level.                *
c          cp         real        specific heat of dry air at constant *
c                                 pressure.                            *
c          rgas       real        constant, = 0.622                    *
c                                                                      *
c called fun./sub. : hyd      fit      cal_ta      cal_energy          *
c                                                                      *
c***********************************************************************

       subroutine cal_cape(p0,zz0,t0,td0,nk,p,z,t,p_lcl,p_ccl,p_lfc,    
     +                                                 p_el,cape,cin)

        parameter ( g=9.81,cp=1004.,rgas=287.,aslon=0.622 )
        dimension p(nk),z(nk),t(nk),wxx(901),wyy(901),qws(901)
        

c-----------------------------------------------------------------------
c* To check the input data of air parcel

        if( (p0.lt.100.).or.(p0.gt.1500.) )then
          print*,'pressure of air parcel is error, 100 < p0 < 1500'
          print*,'pressure of air parcel = ',p0,'(hPa)'
          go to 10
        endif
        if( (t0.lt.-150.).or.(t0.gt.100.) )then
          print*,'temperature of air parcel is error, -150 < t0 < 100'
          print*,'temperature of air parcel = ',t0,'(C)'
          go to 10
        endif
        if( (td0.lt.-150.).or.(td0.gt.100.) )then
        print*,'dew-point tem. of air parcel is error, -150 < td0 < 100'
        print*,'dew-point temperature of air parcel = ',td0,'(C)'
         go to 10
        endif
        if( td0.gt.t0 )then
          print*,'input data position of air parcel is error'
          print*,'temperature of air parcel = ',t0,'(C)'
          print*,'dew-point temperature of air parcel = ',td0,'(C)'
         ! td0 = t0
          go to 10
        endif
        if( nk.lt.2 )then
          print*,'the number of input sounding data must be over 2'
          print*,'nk= ',nk
          go to 10
        endif

c-----------------------------------------------------------------------
c* To check temperature of souding data

        do k=1,nk-1,1
         if( (t(k).lt.-150.).or.(t(k).gt.100.) )then
          print*,'temperature of sounding is error, -150 < t < 100'
          do i=1,k
            print*,' i = ',i,' temperature t(i) = ',t(i),'(C)'
          enddo
          go to 10
         endif
         if( (t(k).lt.t(k+1)).and.((t(k+1)-t(k)).gt.10.) )then
          print*,'input the sounding data of temperature is error'
          print*,'the inverse rate of temperature is too large.'
          do i=1,k+1
            print*,' i = ',i,' temperature t(i) = ',t(i),'(C)'
          enddo
          go to 10
         endif
        enddo


c-----------------------------------------------------------------------
c* To check geopotention of souding data

        if (zz0.lt.0.) then
         tbar=0.5*(t0+t(1)+546.3)
         zz0=hyd(z(1),p(1),tbar,p0)
        endif

c-----------------------------------------------------------------------
c..     unit transform from degree C to degree K

        tt=t0+273.15
        ttd=td0+273.15

c-----------------------------------------------------------------------
c* To compute potential temperature of air parcel.

        theta0=tt*(1000./p0)**(rgas/cp)

c-----------------------------------------------------------------------
c* To compute isentropic condensation temperature(K) of air parcel.

        tc0=1./(1./(ttd-56.)+log(tt/ttd)/800.)+56.

c-----------------------------------------------------------------------
c* To compute mixing ratio of air parcel.

        p_w0=6.112*exp(17.67*td0/(td0+243.5))
        w0=aslon*p_w0/(p0-p_w0)
        q0=1000.*w0

c-----------------------------------------------------------------------
c* To compute equivalent potential temperature of air parcel.

        if( p0.gt.100. )then
          thetae0=theta0*exp(2675.*w0/tc0)
        else
          thetae0=theta0
        endif

c-----------------------------------------------------------------------
c* To compute the height of LCL.

        p_lcl=p0*(tc0/tt)**(cp/rgas)
        t_lcl=tc0-273.15

        if ( p_lcl.ge.p(1) ) then
           k_lcl=1
           go to 40
        else
           do k=1,nk-1
           if( (p_lcl.le.p(k)).and.(p_lcl.gt.p(k+1)) )then
               tbar1=0.5*(t(k)+t_lcl+546.3)
               tbar2=0.5*(t(k+1)+t_lcl+546.3)
               za1=hyd(z(k),p(k),tbar1,p_lcl)
               za2=hyd(z(k+1),p(k+1),tbar2,p_lcl)
               h_lcl=0.5*(za1+za2) 
               k_lcl=k+1
               go to 30
           endif
           enddo
        endif

c-----------------------------------------------------------------------
c* To compute the height of CCL.

 30     do k=1,nk
           p_w=6.112*exp(17.67*t(k)/(t(k)+243.5))
           qws(k)=622.*p_w/(p(k)-p_w)
        enddo

        do k=1,nk-1
          if( (q0.le.qws(k)).and.(q0.gt.qws(k+1)) )then
              deta=(p_lcl-p0)*(t(k+1)-t(k))-(p(k+1)-p(k))*(t_lcl-td0)
          if( abs(deta).gt.0.0001 )then
              dn=( (p(k)-p0)*(t_lcl-td0)-(p_lcl-p0)*(t(k)-td0) )/deta
              t_ccl=t(k)+(t(k+1)-t(k))*dn
              p_ccl=p(k)+(p(k+1)-p(k))*dn
              tbar1=0.5*(t(k)+t_ccl+546.3)
              tbar2=0.5*(t(k+1)+t_ccl+546.3)
              za1=hyd(z(k),p(k),tbar1,p_ccl)
              za2=hyd(z(k+1),p(k+1),tbar2,p_ccl)
              h_ccl=0.5*(za1+za2)
              k_ccl=k+1
ccccc         go to 40
              go to 999
          endif
          endif
        enddo

c-----------------------------------------------------------------------
c* To compute the height of LFC.

 40     if( t_lcl.ge.t(k_lcl) )then
            p_lfc=p_lcl
            h_lfc=h_lcl
            t_lfc=t_lcl
            k_lfc=k_lcl
            cin=-99.99
            go to 60
        endif

           t_air1=t_lcl
        do k=k_lcl+1,nk
           t_air0=t_air1
           call cal_ta(thetae0,p(k),t_air0,0.001,50000,t_air1,ier)

           if( t_air1.ge.t(k) )then

               detx=(t_air1-t_air0)-(t(k)-t(k-1))
               if( abs(detx).lt.0.0000001 ) go to 10
               dm=(t(k-1)-t_air0)/detx
               rpa1=log(p(k))+(log(p(k))-log(p(k-1)))*dm
               pa1=exp(rpa1)
               te1=t_air0+(t_air1-t_air0)*dm
               call cal_ta(thetae0,pa1,t_air0,0.001,50000,ta1,ier)
               if( ier.eq.2 ) go to 10
               detx=(ta1-t_air0)-(te1-t(k-1))
               if( abs(detx).lt.0.0000001 ) go to 10
               dm=(t(k-1)-t_air0)/detx
               rpa2=log(p(k-1))+(rpa1-log(p(k-1)))*dm
               pa2=exp(rpa2)
               te2=t_air0+(ta1-t_air0)*dm
               call cal_ta(thetae0,pa2,t_air0,0.001,50000,ta2,ier)
               p_lfc=pa2
               t_lfc=ta2
               tbar1=0.5*(t(k-1)+t_lfc+546.3)
               tbar2=0.5*(t(k)+t_lfc+546.3)
               za1=hyd(z(k-1),p(k-1),tbar1,p_lfc)
               za2=hyd(z(k),p(k),tbar2,p_lfc)
               h_lfc=0.5*(za1+za2)
               k_lfc=k
               go to 60
           endif
        enddo
        go to 10
 60     if( td0.eq.t0 )then
            p_lcl=p0
            h_lcl=zz0
            t_lcl=t0
            p_ccl=p0
            t_ccl=t0
            p_lfc=p0
            h_lfc=zz0
            t_lfc=t0
            go to 65
        endif

c-----------------------------------------------------------------------
c* To compute CIN.

        if( k_lfc.ge.2 )then
            nn=0
            do k=1,k_lfc
               t_a=theta0*(p(k)/1000.)**(rgas/cp)-273.15
               value=g*(t(k)-t_a)/(t(k)+273.15)
               if( value.ge.0.0 )then
                   nn=nn+1
                   wxx(nn)=z(k)
                   wyy(nn)=value
               endif
            enddo
            nn=nn+1
            wxx(nn)=h_lfc
            wyy(nn)=0.
            call cal_energy(nn,wxx,wyy,cin)
        endif
 65     continue

c-----------------------------------------------------------------------
c* To compute the height of EL.

           t_air1=t_lfc

        do k=k_lfc+1,nk
           t_air0=t_air1
           call cal_ta(thetae0,p(k),t_air0,0.001,50000,t_air1,ier)

           if( t(k).gt.t_air1 )then
               detx=(t_air1-t_air0)-(t(k)-t(k-1))
               if( abs(detx).lt.0.0000001 ) go to 10
               dm=(t(k-1)-t_air0)/detx
               rpa1=log(p(k-1))+(log(p(k))-log(p(k-1)))*dm
               pa1=exp(rpa1)
               te1=t_air0+(t_air1-t_air0)*dm
               call cal_ta(thetae0,pa1,t_air0,0.001,50000,ta1,ier)
               if( ier.eq.2 ) go to 10
               detx=(ta1-t_air0)-(te1-t(k-1))
               if( abs(detx).lt.0.0000001 ) go to 10
               dm=(t(k-1)-t_air0)/detx
               rpa2=log(p(k-1))+(rpa1-log(p(k-1)))*dm
               pa2=exp(rpa2)
               te2=t_air0+(ta1-t_air0)*dm
               call cal_ta(thetae0,pa2,t_air0,0.001,50000,ta2,ier)
               if( ier.eq.2 ) go to 10
               p_el=pa2
               t_el=ta2
               tbar1=0.5*(t(k-1)+t_el+546.3)
               tbar2=0.5*(t(k)+t_el+546.3)
               za1=hyd(z(k-1),p(k-1),tbar1,p_el)
               za2=hyd(z(k),p(k),tbar2,p_el)
               h_el=0.5*(za1+za2)
               k_el=k-1
               go to 75
           endif
        enddo

c-----------------------------------------------------------------------

 75     if( k_lfc.gt.k_el ) go to 10

c-----------------------------------------------------------------------
c* To compute convective available potential energy(CAPE)
c*    from h_lfc to h_el.

           nn=1
           wxx(nn)=h_lfc
           wyy(nn)=0.
           t_air1=t_lfc

        do k=k_lfc,k_el,1
           t_air0=t_air1
           call cal_ta(thetae0,p(k),t_air0,0.001,50000,t_air1,ier)
           value=g*(t_air1-t(k))/(t(k)+273.15)
           if( value.ge.0.0 )then
               nn=nn+1
               wxx(nn)=z(k)
               wyy(nn)=value
           endif
        enddo

           nn=nn+1
           wxx(nn)=h_el
           wyy(nn)=0.

        call cal_energy(nn,wxx,wyy,cape)

        return

 10     cape=-99.99
        p_lcl=-99.99
        p_ccl=-99.99
        p_lfc=-99.99
        p_el=-99.99
        cin=-99.99
 999    return

        end

c***********************************************************************
c subroutine name  : cal_ta                                            *
c                                                                      *
c description      : To find temperature(ta) at constant thetae and    *
c                    pressure(pa).                                     * 
c                                                                      *
c I/O parameters   :                                                   *
c  input:  name       type        description                          *
c          thetae     real        theta E in unit of degree K.         *
c          pa         real        pressure in unit of hPa.             *
c          ti         real        the first guess temperature in unit  *
c                                 of degree C, the first temperature   *
c                                 must be greater than ta.             *
c          det        real        the interval temperature.            *
c          nx         integer     the number of iteration.             *
c  output: name       type        description                          *
c          ta         real        the tempertaure in unit of degree C. *
c          ier        integer     the error message,                   *
c                                 = 0, no error.                       *
c                                 = 2, cannot to get ta.               *
c                                                                      *
c called fun./sub. : none                                              *
c                                                                      *
c***********************************************************************
        subroutine cal_ta(thetae,pa,ti,det,nx,ta,ier)

        ier=0
        if( pa.lt.100. )then
          ta=thetae*(pa/1000.)**(287./1004.) - 273.15
          return
        endif
        ta=ti
        do i=1,nx
          ta=ta-det
          theta=(ta+273.15)*(1000./pa)**(287./1004.)
          p_wa=6.112*exp(17.67*ta/(ta+243.5))
          wa=0.622*p_wa/(pa-p_wa)
          thetae_a=theta*exp(2675.*wa/(ta+273.15))
          if( thetae_a.le.thetae )return
        enddo
        ier=2
        return
        end

c***********************************************************************
c subroutine name  : cal_energy                                        *
c                                                                      *
c description      : To compute energy or integral.                    *
c                                                                      *
c I/O parameters   :                                                   *
c  input:  name       type        description                          *
c          nn         integer     the number of input data.            *
c          x(nn)      real array  the x-coordinate.                    *
c                                 = g*(T(parcel)-T(envir.))/T(envir.)  *
c          y(nn)      real array  the y-coordinate.                    *
c                                 = height in unit of meters           *
c  output: name       type        description                          *
c          energy     real        the energy or integral.              *
c                                                                      *
c called fun./sub. : rsum                                              *
c                                                                      *
c***********************************************************************
        subroutine cal_energy(nn,x,y,energy)
        dimension x(nn),y(nn)
        if( nn.le.2 )then
            energy=-99.99
            return
        endif
        energy=0.
        do 10 k=1,nn-1
          if( k.eq.1 )then
            a1=rsum(x(k),x(k+1),x(k+2),y(k),y(k+1),y(k+2),x(k),x(k+1))
            energy=energy+a1
            go to 10
          endif
          if( k.eq.(nn-1) )then
            a2=rsum(x(k-1),x(k),x(k+1),y(k-1),y(k),y(k+1),x(k),x(k+1))
            energy=energy+a2
            go to 10
          endif
          a1=rsum(x(k),x(k+1),x(k+2),y(k),y(k+1),y(k+2),x(k),x(k+1))
          a2=rsum(x(k-1),x(k),x(k+1),y(k-1),y(k),y(k+1),x(k),x(k+1))
          energy=energy+0.5*(a1+a2)
 10     continue
        return
        end

c***********************************************************************
c function name    : rsum                                              *
c                                                                      *
c description      : To compute integral from point P to Q by          *
c                    parabolic curved fitting.                         *
c                                                                      *
c I/O parameters   :                                                   *
c  input:  name       type        description                          *
c          x1         real        the x-coordinate at point 1.         *
c          x2         real        the x-coordinate at point 2.         *
c          x3         real        the x-coordinate at point 3.         *
c          y1         real        the y-coordinate at point 1.         *
c          y2         real        the y-coordinate at point 2.         *
c          y3         real        the y-coordinate at point 3.         *
c          p          real        the x-coordinate at point P.         *
c          q          real        the x-coordinate at point Q.         *
c  output: name       type        description                          *
c          rsum       real        the integral from P to Q.            *
c                                                                      *
c called fun./sub. : none                                              *
c                                                                      *
c***********************************************************************
        function rsum(x1,x2,x3,y1,y2,y3,p,q)
        x12=x1-x2
        x13=x1-x3
        x23=x2-x3
        if( abs(x12).le.1.0E-10 )x12=1.0E-10
        if( abs(x13).le.1.0E-10 )x13=1.0E-10
        if( abs(x23).le.1.0E-10 )x23=1.0E-10
        a=y1/(x12*x13)
        b=-y2/(x12*x23)
        c=y3/(x13*x23)
        rsum=(a+b+c)*(q*q*q-p*p*p)/3.-
     +       0.5*( (x2+x3)*a+(x1+x3)*b+(x1+x2)*c )*(q*q-p*p)+
     +       ( x2*x3*a+x1*x3*b+x1*x2*c )*(q-p)
        return
        end

c***********************************************************************
c function name    : hyd                                               *
c                                                                      *
c description      : To interpolate height by hydrostatic              *
c                    approximation.                                    *
c                                                                      *
c I/O parameters   :                                                   *
c  input:  name       type        description                          *
c          z          real        the known height in unit of meters.  *
c          p          real        the known pressure in unit of hPa.   *
c          t          real        the known temperature in unit of     *
c                                 degree C.                            *
c          p0         real        the interpolated pressure in init of *
c                                 hPa.                                 *
c  output: name       type        description                          *
c          hyd        real        the height in unit of meters.        *
c                                                                      *
c called fun./sub. : none                                              *
c                                                                      *
c***********************************************************************
        function hyd(z,p,t,p0)
        hyd=z+287.*t*log(p/p0)/9.81
        return
        end 

c***********************************************************************
c function name    : fit                                               *
c                                                                      *
c description      : To interpolate y-value at point x by linear method*
c                                                                      *
c I/O parameters   :                                                   *
c  input:  name       type        description                          *
c          x1         real        the x-coordinate at Point 1.         *
c          y1         real        the y-coordinate at Point 1.         *
c          x2         real        the x-coordinate at Point 2.         *
c          y2         real        the y-coordinate at Point 2.         *
c          x          real        the x-coordinate at Point x.         *
c  output: name       type        description                          *
c          fit        real        the y-value at point x.              *
c                                                                      *
c called fun./sub. : none                                              *
c                                                                      *
c***********************************************************************
        function fit(x1,y1,x2,y2,x)
        fit=y1+(x-x1)*(y2-y1)/(x2-x1)
        return
        end
