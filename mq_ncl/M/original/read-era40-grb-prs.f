c***********************************************************************
c     This program use to read ERA40 daily & monthly datasets
c
c     Yin-Min, Cho (2005/07/07)
c***********************************************************************
c***********************************************************************
c.. Pressure level variable index
c***********************************************************************
c..  10  Z    m**2/s**2    Geopotential (orography at surface)
c..  11  T    K            Temperature
c..  12  U    m/s          U-Velocity
c..  13  V    m/s          V-Velocity
c..  14  Q    kg/kg        Specific Humidity
c..  15  W    Pa/s         Vertical Velocity
c..  20  PV   K*m**2/kg*s  Potential Voticity
c..  21  RH   %            Relative Humidity
c..  22  O3   kg/kg        Ozone Mass Mixing Ratio
c..  30  VOR  1/s          Vorticity (relative)
c..  31  DIV  1/s          Divergence
c***********************************************************************
c.. level index
c***********************************************************************
c..   1000hPa , 925hPa , 850hPa , 775hPa , 700hPa , 600hPa , 500hPa
c..    400hPa , 300hPa , 250hPa , 200hPa , 150hPa , 100hPa ,  70hPa
c..     50hPa ,  30hPa ,  20hPa ,  10hPa ,   7hPa ,   5hPa ,   3hPa
c..      2hPa ,   1hPa
c***********************************************************************
c***********************************************************************

      parameter ( nx=144  , ny=73  , spv=-99.99 )

      dimension dd(nx,ny) , ddm(nx,ny)

      character fn(5)*100 ,  year*4 ,  var*3  ,  vlev*4
     &,         csi*3     ,  cas*3  ,  sey*9  ,  dir*100
     &,         mmm*3     ,  stp1*4 ,  stp2*4 ,  cnt*5

c***********************************************************************

      call read_dir(dir,ndir)

      call find_mac(nx,ny,nxy)

      call find_var(ivar,var,nnv,trans)

      call find_level(lev,vlev,nlev)

      call find_date(isy,iey,inp,nh,ist,iet,csi,cas,sey)

c***********************************************************************
      fn(2) = dir(1:ndir)//var(1:nnv)//'/'//vlev(1:nlev)//'hPa/'//
     &        'ERA40_'//var(1:nnv)//'_'//vlev(1:nlev)//
     &        'hPa_'//sey//'_'//cas//'.bin'

      if (isy.eq.iey) fn(2) = dir(1:ndir)//var(1:nnv)//'/'//
     &         vlev(1:nlev)//'hPa/ERA40_'//var(1:nnv)//'_'//
     &         vlev(1:nlev)//'hPa_'//sey(1:4)//'_'//cas//'.bin'

      open (21,file=fn(2),status='unknown',form='unformatted'
     &                   ,access='direct',recl=nxy)
c***********************************************************************

      icnt = 0

      do 999 ity = isy , iey

         write (year,'(i4)') ity

         fn(1) = dir(1:ndir)//'ERA40_'//var(1:nnv)//'_'//
     &           vlev(1:nlev)//'hPa_'//year//'_'//csi//'.bin'

         open (11,file=fn(1),status='old',form='unformatted'
     &                      ,access='direct',recl=nxy)

         print*,'input filename = ', fn(1)

c**********************************************************************

      if ((nh.eq.1).and.(ity.ne.1957).and.(ity.ne.2002).and.
     &    ((iet.eq.365).or.(iet.eq.366))) then
         if (mod(ity,4).ne.0) iet = 365
         if (mod(ity,4).eq.0) iet = 366
      endif

      do 888 itm = ist , iet

         call time_step(ity,itm,inp,nh,nn,nnm)

         do j = 1 , ny
         do i = 1 , nx
            ddm(i,j) = 0.
         enddo
         enddo

      do 444 n = nnm+1 , nnm+nn

         read(11,rec=n,err=888) ((dd(i,ny-j+1),i=1,nx),j=1,ny)

      do 333 j = 1 , ny
      do 333 i = 1 , nx
         if ((ivar.eq.11).and.(dd(i,j).ne.spv))
     &       dd(i,j) = dd(i,j) - 273.15
         ddm(i,j) = ddm(i,j) + dd(i,j) / (float(nn)*trans)
 333  continue

 444  continue

      icnt = icnt + 1
      if ((nh.eq.0).and.(ity.le.1957).and.(itm.lt.9)) icnt = icnt - 1
      if ((nh.eq.1).and.(ity.le.1957).and.(itm.lt.244)) icnt = icnt - 1
      if ((nh.eq.5).and.(ity.le.1957).and.(itm.lt.49)) icnt = icnt - 1
      write(21,rec=icnt,err=888) ((ddm(i,j),i=1,nx),j=1,ny)

      write(mmm ,'(i3)') itm
      write(stp1,'(i4)') nnm+1
      write(stp2,'(i4)') nnm+nn
      write(cnt ,'(i5)') icnt
      print*,'Year = ',year,' , Time = ',mmm,'  , Step 1 = ',stp1
     &      ,'  , Step 2 = ',stp2,'  , CNT = ',cnt

 888  continue

 999  continue

      close(21)
      print*,fn(2)

      stop
      end

c**********************************************************************
      subroutine read_dir(dir,ndir)
c**********************************************************************
        character dir*100

        print*,'************************************************'
        print*,' Output data directory '
        print*,'************************************************'
        read(*,'(a100)') dir
        read(*,'(i3)') ndir
        print*,'************************************************'

        return
        end

c**********************************************************************
      subroutine find_mac(nx,ny,nxy)
c**********************************************************************

        print*,'************************************************'
        print*,' DEC Alpha ( 0 )   or   LINUX ( 1 ) '
        print*,'************************************************'
        read(*,'(i1)') mac
        print*,'************************************************'

        if (mac.eq.0) nxy = nx*ny*1
        if (mac.ne.0) nxy = nx*ny*4

        return
        end

c***********************************************************************
      subroutine find_var(ivar,var,nnv,trans)
c***********************************************************************
        character var*3

        print*,'************************************************'
        print*,' Please Input Variables '
        print*,'************************************************'
        print*,'  10  => Geopotential [m**2/s**2]              '
        print*,'  11  => Temperature [K]                       '
        print*,'  12  => U-Velocity [m/s]                      '
        print*,'  13  => V-Velocity [m/s]                      '
        print*,'  14  => Specific Humidity [kg/kg]             '
        print*,'  15  => Vertical Velocity [Pa/s]              '
        print*,'  20  => Potential Voticity [K*m**2/kg*s]      '
        print*,'  21  => Relative Humidity [%]                 '
        print*,'  22  => Ozone Mass Mixing Ratio [kg/kg]       '
        print*,'  30  => Vorticity (relative) [1/s]            '
        print*,'  31  => Divergence [1/s]                      '
        print*,'************************************************'
        read(*,'(i2)') ivar
        print*,'************************************************'

        trans = 1.0
        if (ivar.eq.10) trans = 9.8
        if (ivar.eq.14) trans = 10.0**-3
        if (ivar.eq.20) trans = 10.0**-6
        if (ivar.eq.22) trans = 10.0**-6
        if (ivar.eq.30) trans = 10.0**-6
        if (ivar.eq.31) trans = 10.0**-6

        if (ivar.eq.10) var = 'Z'
        if (ivar.eq.11) var = 'T'
        if (ivar.eq.12) var = 'U'
        if (ivar.eq.13) var = 'V'
        if (ivar.eq.14) var = 'Q'
        if (ivar.eq.15) var = 'W'
        if (ivar.eq.20) var = 'PV'
        if (ivar.eq.21) var = 'RH'
        if (ivar.eq.22) var = 'O3'
        if (ivar.eq.30) var = 'VOR'
        if (ivar.eq.31) var = 'DIV'

        if (ivar.le.31) nnv = 3
        if (ivar.le.22) nnv = 2
        if (ivar.le.15) nnv = 1

        return
        end

c***********************************************************************
      subroutine find_level(lev,vlev,nlev)
c***********************************************************************
        character vlev*4

        print*,'******************************************************'
        print*,' Please Input Pressure Levels '
        print*,'******************************************************'
        print*,' 1000 => 1000hPa ,  925 => 925hPa , 850 => 850hPa '
        print*,'                                                  '
        print*,'  775 =>  775hPa ,  700 => 700hPa , 600 => 600hPa '
        print*,'                                                  '
        print*,'  500 =>  500hPa ,  400 => 400hPa , 300 => 300hPa '
        print*,'                                                  '
        print*,'  250 =>  250hPa ,  200 => 200hPa , 150 => 150hPa '
        print*,'                                                  '
        print*,'  100 =>  100hPa ,   70 =>  70hPa ,  50 =>  50hPa '
        print*,'                                                  '
        print*,'   30 =>   30hPa ,   20 =>  20hPa ,  10 =>  10hPa '
        print*,'                                                  '
        print*,'    7 =>    7hPa ,    5 =>   5hPa ,   3 =>   3hPa '
        print*,'                                                  '
        print*,'    2 =>    2hPa ,    1 =>   1hPa , 999 =>    ALL '
        print*,'******************************************************'
        read(*,'(i4)') lev

        write(vlev,'(i4)') lev

        if (lev.ge.1000) write(vlev,'(i4)') lev
        if (lev.lt.1000) write(vlev,'(i3)') lev
        if (lev.lt.100)  write(vlev,'(i2)') lev
        if (lev.lt.10)   write(vlev,'(i1)') lev

        if (lev.ge.1000) nlev = 4
        if (lev.lt.1000) nlev = 3
        if (lev.lt.100)  nlev = 2
        if (lev.lt.10)   nlev = 1

        return
        end

c***********************************************************************
      subroutine find_date(isy,iey,inp,nh,ist,iet,csi,cas,sey)
c***********************************************************************
        character csi*3 , cas*3 , sey*9

        print*,'************************************************'
        print*,' Please Input 1st Year ( ex : 1957 ) '
        print*,'************************************************'
        read(*,'(i4)') isy
        print*,'************************************************'

        print*,'************************************************'
        print*,' Please Input 2nd Year ( ex : 2002 ) '
        print*,'************************************************'
        read(*,'(i4)') iey
        print*,'************************************************'

        print*,'************************************************'
        print*,' Input Data : '
        print*,' Monthly ( 0 ) or 6 hours ( 6 ) '
        print*,'************************************************'
        read(*,'(i1)') inp
        print*,'************************************************'

        print*,'************************************************'
        print*,' Output Data : '
        print*,' Monthly ( 0 ) or daily ( 1 ) or pentad ( 5 ) '
        print*,'************************************************'
        read(*,'(i1)') nh
        print*,'************************************************'

        print*,'************************************************'
        print*,' Output Data : '
        print*,' Please Input 1st Time Step '
        print*,'************************************************'
        read(*,'(i3)') ist
        print*,'************************************************'

        print*,'************************************************'
        print*,' Output Data : '
        print*,' Please Input 2nd Time Step '
        print*,'************************************************'
        read(*,'(i3)') iet
        print*,'************************************************'

        if (inp.eq.0) csi = 'mon'
        if (inp.ne.0) csi = '6hr'

        if (nh.eq.0) cas = 'mon'
        if (nh.eq.1) cas = 'day'
        if (nh.eq.5) cas = 'pnt'

        write(sey,'(i4,a1,i4)') isy,'_',iey

        return
        end

c***********************************************************************
      subroutine time_step(ity,itm,inp,nh,nn,nnm)
c***********************************************************************
        dimension mon(0:12)
        data mon / 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

        nn = nh*4

        if (mod(ity,4).ne.0) mon(2) = 28
        if (mod(ity,4).eq.0) mon(2) = 29

        if (nh.eq.0) nn = mon(itm)*4

        nnm = 0

        do iim = 1 , itm
           if (nh.eq.0) ivm = mon(iim-1)
           if ((nh.eq.0).and.(ity.eq.1957).and.(iim.le.9)) ivm = 0
           if (nh.eq.0) nnm = nnm + ivm
           if (nh.ne.0) nnm = iim - 1
           if ((nh.eq.1).and.(ity.eq.1957)) nnm = iim - 244
           if ((nh.eq.5).and.(ity.eq.1957)) nnm = iim - 49
           if ((nh.eq.5).and.(ity.eq.1957).and.(iim.eq.49)) nn = 8
        enddo

        if (nh.eq.0) nnm = nnm*4
        if (nh.eq.1) nnm = nnm*nn
        if ((nh.eq.1).and.(ity.eq.1957).and.(iim.le.244)) nnm = 0
        if ((nh.eq.5).and.(ity.ne.1957)) nnm = nnm*nn
        if ((nh.eq.5).and.(ity.eq.1957).and.(itm.ne.49)) nn = nh*4
        if ((nh.eq.5).and.(ity.eq.1957)) nnm = 8 + (nnm-1)*nn
        if ((nh.eq.5).and.(ity.eq.1957).and.(iim.le.49)) nnm = 0

        if (inp.eq.0) nn  = 1
        if (inp.eq.0) nnm = itm-1
        if ((inp.eq.0).and.(ity.eq.1957)) nnm = itm-9

c       print*,'nn = ',nn,'nnm = ',nnm

        return
        end
