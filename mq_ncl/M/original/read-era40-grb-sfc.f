c***********************************************************************
c     This program use to read ERA40 daily & monthly datasets
c
c     Yin-Min, Cho (2005/07/07)
c***********************************************************************
c***********************************************************************
c.. Surface level variable index
c***********************************************************************
c..  10  E    m of water   Evaporation (accum.)
c..  20  2T   K            2 m Temperature
c..  21  2D   K            2 m Dewpoint Temperature
c..  22  CP   m of water   Convective Precipitation (accum.)
c..  23  SD   m of water   Snow Depth
c..           equivalent
c..  24  SF   m of water   Snowfall (convective + stratiform) (accum.)
c..           equivalent
c..  25  RO   m of water   Runoff (accum.)
c..  30  10U  m/s          10 m U Wind Component
c..  31  10V  m/s          10 m V Wind Component
c..  32  TCW  kg/m**2      Total Column Water (liquid + ice + vapour)
c..  33  LSP  m of water   Large-Scale Precipitation (accum.)
c..  34  MSL  Pa           Mean Sea Level Pressure
c..  35  BLH  m            Boundary Layer Height
c..  36  SSR  W*s/m**2     Surface Solar Radiation (accum.)
c..  37  STR  W*s/m**2     Surface Thermal Radiation (accum.)
c..  38  TSR  W*s/m**2     Top Solar Radiation (accum.)
c..  39  TTR  W*s/m**2     Top Thermal Radiation (accum.)
c..  40  TCC  0-1          Total Cloud Cover
c..  41  LCC  0-1          Low Cloud Cover
c..  42  MCC  0-1          Medium Cloud Cover
c..  43  HCC  0-1          High Cloud Cover
c..  44  TSN  K            Temperature of Snow Layer
c..  50  TCWV kg/m**2      Total Column Water Vapor
c..  51  TCO3 kg/m**2      Total Column Ozone
c..  52  SSHF W*s/m**2     Surface Sensible Heat Flux (accum.)
c..  53  SLHF W*s/m**2     Surface Latent Heat Flux (accum.)
c..  54  SSRD W*s/m**2     Surface Solar Radiation Downwards (accum.)
c..  55  STRD W*s/m**2     Surface Thermal Radiation Downwards (accum.)
c..  56  EWSS N*s/m**2     East/West Surface Stress (accum.)
c..  57  NSSS N*s/m**2     North/South Surface Stress (accum.)
c..  58  LGWS N*s/m**2     Latitudinal Component of Gravity Wave Stress (accum.)
c..  59  MGWS N*s/m**2     Meridional Component of Gravity Wave Stress (accum.)
c..  60  TSRC W*s/m**2     Top Net Solar Radiation, Clear Sky (accum.)
c..  61  TTRC W*s/m**2     Top Net Thermal Radiation, Clear Sky (accum.)
c..  62  SSRC W*s/m**2     Surface Net Solar Radiation, Clear Sky (accum.)
c..  63  STRC W*s/m**2     Surface Net Thermal Radiation, Clear Sky (accum.)
c..  64  STL1 K            Soil Temperature Level 1
c..  65  STL2 K            Soil Temperature Level 2
c..  66  STL3 K            Soil Temperature Level 3
c..  67  STL4 K            Soil Temperature Level 4
c***********************************************************************
c***********************************************************************

      parameter ( nx=144  , ny=73  , spv=-99.99 )

      dimension dd(nx,ny) , ddm(nx,ny)

      character fn(5)*100 ,  year*4 ,  var*4   ,  csi*3
     &,          cas*3    ,  sey*9  ,  dir*100 ,  mmm*3
     &,          stp1*4   ,  stp2*4 ,  cnt*5

c***********************************************************************

      call read_dir(dir,ndir)

      call find_mac(nx,ny,nxy)

      call find_var(ivar,var,nnv,trans)

      call find_date(isy,iey,inp,nh,ist,iet,csi,cas,sey)

c***********************************************************************
      fn(2) = dir(1:ndir)//var(1:nnv)//'/ERA40_'//var(1:nnv)//
     &        '_SFC_'//sey//'_'//cas//'.bin'

      if (isy.eq.iey) fn(2) = dir(1:ndir)//var(1:nnv)//
     &         '/ERA40_'//var(1:nnv)//'_SFC_'//sey(1:4)//
     &         '_'//cas//'.bin'

      open (21,file=fn(2),status='unknown',form='unformatted'
     &                   ,access='direct',recl=nxy)
c***********************************************************************

      icnt = 0

      do 999 ity = isy , iey

         write (year,'(i4)') ity

         fn(1) = dir(1:ndir)//'ERA40_'//var(1:nnv)//'_SFC_'//
     &           year//'_'//csi//'.bin'

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
         if (((ivar.eq.20).or.(ivar.eq.21).or.(ivar.eq.44).or.
     &       (ivar.eq.64).or.(ivar.eq.65).or.(ivar.eq.66).or.
     &       (ivar.eq.67)).and.(dd(i,j).ne.spv))
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
        character var*4

        print*,'*******************************************************'
        print*,' Please Input Variables '
        print*,'*******************************************************'
        print*,'  10 ==>    E [m]                                    '
        print*,'  20 ==>   2T [K]         ,  21 ==>   2D [K]         '
        print*,'  22 ==>   CP [m]         ,  23 ==>   SD [m]         '
        print*,'  24 ==>   SF [m]         ,  25 ==>   RO [m]         '
        print*,'  30 ==>  10U [m/s]       ,  31 ==>  10V [m/s]       '
        print*,'  32 ==>  TCW [kg/m**2]   ,  33 ==>  LSP [m]         '
        print*,'  34 ==>  MSL [Pa]        ,  35 ==>  BLH [m]         '
        print*,'  36 ==>  SSR [W*s/m**2]  ,  37 ==>  STR [W*s/m**2]  '
        print*,'  38 ==>  TSR [W*s/m**2]  ,  39 ==>  TTR [W*s/m**2]  '
        print*,'  40 ==>  TCC [0-1]       ,  41 ==>  LCC [0-1]       '
        print*,'  42 ==>  MCC [0-1]       ,  43 ==>  HCC [0-1]       '
        print*,'  44 ==>  TSN [K]                                    '
        print*,'  50 ==> TCWV [kg/m**2]   ,  51 ==> TCO3 [kg/m**2]   '
        print*,'  52 ==> SSHF [W*s/m**2]  ,  53 ==> SLHF [W*s/m**2]  '
        print*,'  54 ==> SSRD [W*s/m**2]  ,  55 ==> STRD [W*s/m**2]  '
        print*,'  56 ==> EWSS [N*s/m**2]  ,  57 ==> NSSS [N*s/m**2]  '
        print*,'  58 ==> LGWS [N*s/m**2]  ,  59 ==> MGWS [N*s/m**2]  '
        print*,'  60 ==> TSRC [W*s/m**2]  ,  61 ==> TTRC [W*s/m**2]  '
        print*,'  62 ==> SSRC [W*s/m**2]  ,  63 ==> STRC [W*s/m**2]  '
        print*,'  64 ==> STL1 [K]         ,  65 ==> STL2 [K]         '
        print*,'  66 ==> STL3 [K]         ,  67 ==> STL4 [K]         '
        print*,'*******************************************************'
        read(*,'(i2)') ivar
        print*,'*******************************************************'

        trans = 1.0
        if (ivar.eq.10) trans = 1./(1000.*4.)
        if (ivar.eq.22) trans = 1./(1000.*4.)
        if (ivar.eq.23) trans = 1./(1000.)
        if (ivar.eq.24) trans = 1./(1000.*4.)
        if (ivar.eq.25) trans = 1./(1000.*4.)
        if (ivar.eq.33) trans = 1./(1000.*4.)
        if (ivar.eq.34) trans = 100.
        if (ivar.eq.36) trans = 6.*60.*60.
        if (ivar.eq.37) trans = 6.*60.*60.
        if (ivar.eq.38) trans = 6.*60.*60.
        if (ivar.eq.39) trans =-6.*60.*60.
        if (ivar.eq.51) trans = 2.1415*10.0**-5
        if (ivar.eq.52) trans =-6.*60.*60.
        if (ivar.eq.53) trans =-6.*60.*60.
        if (ivar.eq.54) trans = 6.*60.*60.
        if (ivar.eq.55) trans = 6.*60.*60.
        if (ivar.eq.56) trans = 6.*60.*60.
        if (ivar.eq.57) trans = 6.*60.*60.
        if (ivar.eq.58) trans = 6.*60.*60.
        if (ivar.eq.59) trans = 6.*60.*60.
        if (ivar.eq.60) trans = 6.*60.*60.
        if (ivar.eq.61) trans =-6.*60.*60.
        if (ivar.eq.62) trans = 6.*60.*60.
        if (ivar.eq.63) trans =-6.*60.*60.

        if (ivar.eq.10) var = 'E'
        if (ivar.eq.20) var = '2T'
        if (ivar.eq.21) var = '2D'
        if (ivar.eq.22) var = 'CP'
        if (ivar.eq.23) var = 'SD'
        if (ivar.eq.24) var = 'SF'
        if (ivar.eq.25) var = 'RO'
        if (ivar.eq.30) var = '10U'
        if (ivar.eq.31) var = '10V'
        if (ivar.eq.32) var = 'TCW'
        if (ivar.eq.33) var = 'LSP'
        if (ivar.eq.34) var = 'MSL'
        if (ivar.eq.35) var = 'BLH'
        if (ivar.eq.36) var = 'SSR'
        if (ivar.eq.37) var = 'STR'
        if (ivar.eq.38) var = 'TSR'
        if (ivar.eq.39) var = 'TTR'
        if (ivar.eq.40) var = 'TCC'
        if (ivar.eq.41) var = 'LCC'
        if (ivar.eq.42) var = 'MCC'
        if (ivar.eq.43) var = 'HCC'
        if (ivar.eq.44) var = 'TSN'
        if (ivar.eq.50) var = 'TCWV'
        if (ivar.eq.51) var = 'TCO3'
        if (ivar.eq.52) var = 'SSHF'
        if (ivar.eq.53) var = 'SLHF'
        if (ivar.eq.54) var = 'SSRD'
        if (ivar.eq.55) var = 'STRD'
        if (ivar.eq.56) var = 'EWSS'
        if (ivar.eq.57) var = 'NSSS'
        if (ivar.eq.58) var = 'LGWS'
        if (ivar.eq.59) var = 'MGWS'
        if (ivar.eq.60) var = 'TSRC'
        if (ivar.eq.61) var = 'TTRC'
        if (ivar.eq.62) var = 'SSRC'
        if (ivar.eq.63) var = 'STRC'
        if (ivar.eq.64) var = 'STL1'
        if (ivar.eq.65) var = 'STL2'
        if (ivar.eq.66) var = 'STL3'
        if (ivar.eq.67) var = 'STL4'

        if (ivar.le.70) nnv = 4
        if (ivar.le.44) nnv = 3
        if (ivar.le.25) nnv = 2
        if (ivar.le.10) nnv = 1

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
