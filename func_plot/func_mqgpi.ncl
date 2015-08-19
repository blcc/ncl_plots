load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "func_read_rvke_sst.ncl"
load "func_read_era_all.ncl"
load "func_read_mq.ncl"
load "func_plot_gpi.ncl"
load "func_dots.ncl"

undef("annual_sum3d")  ;; t,y,x
function annual_sum3d(monts[*][*][*],mon[*])
begin
    if(.not.any(mon.eq.0))then
        yrts = monts(::12,:,:)
        yrts = 0
        do m = 0,11
            mm = m+1
            if(any(mm.eq.mon))then
                yrts = yrts + monts(m::12,:,:)
            end if
        end do
        nm = dimsizes(mon)
    else
        yrts = month_to_annual(ts,0)
        nm = 12
    end if
    yrts@nmonth = nm
    return yrts
end

undef("cal_mqgpiclm12")
function cal_mqgpiclm12(year[1]) ;; Genesis potential index, Carmago et al. 2007
begin
    ;; 850hPa absolute vorticity
    vn = "VOR"
    vn@lev = 850
    VOR = read_ERAdataMonClm(year,vn)
    f = VOR
    dims = dimsizes(f)
    f = 0
    pi = 3.14159
    do j = 0, dims(1)-1
        f(:,j,:) = (2.) *0.000072921   * sin(f&lat(j)*pi/180)
              ;;    2    omega           sin(lat)
    end do
    ETA = VOR
    ETA = VOR + f

    ;; 700hPa relative humidity
    vn = "RH"
    vn@lev = 700
    RH = read_ERAdataMonClm(year,vn)

    ;; Mq as a replacement of tropical cyclone potential intensity
    mqfn = creat_mmq_cache(year)
    mqfile = addfile(mqfn,"r")
    Mq = mqfile->mq

    ;; 200-850 hPa wind shear magnitude
    vn = "U"
    vn@lev = 200
    Umon = read_ERAdataMonClm(year,vn)
    U = Umon
    vn@lev = 850
    Umonlo = read_ERAdataMonClm(year,vn)
    U = U - Umonlo
    vn = "V"
    vn@lev = 200
    Vmon = read_ERAdataMonClm(year,vn)
    V = Vmon
    vn@lev = 850
    Vmonlo = read_ERAdataMonClm(year,vn)
    V = V - Vmonlo
    Vshear = V
    Vshear = sqrt(U*U+V*V)

    ;; testing parameters
    mqpar = 70.
    mqind = 3.
    etaind = 1.5
    rhpar = 50.
    rhind = 3.
    vshearind = -2.

    if(isatt(year,"mqind"))then
        mqind = year@mqind
    end if
    if(isatt(year,"mqpar"))then
        mqpar = year@mqpar
    end if
    if(isatt(year,"etaind"))then
        etaind = year@etaind
    end if
    if(isatt(year,"rhind"))then
        rhind = year@rhind
    end if
    if(isatt(year,"rhpar"))then
        rhpar = year@rhpar
    end if
    if(isatt(year,"vshearind"))then
        vshearind = year@vshearind
    end if

    ;; GPI terms
    cETA = ETA
    cRH  = RH
    cMq  = Mq
    cMq  = sqrt(Mq)  ;; Mq^.5 unit: m/s
    cVshear = Vshear


    cETA = abs(100000*ETA)^etaind
    cRH  = (RH/rhpar)^rhind
    cMq  = (cMq/mqpar)^(mqind)
    cVshear = (1+0.1*Vshear)^(vshearind)

    ;; GPI
    mqGPI = cMq  ;; for coordinates
    mqGPI = mqGPI@_FillValue
    mqGPI = cETA(:,{-40:40},:) * cRH(:,{-40:40},:) * cMq(:,{-40:40},:) * cVshear(:,{-40:40},:)

    mqGPI@long_name = "Mq Genesis potential index"
    mqGPI@cETA = cETA(:,{-40:40},:)
    mqGPI@cRH  = cRH(:,{-40:40},:)
    mqGPI@cMq  = cMq(:,{-40:40},:)
    mqGPI@cVshear = cVshear(:,{-40:40},:)

    if(isatt(year,"savencfile"))then
        df = addfile("testmqgpi.nc","c")
        df->mqGPI = mqGPI
        df->cETA = cETA(:,{-40:40},:)
        df->cRH  = cRH(:,{-40:40},:)
        df->cMq  = cMq(:,{-40:40},:)
        df->cVshear = cVshear(:,{-40:40},:)
    end if

    ;; mask land
    oSST = read_ersst(2008)
    SST = linint2_Wrap(oSST&lon,oSST&lat,oSST,True,mqGPI&lon,mqGPI&lat,0)
    mqGPI = mask(mqGPI,ismissing(SST),False)
    return mqGPI
end

undef("creat_mqgpi_cache")
function creat_mqgpi_cache(year[1]:integer)
begin
    fn = "gpi/mqgpi"+year+".nc"
    if(isatt(year,"renew").and.year@renew)then
        system("rm -f "+fn)
        print("renew "+fn)
    else
        if(isfilepresent(fn))then
            return fn
        end if
    end if
    mqgpi = cal_mqgpiclm12(year)
    cETA  = mqgpi@cETA
    cRH   = mqgpi@cRH
    cMq   = mqgpi@cMq
    cVshear = mqgpi@cVshear
    copy_VarCoords(mqgpi,cETA)
    copy_VarCoords(mqgpi,cRH)
    copy_VarCoords(mqgpi,cMq)
    copy_VarCoords(mqgpi,cVshear)
      
    ff = addfile(fn,"c")
    ff->mqGPI   = mqgpi
    ff->cETA    = cETA
    ff->cRH     = cRH
    ff->cMq     = cMq
    ff->cVshear = cVshear

    return fn
end

undef("read_mqgpi_datafile")
function read_mqgpi_datafile(years[*]:integer)
begin
    ny = dimsizes(years)
    do y = 0,ny-1
        fn = creat_mqgpi_cache(years(y))
        df = addfile(fn,"r")
        mqgpi1 = df->mqGPI
        if(.not.isvar("mqgpi"))then
            dims = dimsizes(mqgpi1)
            dims(0) = ny*12
            mqgpi = new(dims,typeof(mqgpi1))
            mqgpi(0:11,:,:) = mqgpi1
        else
            mqgpi(y*12:y*12+11,:,:) = mqgpi1
        end if
        delete(mqgpi1)
    end do
    time = new(ny*12,"integer")
    do m = 1,12
        time(m-1::12) = years*100 + m
    end do
    mqgpi!0 = "time"
    mqgpi&time = time

    return mqgpi
end

undef("plot_race_all_mqgpits")
function plot_race_all_mqgpits(yb[1]:integer,ye[1]:integer,title[1]:string,filename[1]:string) 
begin
    years = ispan(yb,ye,1)
    copy_VarAtts(yb,years)
    if(isatt(yb,"mons"))then
        mons = yb@mons
    else
        mons = ispan(1,12,1)
    end if
    if(isatt(yb,"smoothing"))then
        smoothing = yb@smoothing
    else
        smoothing = "none"
    end if
    tycache = "tycache.nc"
    if(isfilepresent(tycache))then
        f = addfile(tycache,"r")
        gen = f->gen
        race = f->race
    else
        gen  = readmonTygen(yb,ye)
        ;hr   = readmonTyhr(yb,ye)
        ;ace  = readmonACE(yb,ye)
        race = readmonRACE(yb,ye)
        f = addfile(tycache,"c")
        f->gen = gen
        f->race = race
    end if

    mqgpi= read_mqgpi_datafile(years)
    gpi  = read_gpi_datafile(years)

    if(isatt(years,"region").and.years@region.eq."MGR")then
        print("Region: "+years@region)
        ;genmonts  = wgt_areaave_Wrap(  gen(:,{10:25},{110:160}),1.,1.,0)
        genmonts  = dim_sum_n_Wrap(  gen(:,{10:25},{110:160}),(/1,2/))
        ;hrmonts   = wgt_areaave_Wrap(   hr(:,{10:25},{110:160}),1.,1.,0)
        ;acemonts  = wgt_areaave_Wrap(  ace(:,{10:25},{110:160}),1.,1.,0)
        racemonts = wgt_areaave_Wrap( race(:,{10:25},{110:160}),1.,1.,0)
        mqgpimonts= wgt_areaave_Wrap(mqgpi(:,{10:25},{110:160}),1.,1.,0)
        gpimonts  = wgt_areaave_Wrap(  gpi(:,{10:25},{110:160}),1.,1.,0)
    else
        ;genmonts  = wgt_areaave_Wrap(  gen(:,{0:40},{110:180}),1.,1.,0)
        genmonts  = dim_sum_n_Wrap(  gen(:,{0:40},{110:180}),(/1,2/))
        ;hrmonts   = wgt_areaave_Wrap(   hr(:,{0:40},{110:180}),1.,1.,0)
        ;acemonts  = wgt_areaave_Wrap(  ace(:,{0:40},{110:180}),1.,1.,0)
        racemonts = wgt_areaave_Wrap( race(:,{0:40},{110:180}),1.,1.,0)
        mqgpimonts= wgt_areaave_Wrap(mqgpi(:,{0:40},{110:180}),1.,1.,0)
        gpimonts  = wgt_areaave_Wrap(  gpi(:,{0:40},{110:180}),1.,1.,0)
    end if


    if(dimsizes(mons).eq.1 .and. mons.eq.0)then
    else
        do m = 0,11
            mm = m+1
            if(.not.any(mm.eq.mons))then
                genmonts(m::12)  =  genmonts@_FillValue
                hrmonts(m::12)   =  hrmonts@_FillValue
                acemonts(m::12)  =  acemonts@_FillValue
                racemonts(m::12) =  racemonts@_FillValue
                mqgpimonts(m::12)=  mqgpimonts@_FillValue
                gpimonts(m::12)  =  gpimonts@_FillValue
            end if
        end do
    end if
    genyrts  = month_to_annual(genmonts,0)
    ;hryrts   = month_to_annual(hrmonts,0)
    ;aceyrts  = month_to_annual(acemonts,0)
    raceyrts = month_to_annual(racemonts,0)
    mqgpiyrts= month_to_annual(mqgpimonts,0)
    gpiyrts  = month_to_annual(gpimonts,0)

    allyrts = new((/4,dimsizes(mqgpiyrts)/),typeof(genyrts))

    ;print("mqgpi 1970 "+mqgpiyrts(5))
    ;print("  gpi 1970 "+  gpiyrts(5))
    ;print("  NTC 1970 "+  genyrts(5))
    corgen_mqgpi = escorc(genyrts,mqgpiyrts)
    print("cor NTC-mqgpi = "+corgen_mqgpi)

    sumerr = sum(abs(genyrts - mqgpiyrts))
    print("NTC-mqgpi error sum = "+sumerr)
    sumerr = sum(abs(mqgpiyrts - gpiyrts))
    print("gpi-mqgpi error sum = "+sumerr)
    sumerr = sum(abs(mqgpiyrts - raceyrts))
    print("RACE-mqgpi error sum = "+sumerr)

    print("gen stddev: "+stddev(genyrts))
    print("mqgpi stddev: "+stddev(mqgpiyrts))
    mgratio = stddev(genyrts)/stddev(mqgpiyrts)
    print("gen/mqgpiis about: "+mgratio)

    allyrts(0,:) = dim_standardize(genyrts,1)
    ;;allyrts(1,:) = dim_standardize(hryrts,1)
    ;;allyrts(2,:) = dim_standardize(aceyrts,1)
    allyrts(1,:) = dim_standardize(gpiyrts,1)
    allyrts(2,:) = dim_standardize(raceyrts,1)
    allyrts(3,:) = dim_standardize(mqgpiyrts,1)

    allyrts!0 = "tyind"
    allyrts!1 = "year"
    allyrts&year = years
    ;allyrts&tyind = (/"NTC","TCD","ACE","RACE","GPI"/)
    allyrts&tyind = (/"NTC","GPI","RACE","MqGPI"/)

    if(smoothing.eq."9yr")then
        allyrts = runave_n_Wrap(allyrts,9,0,1)
    end if
    if(isatt(yb,"plotmode"))then
        allyrts@plotmode = yb@plotmode
    end if

    a = plot_allts(allyrts,title,filename) ;; 
    return a

        ;;    wks = gsn_open_wks("ps",filename)
        ;;    res = True
        ;;
        ;;    res@vpWidthF   = 0.8
        ;;    res@vpHeightF   = 0.3
        ;;    res@tiMainFontHeightF = 0.018
        ;;    res@gsnLeftString = title
        ;;    res@gsnLeftStringFontHeightF = 21
        ;;
        ;;    res@pmLegendDisplayMode    = "Always" 
        ;;    res@pmLegendParallelPosF   = .65                   ; move units right
        ;;    res@pmLegendOrthogonalPosF = -0.5                  ; more neg = down
        ;;      
        ;;    res@pmLegendWidthF         = 0.10                  ; Change width and
        ;;    res@pmLegendHeightF        = 0.10                  ; height of legend.
        ;;    res@lgPerimOn              = False                 ; no box around
        ;;    res@lgLabelFontHeightF     = .01                   ; change font height
        ;;    res@xyExplicitLegendLabels = allyrts&tyind
        ;;
        ;;    res@trXMinF = yb
        ;;    res@trXMaxF = ye
        ;;
        ;;    res@xyLineThicknesses = (/ 2.0   , 2.0 , 2.0    , 3.0    /)
        ;;    res@xyLineColors      = (/"black","red", "green", "black"/)   
        ;;
        ;;    res@gsnYRefLine  = 0.
        ;;
        ;;
        ;;    print("plot xy lines: "+filename)
        ;;    plot  = gsn_csm_xy (wks,allyrts&year,allyrts,res)
        ;;
        ;;
        ;;    return allyrts
end

undef("plot_mqgpi_dots")
function plot_mqgpi_dots(years[*]:integer,mons[*]:integer,title[1]:string,filename[1]:string) 
begin
    mqgpi= read_mqgpi_datafile(years)
    ;; from 2.5x2.5 -> 5x5
    gpilo = linint2_Wrap(mqgpi&lon,mqgpi&lat,mqgpi,False,mqgpi&lon(1::2),mqgpi&lat(1::2),0)
    gpiloclm = clmMonTLL(gpilo)
    gpisumlo = annual_sum3d(gpiloclm,mons)
    gpisumlo = gpisumlo/gpisumlo@nmonth  ;; as average
    ;;printVarSummary(gpisumlo)
    print(max(gpisumlo(0,{0:40},{110:180})))

    res = True
    res@title = title
    res@filename = filename
    res@overmark = max(abs(gpisumlo(0,{0:40},{110:180})))/2
    res@datamax = 8.

    a =plot_dot(gpisumlo(0,:,:),res)

    return a
end


    ;; test
    ;a = cal_mqgpiclm12(1995)
    ;a = read_mqgpi_datafile(ispan(1965,2001,1))
    ;a = plot_race_all_mqgpits(1965,2008,"Standardized MqGPI, GPI and other metrics","mqgpi_ts")
    ;;printVarSummary(a)
