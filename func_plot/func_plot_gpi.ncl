load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/home/pgchiu/walker1/func_read/func_misc.ncl"
load "func_dots.ncl"
load "func_plot_time_series.ncl" ;;plot_time_series(dts,title,filename)
load "func_read_rvke_sst.ncl"
load "func_gpi.ncl"
load "func_plot_trend.ncl"
load "func_dailygpi.ncl"
;;load "func_avginmons.ncl"

function plotsquare(wks,plot)
begin
    ;; bottom left corner start clockwise
        ;;ts = wgt_areaave(gpic(:,{10:25},{110:160}),1.,1.,0) ;; MGR
    linex = (/110.,110.,160.,160.,110./)
    liney = (/ 10., 25., 25., 10., 10./)
    mdrres = True
    mdrres@gsLineColor  ="red" 
    mdrres@gsLineThicknessF = 5.0    
    a = gsn_add_polyline(wks,plot,linex,liney,mdrres)
    return a
end
undef("annual_sum")
function annual_sum(monts[*],mon[*])
begin
    if(.not.any(mon.eq.0))then
        yrts = monts(::12)
        yrts = 0
        do m = 0,11
            mm = m+1
            if(any(mm.eq.mon))then
                yrts = yrts + monts(m::12)
            end if
        end do
        nm = dimsizes(mon)
    else
        yrts = month_to_annual(monts,0)
        nm = 12
    end if
    yrts@nmonth = nm
    return yrts
end
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

undef("read_gpic_datafile")
function read_gpic_datafile(years[*],var)
begin
    ;; var include: cETA, cRH, cVpot, cVshear, gpimax
    ;;              MqGPI, cMq
    mqvar = (/"mqGPI","cMq"/)
    if(isatt(years,"dataset"))then
        dataset = ""+years@dataset
    else
        dataset = ""
    end if
    if(dataset .eq."eraim")then
        dataset = ""
    end if

    ny = dimsizes(years)
    if(any(var.eq.mqvar))then
        f = addfile("gpi/mqgpi"+dataset+years(0)+".nc","r")
    else
        f = addfile("gpi/gpi"+dataset+years(0)+"comp.nc","r")
    end if
        
    gpiconeyear = f->$var$
    dims = dimsizes(gpiconeyear)
    dims(0) = ny*12
    gpicmons = new(dims,typeof(gpiconeyear))
    times   = new(ny*12,"integer")
    gpicmons!0 = "time"
    gpicmons!1 = gpiconeyear!1
    gpicmons&$gpicmons!1$ = gpiconeyear&$gpiconeyear!1$
    gpicmons!2 = gpiconeyear!2
    gpicmons&$gpicmons!2$ = gpiconeyear&$gpiconeyear!2$

    gpicmons(0:11,:,:) = (/gpiconeyear/)
    times(0:11) = yyyymm_time(years(0),years(0),"integer")

    do y = 1, ny-1
        if(any(var.eq.mqvar))then
            f = addfile("gpi/mqgpi"+dataset+years(y)+".nc","r")
        else
            f = addfile("gpi/gpi"+dataset+years(y)+"comp.nc","r")
        end if
        times(y*12:y*12+11) = yyyymm_time(years(y),years(y),"integer")
        gpicmons(y*12:y*12+11,:,:) = (/f->$var$/)
    end do

    gpicmons&time = times

    gpicmons@dataset = dataset
    return gpicmons
end

undef("read_gpi_datafile")
function read_gpi_datafile(years[*])
begin
    ny = dimsizes(years)
    if(isatt(years,"dataset"))then
        dataset = years@dataset
        print("data use: "+dataset)
    else
        dataset = ""
    end if
    f = addfile("gpi/gpi"+dataset+years(0)+".nc","r")
    gpioneyear = f->gpi
    dims = dimsizes(gpioneyear)
    dims(0) = ny*12
    gpimons = new(dims,typeof(gpioneyear))
    times   = new(ny*12,"integer")
    gpimons!0 = "time"
    gpimons!1 = gpioneyear!1
    gpimons&$gpimons!1$ = gpioneyear&$gpioneyear!1$
    gpimons!2 = gpioneyear!2
    gpimons&$gpimons!2$ = gpioneyear&$gpioneyear!2$

    gpimons(0:11,:,:) = (/gpioneyear/)
    times(0:11) = yyyymm_time(years(0),years(0),"integer")

    do y = 1, ny-1
        f = addfile("gpi/gpi"+dataset+years(y)+".nc","r")
        times(y*12:y*12+11) = yyyymm_time(years(y),years(y),"integer")
        gpimons(y*12:y*12+11,:,:) = (/f->gpi/)
    end do

    gpimons&time = times

    return gpimons
end

undef("read_gpiClm")
function read_gpiClm(years[*],mon)
begin
    if(isatt(years,"calclm").and.years@calclm)then
        gpi = cal_gpiclm(years,mon)
        return gpi
    else
        if(any(mon.eq.0))then
            mons = ispan(1,12,1)
        else 
            mons = mon
        end if
        gpimon = read_gpi_datafile(years)
        gpiclm = clmMonTLL(gpimon)
        gpiclm&month = ispan(1,12,1)
        if(dimsizes(mons).gt.1)then
            gpi = dim_avg_n_Wrap(gpiclm({mons},:,:),0)
        else
            gpi = gpiclm({mons},:,:)
        end if
        return gpi
    end if
    print("read_gpiClm() error")
    exit
end

undef("plot_gpi")
function plot_gpi(years,mon,title,filename)
begin
    gpi = read_gpiClm(years,mon) ;; read averaged gpi/month
    datamax = 8.
    if(isatt(years,"ano").and.years@ano)then
        load "res_years.ncl"
        gpiall = read_gpiClm(allyears,mon)
        gpi = gpi - gpiall
        datamax =  4.0
    end if
    gpisumlo = linint2_Wrap(gpi&lon,gpi&lat,gpi,False,gpi&lon(1::2),gpi&lat(1::2),0)

    res = True
    res@filename =  filename
    res@title = title
    
    if(any(mon.eq.0))then ;; allmon get lower GPI
        datamax = 4.0
    end if

    res@datamax =  datamax

    ;res@plotarea = "Pac"

    print(max(gpisumlo({0:40},{110:180})))

    res@overmark = max(abs(gpisumlo({0:40},{110:180})))/2
    a = plot_dot(gpisumlo,res)

    return True
end

undef("plot_dailygpi")
function plot_dailygpi(years,mon[*]:integer,title:string,filename:string)
begin
    m1 = (/1,32,60,91,121,152,182,213,244,274,305,335/) ;; 12 day of year
    m1 = m1-1
    m2 = (/31,59,90,120,151,181,212,243,273,304,334,365/)
    m2 = m2-1
    dgpi = cal_dailygpi(1995)
    dgpiclm = dim_avg_n_Wrap(dgpi,0)
    dgpiclm = 0.

    nm = dimsizes(mon)
    do m = 0,nm-1
        dmon = mon(m)-1
        dgpiclm = dgpiclm + dim_avg_n(dgpi(m1(dmon):m2(dmon),:,:),0)
    end do
    dgpiclm = dgpiclm/nm

    res = True
    res@filename =  filename
    res@title = title
    
    res@datamax =  8.

    lat  = fspan(-87.5,87.5,36)
    lon  = fspan(2.5,357.5,72)

    ;;dgpiclmlo = linint2_Wrap(dgpiclm&lon,dgpiclm&lat,dgpiclm,False,lon,lat,0)
        ;;dgpiclmlo!0 = "lat"
        ;;dgpiclmlo!1 = "lon"
        ;;dgpiclmlo&lat = lat
        ;;dgpiclmlo&lon = lon
    dgpiclmlo = dgpiclm(1::2,1::2)
    dgpiclmlo({17.5},{122.5}) = dgpiclmlo@_FillValue
    print(max(dgpiclmlo({0:40},{110:180})))
    res@overmark = max(abs(dgpiclmlo({0:40},{110:180})))/2
    a = plot_dot(dgpiclmlo({0:40},{110:180}),res)

    return True
end

undef("plot_gpic")
function plot_gpic(years[*],mon[*],var[1],title[1],filename[1])
begin
    ogpic = read_gpic_datafile(years,var)
    gpic = ogpic
    gpic = where(ogpic.eq.0,ogpic@_FillValue,ogpic)
    gpic!1 = "lat"
    gpic!2 = "lon"

    gpiclo = linint2_Wrap(gpic&lon,gpic&lat,gpic,False,gpic&lon(1::2),gpic&lat(1::2),0)
    gpiclo = log(gpiclo)/log(2.)  ;; log2 to emphasis result
    ;;gpiclo = log(gpiclo)
  
    gpicclm = clmMonTLL(gpiclo)

    if(.not.any(mon.eq.0))then
        do m = 0, 11
            mm = m+1
            if(.not.any(mm.eq.mon))then
                gpicclm(m,:,:) = gpicclm@_FillValue
            end if
        end do
    end if
    pgpicclm = dim_avg_n_Wrap(gpicclm,0)

    res = True
    res@filename =  filename
    res@title = title
    
    if(any(mon.eq.0))then
        nm = 12
    else
        nm = dimsizes(mon)
    end if

    res@datamax = 5.

    ;res@plotarea = "Pac"

    print("value from "+min(pgpicclm({0:40},{110:180}))+" to "+max(pgpicclm({0:40},{110:180})))
    res@overmark = max(abs(pgpicclm({0:40},{110:180})))/2
    if(filename .ne. "")then
        a = plot_dot(pgpicclm,res)
        return True
    end if

            ;;wks = gsn_open_wks("ps",var)
            ;;res = True
            ;;res@vpWidthF   = 0.80                          ; make map bigger
            ;;res@vpHeightF  = 0.80
            ;;res@mpMaxLatF  =  60.                           ; select subregion
            ;;res@mpMinLatF  = -60.
            ;;res@mpMinLonF  =  00.
            ;;res@mpMaxLonF  = 360.
            ;;res@mpCenterLonF = 180.
            ;;res@mpFillDrawOrder = "PreDraw"
            ;;res@cnLevelSpacingF  = 10
            ;;p = gsn_csm_contour_map_ce(wks,ogpic(7,:,:),res)

    return pgpicclm
end

undef("plot_gpic_plus")
function plot_gpic_plus(years[*],mon[*],vn1[1],vn2[1],title[1],filename[1])
begin
    var1 = plot_gpic(years,mon,vn1,"","")
    var2 = plot_gpic(years,mon,vn2,"","")

    pgpicclm = var1
    pgpicclm = var1+var2

    res = True
    res@filename =  filename
    res@title = title
    res@datamax = 4.
    res@overmark = max(abs(pgpicclm({0:40},{110:180})))/2
    a = plot_dot(pgpicclm,res)
    return a
end

undef("plot_gpic_ts")
function plot_gpic_ts(yb,ye,mon,var,title,filename)
begin
    ;title@RightMean  = True
    years = ispan(yb,ye,1)
    copy_VarAtts(yb,years)
    ym    = yyyymm_time(yb,ye,"integer")
    gpic = read_gpic_datafile(years,var)
    gpic = where(gpic.le.0,gpic@_FillValue,gpic)
    if(isatt(yb,"log").and.yb@log)then
        print("take log")
        gpic = log(gpic)/log(2.)
    end if

    if(isatt(yb,"range").and.yb@range.eq."MDR")then
        ts = wgt_areaave(gpic(:,{10:30},{110:150}),1.,1.,0) ;; MDR
    end if
    if(isatt(yb,"range").and.yb@range.eq."MGR")then
        ts = wgt_areaave(gpic(:,{10:25},{110:160}),1.,1.,0) ;; MGR
    end if
    if(.not.isdefined("ts"))then
        ts = wgt_areaave(gpic(:,{0:40},{110:180}),1.,1.,0)  ;; WNP
    end if
    if(isatt(yb,"ano"))then
        ano = yb@ano
    else
        ano = True
    end if
    ts!0 = "yyyymm"
    ts&yyyymm = ym
    ;;print(ts&$ts!0$+" "+ts)
    yrts = annual_sum(ts,mon)
    yrts = yrts/yrts@nmonth
    yrts!0 = "year"
    yrts&$yrts!0$ = years
    yrts@yref = 0. ;avg(yrts)
    if(ano)then
        yrts@ymax = 1.2 
        yrts@ymin = yrts@ymax*(-1)
        yrts = yrts - avg(yrts)
    else 
        yrts@ymax = 5.0 
        yrts@ymin = -3.0
    end if
    if(isatt(yb,"plotmode"))then
        yrts@plotmode = yb@plotmode
    end if
    p = plot_time_series(yrts,title,filename)
    return yrts
end

undef("plot_gpi_ts")
function plot_gpi_ts(yb,ye,mon,title,filename)
begin
    ;title@RightMean  = True
    years = ispan(yb,ye,1)
    copy_VarAtts(yb,years)
    gpi = read_gpi_datafile(years)
    gpisumlo = linint2_Wrap(gpi&lon,gpi&lat,gpi,False,gpi&lon(1::2),gpi&lat(1::2),0)
    ;;ts = dim_sum_Wrap(dim_sum_Wrap(gpisumlo(:,{0:40},{110:180})))
    ts = wgt_areaave(gpisumlo(:,{0:40},{110:180}),1.,1.,0)
    yrts = annual_sum(ts,mon)
    yrts = yrts / yrts@nmonth
    yrts!0 = "year"
    yrts&$yrts!0$ = years

    yrts@yref =    0.
    yrts@ymax =  6. ;4000.

    p = plot_time_series(yrts,title,filename)
    return yrts
end

undef("cor_gpi_race")
function cor_gpi_race(gpits,mons,var)
begin
    ;; cal correlation
    racename = (/"RACE","race","rvke","RvKE"/)
    acename  = (/"ACE" ,"ace","TyKE","tyke"/)
    tyhrname = (/"TyHr","hr"  ,"Hr","TLS","TCD","tcd"/)
    tygenname= (/"TyGen","Tygen","gen","TGN"/)
    pdiname  = (/"pdi","PDI"/)
    yb = min(gpits&year)
    if(isatt(gpits,"emanuel2005Correcting").and.gpits@emanuel2005Correcting)then
        yb@emanuel2005Correcting = True
    end if
    ye = max(gpits&year)

    if(any(var .eq. racename))then
        tyind = readmonRACE(yb,ye)
    end if
    if(any(var .eq. tygenname))then
        tyind = readmonTygen(yb,ye)
    end if
    if(any(var .eq. tyhrname))then
        tyind = readmonTyhr(yb,ye)
    end if
    if(any(var .eq. acename))then
        tyind = readmonACE(yb,ye)
    end if
    if(any(var .eq. pdiname))then
        print("read PDI not ready yet.")
        exit
        tyind = readmonACE(yb,ye)
    end if
    if(.not.isvar("tyind"))then
        print("no ty index named"+var)
        exit
    end if

    ;;tyindmonts = dim_sum_Wrap(dim_sum_Wrap(tyind(:,{0:40},{110:180})))
    tyindmonts = wgt_areaave_Wrap(tyind(:,{0:40},{110:180}),1.,1.,0)
    if(dimsizes(mons).eq.1 .and. mons.eq.0)then
    else
        do m = 0,11
            mm = m+1
            if(.not.any(mm.eq.mons))then
                tyindmonts(m::12) =  0 ;tyindmonts@_FillValue
            end if
        end do
    end if
    tyindyrts = month_to_annual(tyindmonts,0)
    cor = escorc(gpits,tyindyrts)

    return cor
end


undef("plot_race_all_gpits")
function plot_race_all_gpits(yb[1]:integer,ye[1]:integer,title[1]:string,filename[1]:string)
begin
    years = ispan(yb,ye,1)
    copy_VarAtts(yb,years)
    if(isatt(yb,"mons"))then
        mons = yb@mons
    else
        mons = ispan(1,12,1)
    end if
    if(isatt(yb,"range"))then
        range = yb@range
    else 
        range = "WNP"
    end if

    gen  = readmonTygen(yb,ye)
    hr   = readmonTyhr(yb,ye)
    ace  = readmonACE(yb,ye)
    race = readmonRACE(yb,ye)
    gpi  = read_gpi_datafile(years)

    if(range.eq."MGR")then
        genmonts  = wgt_areaave_Wrap( gen(:,{10:25},{110:160}),1.,1.,0)
        hrmonts   = wgt_areaave_Wrap(  hr(:,{10:25},{110:160}),1.,1.,0)
        acemonts  = wgt_areaave_Wrap( ace(:,{10:25},{110:160}),1.,1.,0)
        racemonts = wgt_areaave_Wrap(race(:,{10:25},{110:160}),1.,1.,0)
        gpimonts  = wgt_areaave_Wrap( gpi(:,{10:25},{110:160}),1.,1.,0)
    else
        genmonts  = wgt_areaave_Wrap( gen(:,{0:40},{110:180}),1.,1.,0)
        hrmonts   = wgt_areaave_Wrap(  hr(:,{0:40},{110:180}),1.,1.,0)
        acemonts  = wgt_areaave_Wrap( ace(:,{0:40},{110:180}),1.,1.,0)
        racemonts = wgt_areaave_Wrap(race(:,{0:40},{110:180}),1.,1.,0)
        gpimonts  = wgt_areaave_Wrap( gpi(:,{0:40},{110:180}),1.,1.,0)
    end if

    if(dimsizes(mons).eq.1 .and. mons.eq.0)then
        print("allmon")
    else
        do m = 0,11
            mm = m+1
            if(.not.any(mm.eq.mons))then
                genmonts(m::12)  =  0
                hrmonts(m::12)   =  0
                acemonts(m::12)  =  0
                racemonts(m::12) =  0
                gpimonts(m::12)  =  0
            end if
        end do
    end if
    genyrts  = month_to_annual(genmonts,0)
    hryrts   = month_to_annual(hrmonts,0)
    aceyrts  = month_to_annual(acemonts,0)
    raceyrts = month_to_annual(racemonts,0)
    gpiyrts  = month_to_annual(gpimonts,0)

    allyrts = new((/5,dimsizes(gpiyrts)/),typeof(genyrts))

    allyrts(0,:) = dim_standardize(genyrts,1)
    allyrts(1,:) = dim_standardize(hryrts,1)
    allyrts(2,:) = dim_standardize(aceyrts,1)
    allyrts(3,:) = dim_standardize(raceyrts,1)
    allyrts(4,:) = dim_standardize(gpiyrts,1)

    allyrts!0 = "tyind"
    allyrts!1 = "year"
    allyrts&year = years
    allyrts&tyind = (/"NTC","TCD","ACE","RACE","GPI"/)

    if(isatt(yb,"plotmode"))then
        allyrts@plotmode = yb@plotmode
    end if

    a = plot_allts(allyrts,title,filename)
    return allyrts
end

undef("plot_ntc_gpits")
function plot_ntc_gpits(yb[1]:integer,ye[1]:integer,title[1]:string,filename[1]:string)
begin
    years = ispan(yb,ye,1)
    copy_VarAtts(yb,years)
    if(isatt(yb,"mons"))then
        mons = yb@mons
    else
        mons = ispan(1,12,1)
    end if
    if(isatt(yb,"range"))then
        range = yb@range
    else 
        range = "WNP"
    end if
    if(isatt(yb,"nostand"))then
        nostand = True
    else
        nostand = False
    end if


    gen  = readmonTygen(yb,ye)
    gpi  = read_gpi_datafile(years)

    if(range.eq."MGR")then
        if(nostand)then
            genmonts  = dim_sum_n_Wrap( gen(:,{10:25},{110:160}),(/1,2/))
        end if
        genmonts  = wgt_areaave_Wrap( gen(:,{10:25},{110:160}),1.,1.,0)
        gpimonts  = wgt_areaave_Wrap( gpi(:,{10:25},{110:160}),1.,1.,0)
    else
        if(nostand)then
            genmonts  = dim_sum_n_Wrap( gen(:,{0:40},{110:180}),(/1,2/))
        else
            genmonts  = wgt_areaave_Wrap( gen(:,{0:40},{110:180}),1.,1.,0)
        end if
        gpimonts  = wgt_areaave_Wrap( gpi(:,{0:40},{110:180}),1.,1.,0)
    end if

    if(dimsizes(mons).eq.1 .and. mons.eq.0)then
    else
        do m = 0,11
            mm = m+1
            if(.not.any(mm.eq.mons))then
                genmonts(m::12)  =  0
                gpimonts(m::12)  =  0
            end if
        end do
    end if
    genyrts  = month_to_annual(genmonts,0)
    gpiyrts  = month_to_annual(gpimonts,0)

    allyrts = new((/2,dimsizes(gpiyrts)/),typeof(genyrts))

    if(nostand)then
        allyrts(0,:) = genyrts
        allyrts(1,:) = gpiyrts
    else
        allyrts(0,:) = dim_standardize(genyrts,1)
        allyrts(1,:) = dim_standardize(gpiyrts,1)
    end if

    allyrts!0 = "tyind"
    allyrts!1 = "year"
    allyrts&year = years
    allyrts&tyind = (/"NTC","GPI"/)

    if(isatt(yb,"plotmode"))then
        allyrts@plotmode = yb@plotmode
    end if

    a = plot_allts(allyrts,title,filename)
    return allyrts
end

undef("plot_gpists")
function plot_gpists(mons,title[1]:string,filename[1]:string)
begin
    years = ispan(1965,2008,1)
    eyears = ispan(1965,2001,1)
    cyears = ispan(1979,2008,1)  ;; 1979-2009

    gpi  = read_gpi_datafile(years)
    eyears@dataset = "era40"
    egpi = read_gpi_datafile(eyears)
    cyears@dataset = "cfsr"
    cgpi = read_gpi_datafile(cyears)

    gpimonts    = wgt_areaave_Wrap( gpi(:,{0:40},{110:180}),1.,1.,0)
    egpimonts   = wgt_areaave_Wrap(egpi(:,{0:40},{110:180}),1.,1.,0)
    cgpimonts   = wgt_areaave_Wrap(cgpi(:,{0:40},{110:180}),1.,1.,0)

    if(dimsizes(mons).eq.1 .and. mons.eq.0)then
    else
        do m = 0,11
            mm = m+1
            if(.not.any(mm.eq.mons))then
                gpimonts(m::12)   =  gpimonts@_FillValue
                egpimonts(m::12)  =  egpimonts@_FillValue
                cgpimonts(m::12)  =  cgpimonts@_FillValue
            end if
        end do
    end if
    gpiyrts  = month_to_annual(gpimonts,0)
    egpiyrts  = month_to_annual(egpimonts,0)
    cgpiyrts  = month_to_annual(cgpimonts,0)

    allyrts = new((/3,dimsizes(gpiyrts)/),typeof(gpiyrts))

    allyrts!0 = "GPI data"
    allyrts!1 = "year"
    allyrts&$allyrts!1$ = years
    allyrts&$allyrts!0$ = (/"ERAIM+ERA40","ERA40","CFSR"/)

    allyrts(0,:) = dim_standardize(gpiyrts,1)
    allyrts(1,{eyears}) = dim_standardize(egpiyrts,1)
    allyrts(2,{cyears}) = dim_standardize(cgpiyrts,1)

    allyrts(0,:) = (/gpiyrts/)
    allyrts(1,{eyears}) = (/egpiyrts/)
    allyrts(2,{cyears}) = (/cgpiyrts/)
    ;;print(allyrts(0,:)+" "+allyrts(1,:)+" "+allyrts(2,:))

    a = plot_allts(allyrts,title,filename)
    return allyrts
end

undef("plot_gpisSSTts")
function plot_gpisSSTts(mons,title[1]:string,filename[1]:string)
begin
    years = ispan(1965,2008,1)

    gpi  = read_gpi_datafile(years)
    years@dataset = "eraimsstP1"
    gpip1  = read_gpi_datafile(years)
    years@dataset = "eraimsstM1"
    gpim1  = read_gpi_datafile(years)
    years@dataset = "eraimsstP2"
    gpip2  = read_gpi_datafile(years)
    years@dataset = "eraimsstM2"
    gpim2  = read_gpi_datafile(years)


    gpimonts    = wgt_areaave_Wrap(  gpi(:,{0:40},{110:180}),1.,1.,0)
    gpim1monts  = wgt_areaave_Wrap(gpim1(:,{0:40},{110:180}),1.,1.,0)
    gpip1monts  = wgt_areaave_Wrap(gpip1(:,{0:40},{110:180}),1.,1.,0)
    gpim2monts  = wgt_areaave_Wrap(gpim2(:,{0:40},{110:180}),1.,1.,0)
    gpip2monts  = wgt_areaave_Wrap(gpip2(:,{0:40},{110:180}),1.,1.,0)

    if(dimsizes(mons).eq.1 .and. mons.eq.0)then
    else
        do m = 0,11
            mm = m+1
            if(.not.any(mm.eq.mons))then
                gpimonts(m::12)   =  gpimonts@_FillValue
                gpim1monts(m::12)   =  gpim1monts@_FillValue
                gpip1monts(m::12)   =  gpip1monts@_FillValue
                gpim2monts(m::12)   =  gpim2monts@_FillValue
                gpip2monts(m::12)   =  gpip2monts@_FillValue
            end if
        end do
    end if
    gpiyrts  = month_to_annual(gpimonts,0)
    gpim1yrts  = month_to_annual(gpim1monts,0)
    gpip1yrts  = month_to_annual(gpip1monts,0)
    gpim2yrts  = month_to_annual(gpim2monts,0)
    gpip2yrts  = month_to_annual(gpip2monts,0)

    allyrts = new((/5,dimsizes(gpiyrts)/),typeof(gpiyrts))

    allyrts!0 = "GPI data"
    allyrts!1 = "year"
    allyrts&$allyrts!1$ = years
    allyrts&$allyrts!0$ = (/"SST-2","SST-1","SST","SST+1","SST+2"/)

    allyrts(2,:) = (/gpiyrts/)
    allyrts(3,:) = (/gpip1yrts/)
    allyrts(1,:) = (/gpim1yrts/)
    allyrts(4,:) = (/gpip2yrts/)
    allyrts(0,:) = (/gpim2yrts/)

    a = plot_allts(allyrts,title,filename)
    return allyrts
end

undef("read_vpot_monly")
function read_vpot_monly(years[*]:integer)
begin
    ny = dimsizes(years)
    dataset = ""
    if(isatt(years,"dataset"))then
        dataset = years@dataset
    end if
    f = addfile("gpi/Vpot"+dataset+years(0)+".nc","r")
    vpotoneyear = f->Vpot
    dims = dimsizes(vpotoneyear)
    dims(0) = ny*12
    vpotmons = new(dims,typeof(vpotoneyear))
    times   = new(ny*12,"integer")
    lon = lonGlobeF(144, "lon", "longitude", "degrees_east")

    vpotmons!0 = "time"
    vpotmons!1 = vpotoneyear!1
    vpotmons&$vpotmons!1$ = vpotoneyear&$vpotoneyear!1$
    vpotmons!2 = vpotoneyear!2
    vpotmons&$vpotmons!2$ = lon ;vpotoneyear&$vpotoneyear!2$

    vpotmons(0:11,:,:) = (/vpotoneyear/)
    times(0:11) = yyyymm_time(years(0),years(0),"integer")

    do y = 1, ny-1
        f = addfile("gpi/Vpot"+dataset+years(y)+".nc","r")
        times(y*12:y*12+11) = yyyymm_time(years(y),years(y),"integer")
        vpotmons(y*12:y*12+11,:,:) = (/f->Vpot/)
    end do

    vpotmons&time = times

    return vpotmons
end
undef("read_vpot_clm")
function read_vpot_clm(years[*]:integer,imons[*]:integer)
begin
    if(any(imons.eq.0))then
        mons = ispan(1,12,1)
    else
        mons = imons
    end if
    vpotmonly = read_vpot_monly(years)
    vpotclm = clmMonTLL(vpotmonly)
    vpot = dim_avg_n_Wrap(vpotclm({mons-1},:,:),0)

    return vpot
end

undef("plot_vpot") ;; contour
function plot_vpot(years[*]:integer,mons[*]:integer,title[1]:string,filename[1]:string)
begin
    vpot = read_vpot_clm(years,mons)

    wks = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

        res                     = True
    
        res@vpHeightF            = 0.40           ; Changes the aspect ratio
        res@vpWidthF             = 0.85
        res@tmXTOn              = False
        res@tmYROn              = False
        res@mpMinLatF            = -40.
        res@mpMaxLatF            =  40.
        res@mpMinLonF            =   0.
        res@mpMaxLonF            = 360.
    
        res@mpCenterLonF         = 180.
    
        res@cnLinesOn           = True
        res@cnFillOn            = True            ; turn on color fill
        res@tiMainString        = title
    
        res@tiMainFont     = 29
        res@txFont         = 29
        res@tmXBLabelFont  = 29
        res@tmYLLabelFont  = 29
    
        res@gsnMajorLatSpacing  = 10
        res@gsnMajorLonSpacing  = 60
    
        res@tiMainFontHeightF       = 0.015
        res@gsnStringFontHeightF    = 0.012
        res@tmXBLabelFontHeightF    = 0.010
        res@tmYLLabelFontHeightF    = 0.010
        res@tmXBLabelFontHeightF    = 0.010
        res@tmYLLabelFontHeightF    = 0.010
    
        res@lbLabelBarOn             = True        ; turn on individual lb's
        res@lbLabelFont              = 21
        res@lbLabelFontHeightF       = 0.010
        res@pmLabelBarOrthogonalPosF = 0.2
        res@pmLabelBarWidthF         = 0.85
        res@pmLabelBarHeightF        = 0.04
    
        res@cnLevelSelectionMode  = "ExplicitLevels"
        ;res@cnLevels              = (/ 00,  35,  40,  45,  50,  55,  60,  65/)
        res@cnLevels              = ispan(35, 80,5)


        res@cnFillColors          = (/ 50,  62,  76,   83,  86,  92,  95,  97,  100,  102 ,102/)

        ;;    res2                      = True
        ;;    res2@cnLinesOn            = True
        ;;
        ;;    res2@cnLevelSelectionMode = "ManualLevels"
        ;;    res2@cnMinLevelValF       = 8
        ;;    res2@cnMaxLevelValF       = 18
        ;;    res2@cnLevelSpacingF      = 1           ; set contour spacing
        ;;
        ;;    res2@cnLineLabelFontHeightF = 0.006
        ;;    res2@cnLineThicknessF     = 2.
        ;;    res2@gsnContourZeroLineThicknessF = 0      ; eliminates zero contour
        ;;    res2@gsnContourNegLineDashPattern = 1
        ;;    res2@gsnContourZeroLineThicknessF = 0      ; eliminates zero contour
        ;;    res2@gsnContourNegLineDashPattern = 1
        ;;    res2@cnInfoLabelOn        = False
    plot = gsn_csm_contour_map_ce(wks,vpot,res)
    return vpot
end


;undef("plot_vpot_clm")
function plot_vpot_clm(years[*],imons[*],title[1]:string,filename[1]:string)
begin
    if(any(imons.eq.0))then
        mons = ispan(1,12,1)
        print("mons content 0, set to 1-12")
    else
        mons = imons
    end if
    if(isatt(years,"region"))then
        region = years@region
    else
        region = "WNP"
    end if
    if(isatt(years,"ano"))then
        ano = years@ano
    else
        ano = False
    end if
    if(isatt(years,"plotsst"))then
        plotsst = years@plotsst
    else
        plotsst = False
    end if

    if(plotsst)then
        SST12 = clmMonTLL(read_ersst(years))
        SST12&month = ispan(1,12,1)
        SST = dim_avg_n_Wrap(SST12(mons,:,:),0)
        delete(SST@long_name)
        delete(SST@units)
    else
        vpot = read_vpot_clm(years,mons)
    end if
    if(ano)then
        load "res_years.ncl"
        if(plotsst)then
            clmSST12 = clmMonTLL(read_ersst(allyears))
            clmSST12&month = ispan(1,12,1)
            clmSST = dim_avg_n_Wrap(clmSST12(mons,:,:),0)
            SST = SST - clmSST
        else
            clmvpot = read_vpot_clm(allyears,mons)
            vpot = vpot - clmvpot
        end if
    end if

    wks = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

    res                     = True
        res@gsnDraw = False
        res@gsnFrame = False
        res@vpHeightF            = 0.40           ; Changes the aspect ratio
        res@vpWidthF             = 0.85
        res@tmXTOn              = False
        res@tmYROn              = False
        if(region.eq."tropical")then  ;; default
            res@mpMinLatF            = -40.
            res@mpMaxLatF            =  40.
            res@mpMinLonF            =   0.
            res@mpMaxLonF            = 360.
        end if
        if(region .eq."WNP")then
            res@mpMinLatF            =   0.
            res@mpMaxLatF            =  40.
            res@mpMinLonF            = 110.
            res@mpMaxLonF            = 180.
        end if
    
        res@mpCenterLonF         = 180.
    
        res@cnLinesOn           = True
        res@cnFillOn            = True            ; turn on color fill
        res@tiMainString        = title
    
        res@tiMainFont     = 29
        res@txFont         = 29
        res@tmXBLabelFont  = 29
        res@tmYLLabelFont  = 29
    
        ;res@gsnMajorLatSpacing  = 10
        ;res@gsnMajorLonSpacing  = 60
    
        res@tiMainFontHeightF       = 0.015
        res@gsnStringFontHeightF    = 0.012
        res@tmXBLabelFontHeightF    = 0.010
        res@tmYLLabelFontHeightF    = 0.010
        res@tmXBLabelFontHeightF    = 0.010
        res@tmYLLabelFontHeightF    = 0.010
    
        res@lbLabelBarOn             = True        ; turn on individual lb's
        res@lbLabelFont              = 21
        res@lbLabelFontHeightF       = 0.010
        res@pmLabelBarOrthogonalPosF = 0.2
        res@pmLabelBarWidthF         = 0.85
        res@pmLabelBarHeightF        = 0.04
    

        res@cnFillDrawOrder  = "PreDraw"
        res@cnLevelSelectionMode  = "ExplicitLevels"
        ;res@cnLevels              = ispan(800,2400,200)
        res@cnLevels              = ispan(35, 80,5)
        res@cnFillColors          = (/ 50,  62,  76,   83,  86,  92,  95,  97,  100,  102 ,102/)
        if(ano)then
            delete(res@cnLevels)
            delete(res@cnFillColors)
            res@cnLevelSelectionMode = "ManualLevels"
            res@gsnSpreadColors  = True
            if(plotsst)then
                res@cnMinLevelValF   = -1.
                res@cnLevelSpacingF  =   .1
                res@cnMaxLevelValF   =  1.
            else
                res@cnMinLevelValF   = -10.
                res@cnLevelSpacingF  =   1.
                res@cnMaxLevelValF   =  10.
            end if
        end if

    print("plot: "+filename)
    if(plotsst)then
        plot = gsn_csm_contour_map_ce(wks,SST,res)
    else
        plot = gsn_csm_contour_map_ce(wks,vpot,res)
    end if
       ;a = plotsquare(wks,plot) ;; plot MGR square
    draw(plot)
    frame(wks)
    return plot
end

undef("plot_coloredScatter")
function plot_coloredScatter(x[*],y[*],c[*],title[1]:string,filename[1]:string)
begin
    res = True
    res@tiMainString      = title
    res@tiMainJust = "TopLeft"
    res@tiMainFont      = 21
    res@tiMainFontHeightF = 0.018
    res@tiMainPosition      = "Left"
    res@xyMarker  = 1
    res@xyMarkLineMode    = "Markers"  
    res@gsnFrame = False
    res@gsnDraw = False
    
    wks = gsn_open_wks("ps",filename)
    ;gsn_define_colormap(wks,"temp_19lev")
    gsn_define_colormap(wks,"so4_21")
    nc = 21 ; 
    print("plot "+filename)

    xmaxmin = nice_mnmxintvl(min(x),max(x),nc,True) ;; nc whatever
    ymaxmin = nice_mnmxintvl(min(y),max(y),nc,True)
    cmaxmin = nice_mnmxintvl(min(c),max(c),nc,True)
        ;; set x,y and color range if not setted
        xmin = xmaxmin(0)
        if(isatt(x,"min"))then
            xmin = x@min
        end if
        xmax = xmaxmin(1)
        if(isatt(x,"max"))then
            xmax = x@max
        end if
        ymin = ymaxmin(0)
        if(isatt(y,"min"))then
            ymin = y@min
        end if
        ymax = ymaxmin(1)
        if(isatt(y,"max"))then
            ymax = y@max
        end if
        cmin = cmaxmin(0)
        if(isatt(c,"min"))then
            cmin = c@min
        end if
        cmax = cmaxmin(1)
        if(isatt(c,"max"))then
            cmax = c@max
        end if
        cint = (cmax - cmin)/nc

    np = dimsizes(c)

    plot = gsn_csm_xy(wks,(/xmin,xmax/),(/ymin,ymax/),res)   

    pres =True
    pres@gsMarkerIndex = 16 
    do i =0, np-1
        if(.not.(ismissing(x(i)).or.ismissing(y(i)).or.ismissing(c(i))))then
            pres@gsMarkerColor  = floattoint((c(i)-cmin)/cint)+1
            if(pres@gsMarkerColor .lt. 2)then
                pres@gsMarkerColor = 2
            end if
            if(pres@gsMarkerColor .gt. nc)then
                pres@gsMarkerColor = nc
            end if
            gsn_polymarker(wks,plot,x(i),y(i),pres)
        end if
    end do
    lbres                           = True
    lbres@lbAutoManage              = False         ; we control label bar
    lbres@lbFillColors              = ispan(2,nc,1); use nice strided colors  
    lbres@lbPerimOn                 = False       ; Turn off labelbar perimeter.
    lbres@lbMonoFillPattern         = True          ; one pattern, all solid
    lbres@vpWidthF             = 0.1                 ; size
    lbres@vpHeightF            = 0.6
    lbres@lbLabelFontHeightF   = 0.012               ; label font height

    lblabels = fspan(cmin,cmax,nc)
    lbx = 0.8
    lby = 0.8
    gsn_labelbar_ndc(wks,nc,lblabels,lbx,lby,lbres) ;; warning message, but ok.
    ;;drawNDCGrid(wks)


    draw(plot)
    frame(wks)
    return True
end

undef("plot_cape_vpot") ;; scatter plot for senVpot*
function plot_cape_vpot(years[*]:integer,mons[*]:integer,title[1]:string,filename[1]:string)
begin
    dataset = ""
    if(isatt(years,"dataset"))then
        dataset = years@dataset
    end if
    range = "MDR"
    if(isatt(years,"range"))then
        range = years@range
    end if
    colorVar = "sst"
    if(isatt(years,"colorVar"))then
        colorVar = years@colorVar
    end if
    xaxis = "CAPE"
    if(isatt(years,"xaxis"))then
        xaxis = years@xaxis
    end if
    swapXS = False      ;; bad coding, but for quick plot...
    if(isatt(years,"swapXS"))then 
        swapXS = years@swapXS
    end if
    swapXYS2SXY = False
    if(isatt(years,"swapXYS2SXY"))then 
        swapXYS2SXY = years@swapXYS2SXY
    end if

    ny = dimsizes(years)
    nm = dimsizes(mons)

    selyymm = yyyymm_years(years)

    fns = "gpi/Vpot"+dataset+years+".nc"
    fs  = addfiles(fns,"r")
    ;capes = fs[:]->CAPEE({selyymm},:,:)
    capes = fs[:]->CAPEM({selyymm},:,:)
    vpots = fs[:]->Vpot({selyymm},:,:)
        capes = where(capes.gt.1e10,capes@_FillValue,capes)
    if(xaxis .eq. "dCAPE")then  ;; CAPE^\ast - CAPE|m
        scape = fs[:]->CAPE_ast({selyymm},:,:)
        scape = where(scape.gt.1e10,scape@_FillValue,scape)
        capes = scape - capes
    end if

    if(colorVar .eq."slp")then   ;; color var is sea level p
        vn = "P"
        cVar = read_ERAdataMonly(years,vn)  
        cVar  = cVar/100
        cVar@units = "hPa"
        cVar@max   =  1020.
        cVar@min   =  1000.
    end if

    if(colorVar .eq."sst")then   ;; color var is SST
        ;oSST = read_hadsst2(years)
        oSST = read_ersst(years)
        lon = vpots&$vpots!2$
        lat = vpots&$vpots!1$
        cVar = linint2_Wrap(oSST&lon,oSST&lat,oSST,True,lon,lat,0)
        cVar@max = 30.
        cVar@min = 20.
    end if

    ;; month filter
    do m = 0,11
        if(any(m+1.eq.mons))then
        else
            cVar(m::12,:,:) = cVar@_FillValue
            capes(m::12,:,:) = capes@_FillValue
            vpots(m::12,:,:) = vpots@_FillValue
        end if
    end do

    if(range.eq."MGR")then ;; MGR, Major Genesis Region (Chan and Liu 2004?)
        cape1 = ndtooned(capes(:,{10:25},{110:160}))
        vpot1 = ndtooned(vpots(:,{10:25},{110:160}))
        cVar1 = ndtooned( cVar(:,{10:25},{110:160}))
    end if
    if(.not.isvar("cape1"))then ;; MDR, Major Develop. Region
        cape1 = ndtooned(capes(:,{10:30},{110:150}))
        vpot1 = ndtooned(vpots(:,{10:30},{110:150}))
        cVar1 = ndtooned( cVar(:,{10:30},{110:150}))
    end if

    cVar1@max = cVar@max
    cVar1@min = cVar@min
    cape1@max = 12000.
    cape1@min =  0000.
    if(xaxis .eq."dCAPE")then
        cape1@max =  12000.
        cape1@min =   0000.
    end if
    vpot1@max =   120.
    vpot1@min =     0.

    if(swapXS)then
        b = plot_coloredScatter(cVar1,vpot1,cape1,title,filename)
        return b
    end if
    if(swapXYS2SXY)
        b = plot_coloredScatter(cVar1,cape1,vpot1,title,filename)
        return b
    end if

    b = plot_coloredScatter(cape1,vpot1,cVar1,title,filename)
    return b
    ;;test;;     vpot1(0:30) = fspan(10,110,31)
    ;;test;;     slp1 (0:30) = fspan(990,1010,31)
    ;;test;;     cape1(0:30) = 1000.
    ;;test;;     b = plot_coloredScatter(cape1(0:30),vpot1(0:30),slp1(0:30),title,filename)

end

undef("plot_vpot_sst_ts")
function plot_vpot_sst_ts(years[*]:integer,mons[*]:integer,title[1]:string,filename[1]:string)
begin
    dataset = ""
    if(isatt(years,"dataset"))then
        dataset = years@dataset
    end if
    range = "MDR"
    if(isatt(years,"range"))then
        range = years@range
    end if

    selyymm = yyyymm_years(years)
    fns = "gpi/Vpot"+dataset+years+".nc"
    fs  = addfiles(fns,"r")
    vpots = fs[:]->Vpot({selyymm},:,:)
    if(range.eq."MGR")then
        vpotts = wgt_areaave(vpots(:,{10:25},{110:160}),1.,1.,0) ;; MGR
    else
        vpotts = wgt_areaave(vpots(:,{10:30},{110:150}),1.,1.,0) ;; MDR
    end if
    vpotyrts = annual_sum(vpotts,mons)
    vpotyrts = vpotyrts/vpotyrts@nmonth

    capes = fs[:]->CAPEM({selyymm},:,:)
    scape = fs[:]->CAPE_ast({selyymm},:,:)
    capes@_FillValue = -999.99
    scape@_FillValue = -999.99
    capes = where(capes.gt.1e30,capes@_FillValue,capes)
    scape = where(scape.gt.1e30,scape@_FillValue,scape)
    
    ;capes = scape  - capes
    if(range.eq."MGR")then
        capemts = wgt_areaave(capes(:,{10:25},{110:160}),1.,1.,0) ;; MGR
        capests = wgt_areaave(scape(:,{10:25},{110:160}),1.,1.,0) ;; MGR
    else
        capemts = wgt_areaave(capes(:,{10:30},{110:150}),1.,1.,0) ;; MDR
        capests = wgt_areaave(scape(:,{10:30},{110:150}),1.,1.,0) ;; MDR
    end if
    capemyrts = annual_sum(capemts,mons)
    capesyrts = annual_sum(capests,mons)
    capemyrts = capemyrts/capemyrts@nmonth
    capesyrts = capesyrts/capesyrts@nmonth
    print(capemyrts+" "+capesyrts)

    ;oSST = read_hadsst2(years)
    oSST = read_ersst(years)
    lon = vpots&$vpots!2$
    lat = vpots&$vpots!1$
    sst = linint2_Wrap(oSST&lon,oSST&lat,oSST,True,lon,lat,0)
    if(range.eq."MGR")then
        sstts = wgt_areaave(sst(:,{10:25},{110:160}),1.,1.,0) ;; MGR
    else
        sstts = wgt_areaave(sst(:,{10:30},{110:150}),1.,1.,0) ;; MDR
    end if
    sstyrts = annual_sum(sstts,mons)
    sstyrts = sstyrts/sstyrts@nmonth

    allyrts = new((/5,dimsizes(vpotyrts)/),typeof(vpotyrts))
    allyrts!0 = "name"
    allyrts&name = (/"PI","SST","CAPEM","CAPE*","dCAPE"/)
    allyrts!1 = "time"
    allyrts&time = years
    allyrts(0,:) = (/vpotyrts/)
    allyrts(1,:) = (/sstyrts/)
    allyrts(2,:) = (/capemyrts/)
    allyrts(3,:) = (/capesyrts/)
    allyrts(4,:) = capesyrts - capemyrts
    allyrts = dim_standardize(allyrts,1)
    
    allyrts@ymin = -3.
    a = plot_allts(allyrts,title,filename)
    return allyrts
end

undef("plot_trend_map")
function plot_trend_map(iyears[*]:integer,imons[*]:integer,var[1]:string,title[1]:string,filename[1]:string)
begin
    ;; only min(years) and max(years) is taken.
    yb = min(iyears)
    ye = max(iyears)
    years = ispan(yb,ye,1)
    copy_VarAtts(iyears,years)
    ny = dimsizes(years)

    if(any(imons .eq. 0))then
        mons = ispan(1,12,1)
    else
        mons = imons
    end if

    dataset = ""
    if(isatt(iyears,"dataset"))then
        dataset = iyears@dataset
    end if
    
    ;;sstdt = dtrend_msg_n(ispan(1,ny,ny),sst,True,True,0)
    if(var.eq."SST")then
        ;oSST = read_hadsst2(years)
        oSST = read_ersst(years)
        ;lon = vpots&$vpots!2$
        ;lat = vpots&$vpots!1$
        lat  = latGlobeF(73 , "lat", "latitude", "degrees_north")
        lon  = lonGlobeF(144, "lon", "longitude", "degrees_east")
        sst = linint2_Wrap(oSST&lon,oSST&lat,oSST,True,lon,lat,0)
        sstyrts = annual_sum3d(sst,mons)
        sstyrts = sstyrts/sstyrts@nmonth
        yrts = sstyrts
        yrts@_FillValue = -999
        cmin = -2.5
        cmax =  2.5
        cint =  0.1
        sscale = 100.
    end if
    if(var.eq."RH")then
        ;; 700hPa relative humidity
        vn = "RH"
        vn@lev = 700
        if(dataset.eq."cfsr")then
            print("plot_trend_map() for CFSR not ready yet")
            exit
            tmpall = read_cfsr_monly(years,vn)
            RH    = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
            delete(tmpall)
        else
            RH = read_ERAdataMonly(years,vn)
        end if
        yrts = annual_sum3d(RH,mons)
        yrts = yrts/yrts@nmonth
        cmin = -5.0
        cmax =  5.0
        cint =  0.4
        sscale =10.
    end if
    if(var.eq."Vpot")then
        vpot  = read_vpot_monly(years)
        vpotyrts = annual_sum3d(vpot,mons)
        vpotyrts = vpotyrts/vpotyrts@nmonth
        yrts = vpotyrts
        yrts@_FillValue = -999
        cmin =-10.0
        cmax = 10.0
        cint =  0.1
        sscale = 10.
    end if
    if(var.eq."GPI")then
        gpi  = read_gpic_datafile(years,"cgpi")
        gpiyrts = annual_sum3d(gpi,mons)
        gpiyrts = gpiyrts/gpiyrts@nmonth
        yrts = gpiyrts
        yrts@_FillValue = -99.99
        cmin = -2.5
        cmax =  2.5
        cint =  0.1
        sscale = 10.
    end if
    if(var.eq."cGPI")then
        gpic  = read_gpic_datafile(years,"cgpi")
        gpic = where(gpic.le.0, gpic@_FillValue, gpic)
        gpic  = log(gpic)
        gpicyrts = annual_sum3d(gpic,mons)
        gpicyrts = gpicyrts/gpicyrts@nmonth
        yrts = gpicyrts
        yrts@_FillValue = -999
        cmin = -2.5
        cmax =  2.5
        cint =  0.1
        sscale = 50.
    end if
    if(var.eq."cVshear")then
        gpic  = read_gpic_datafile(years,"cVshear")
        gpic = where(gpic.le.0, gpic@_FillValue, gpic)
        gpic  = log(gpic)
        gpicyrts = annual_sum3d(gpic,mons)
        gpicyrts = gpicyrts/gpicyrts@nmonth
        yrts = gpicyrts
        yrts@_FillValue = -999
        cmin = -2.5
        cmax =  2.5
        cint =  0.1
        sscale = 50.
    end if
    if(var.eq."cETA")then
        gpic  = read_gpic_datafile(years,"cETA")
        gpic = where(gpic.le.0, gpic@_FillValue, gpic)
        gpic  = log(gpic)
        gpicyrts = annual_sum3d(gpic,mons)
        gpicyrts = gpicyrts/gpicyrts@nmonth
        yrts = gpicyrts
        yrts@_FillValue = -999
        cmin = -2.5
        cmax =  2.5
        cint =  0.1
        sscale = 50.
    end if
    if(var.eq."cVpot")then
        gpic  = read_gpic_datafile(years,"cVpot")
        gpic = where(gpic.le.0, gpic@_FillValue, gpic)
        gpic  = log(gpic)
        gpicyrts = annual_sum3d(gpic,mons)
        gpicyrts = gpicyrts/gpicyrts@nmonth
        yrts = gpicyrts
        yrts@_FillValue = -999
        cmin = -2.5
        cmax =  2.5
        cint =  0.1
        sscale = 50.
    end if
    if(var.eq."cRH")then
        gpic  = read_gpic_datafile(years,"cRH")
        gpic = where(gpic.le.0, gpic@_FillValue, gpic)
        gpic  = log(gpic)
        gpicyrts = annual_sum3d(gpic,mons)
        gpicyrts = gpicyrts/gpicyrts@nmonth
        yrts = gpicyrts
        yrts@_FillValue = -999
        cmin = -2.5
        cmax =  2.5
        cint =  0.1
        sscale = 50.
    end if
    title = title+" x"+sscale
    dtrendyr = dtrendR2(yrts)
    dtrendyr@sscale = sscale
    if(isvar("cmin"))then
        dtrendyr@cmin = cmin
    end if
    if(isvar("cmax"))then
        dtrendyr@cmax = cmax
    end if
    dtrendyr@cint = cint
    a = do_trend_plot(dtrendyr,title,filename)
    return a
end

;; test
    ;a = plot_cape_vpot(ispan(1965,1988,1),8,"CAPE-Vpot scatter mon="+8,"t"+1)
    ;do m =1,12
    ;  a = plot_cape_vpot(ispan(1965,2008,1),m,"CAPE-Vpot scatter mon="+m,"t"+m)
    ;end do
          ;begin
          ;;    load "res_years.ncl"
          ;;    ;a = read_gpic_datafile(allyears,"cETA")
          ;;    a = plot_gpic_ts(1965,2008,JJA,"cVshear","ts","ts")
          ;    a=plot_race_all_gpits(1965,2008,"t","t")
          ;a = plot_gpists(ispan(1,12,1),"t","t")
          ; a = plot_gpisSSTts(ispan(1,12,1),"GPI SST +- 1","GPI_SST_ts")
          ;end
    ;a = plot_dailygpi(1995,(/9,10,11/),"daily gpi","daily_gpi")
    ;yys = ispan(1998,2008,1)
    ;yys@plotsst = True
    ;yys@ano = True
    ;a = plot_vpot_clm(yys,(/9,10,11/),"test_vpot","vpot")
