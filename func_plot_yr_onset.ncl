load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "func_dailytopentad.ncl"
load "func_smsi.ncl"
;load "func_read_era_all.ncl"
load "func_read_reanalysis.ncl"
undef("plot_onset")
function plot_onset(u[*][*],v[*][*],rain[*],ores)
    begin
    ;; u,v is time-p,  rain is time(penta)
    title = ores@title
    filename = ores@filename
    if(isatt(ores,"race"))then
        race = ores@race
    end if

    times = u&$u!0$
    levs  = u&$u!1$
    nlev  = dimsizes(levs)
    nt    = dimsizes(times)

    smsindex = smsi(rain,u,v)

    res = True
    res@gsnFrame = False
    res@gsnDraw  = False
    res@vpXF           = 0.1
    res@vpHeightF      = 0.2               ; change aspect ratio of plot
    res@vpWidthF       = 0.7
    res@trXMaxF  = 73
    res@trXMinF  =  1
    res@trXMaxF  = 72
    res@trXMinF  = 25
    ;; fonts setup
    res@tiMainFont      = 21
    res@tiMainPosition      = "Left"
    res@tiMainJust      = "centerLeft"
    res@tiMainFontHeightF = 0.015
    res@tiXAxisFont = 21
    res@tiXAxisFontHeightF  = 0.015
    res@tiYAxisFont = 21
    res@tiYAxisFontHeightF  = 0.015
    res@tmXBLabelFont = 21
    res@tmXBLabelFontHeightF = 0.015
    res@tmYLLabelFont = 21
    res@tmYLLabelFontHeightF = 0.015
    res@tmYRLabelFont = 21
    res@tmYRLabelFontHeightF = 0.015
    ;; tickmarks
    res@tmXBAutoPrecision = False
    res@tmXBFormat = "f"
    res@tmXBMaxTicks = 18
    res@tmXBMinorOn = True
    res@tmXBMajorOutwardLengthF = 0.006
    res@tmXBMajorLengthF = 0.006
    res@tmXBMinorOutwardLengthF = 0.003
    res@tmXBMinorLengthF = 0.003
    res@tmYLMajorThicknessF  = 4.0
    res@tmYLMajorOutwardLengthF = 0.006
    res@tmYLMajorLengthF = 0.006
    res@tmYLMinorOutwardLengthF = 0.003
    res@tmYLMode = "Explicit"
    res@tmYLValues = fspan(0,nlev-1,nlev-1)
    res@tmYLLabels = levs
    ;; borders
    res@tmBorderThicknessF = 3.0

    res@tmYRMajorThicknessF  = 4.0
    res@tmYRMajorOutwardLengthF = 0.006
    res@tmYRMajorLengthF = 0.006
    res@tmYRMinorOutwardLengthF = 0.003

    xyres = res
    xyres@xyMarkLineMode = "MarkLines"
    xyres@xyMarker         =  16                      ; choose type of marker  
    xyres@xyLineThicknessF = 3.0
    xyres@xyLineColor      = "blue4"
    xyres@xyMarkerColor    = xyres@xyLineColor
    xyres@tmYLOn = False
    xyres@trYMaxF  = 20.
    xyres2= xyres
    xyres2@xyLineColor =  "transparent"
    xyres2@xyMarker         =   0                      ; choose type of marker  
    ;xyres2@xyMarkerColor    =  0
    xyres@xyCurveDrawOrder = "PreDraw"

    res@tiMainString = title
    res@vcGlyphStyle = "WindBarb"
    ;res@trYMaxF  = ind(levs.eq.100)
    res@trYMinF  = 0
    res@vcRefAnnoOn =  False
    res@vcRefLengthF            = 0.025   
    res@vcVectorDrawOrder = "PostDraw"
    res@tmYROn = False
    
    empty = rain
    empty = empty@_FillValue

    wks = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"prcp_3")
    plotbarb = gsn_vector(wks,u(lev|:,time|:),v(lev|:,time|:),res)
    plotxy   = gsn_csm_xy2(wks,times,rain,rain,xyres2,xyres)
    delete_VarAtts(xyres,"gsnFrame")
    delete_VarAtts(xyres,"gsnDraw")
    if(isatt(ores,"refline"))then
        refline = ores@refline
        xyres@gsLineColor = "red"
        xyres@gsLineThicknessF = 2.0     
        gsn_polyline(wks,plotxy,(/0,nt/),(/refline,refline/),xyres)
    end if
    if(isatt(ores,"fftrain"))then
        fftrain = fftfilter(rain,ores@fftrain)
        xyres@gsLineColor = "darkgreen"
        xyres@gsLineThicknessF = 3.0     
        xyres@gsLineDashPattern  = 2  
        gsn_polyline(wks,plotxy,times,fftrain,xyres)
    end if
    if(isatt(ores,"onset"))then
        fftrain = fftfilter(rain,ores@fftrain)
        xyres@gsLineColor = "goldenrod4"
        xyres@gsLineThicknessF = 3.0     
        xyres@gsLineDashPattern  = 1  
        onoffset = wnponset(rain,u,v)
        onset  = onoffset(0)+1
        offset = onoffset(1)+1
        if(any(ismissing(onoffset)))then
            print("find onset error")
        else
            gsn_polyline(wks,plotxy,(/onset,onset/),(/0,100/),xyres)
            ;gsn_polyline(wks,plotxy,(/offset,offset/),(/0,100/),xyres)
        end if
    end if
    if( isvar("race"))then
        xyres@gsLineColor = "mediumpurple4"
        xyres@gsLineThicknessF = 3.0     
        xyres@gsLineDashPattern  = 8  
        print(max(race))
        maxrace = 0.03
        maxrace = 1.80
        race = race/maxrace*xyres@trYMaxF
        gsn_polyline(wks,plotxy,times,race,xyres)
    end if
    overlay(plotbarb,plotxy)
    draw(plotbarb)
    frame(wks)
    return True
end

undef("plot_yr_onset")
function plot_yr_onset(year,filename)
begin
    ;; read u,v,p
    pentau = read_ncep_cmap_penta(year,"u")
    printVarSummary(pentau)
    pentav = read_ncep_cmap_penta(year,"v")
    pentar = read_ncep_cmap_penta(year,"r")
    res = True
    res@title = "WNP(120E-150E,5N-20N) Wind barb & CMAP "+year
    res@filename = filename
    res@refline  =  6
    res@fftrain  = 12
    res@onset    = True
    ;; add RACE data
    if (isatt(filename,"race"))then
        varname = "race"
        res@race = read_jtwc_race_penta(year,varname)
        print("add race line")
    end if

    p = plot_onset(pentau,pentav,pentar,res)
    return True
end

undef("read_ncep_cmap_penta")
function read_ncep_cmap_penta(year,varname)
begin
    datad = "../ncep_v2_daily/"
    if(isatt(varname,"odims"))then
        odims = varname@odims
    else
        odims = "TZ"
    end if

    if(varname.eq."u")then
        df = addfile(datad+"uwnd."+year+".nc","r")
        uwnd = df->uwnd
        dailyu = short2flt(uwnd(:,:,::-1,:))
        if (odims.eq."TZ")then
            dailyu = dailyu*1.94
            dailyu@units = "knots"
            pentau = yrdailytopentadTL(dim_avg_Wrap(dim_avg_Wrap(dailyu(:,:,{5:20},{120:150}))))
            pentau!1 = "lev"
            pentau&lev = uwnd&level
        end if
        if (odims.eq."TZYX")then
            pentau = yrdailytopentadTLn(dailyu(:,:,:,:))
            pentau!1 = "lev"
            pentau&lev = uwnd&level
            pentau!2 = "lat"
            pentau!3 = "lon"
            pentau&lat = dailyu&$dailyu!2$
            pentau&lon = dailyu&$dailyu!3$
        end if
        return pentau
    end if

    if(varname.eq."v")then
        df = addfile(datad+"vwnd."+year+".nc","r")
        vwnd = df->vwnd
        dailyv = short2flt(vwnd(:,:,::-1,:))
        if (odims.eq."TZ")then
            dailyv = dailyv*1.94
            dailyv@units = "knots"
            pentav = yrdailytopentadTL(dim_avg_Wrap(dim_avg_Wrap(dailyv(:,:,{5:20},{120:150}))))
            pentav!1 = "lev"
            pentav&lev = vwnd&level
        end if
        if (odims.eq."TZYX")then
            pentav = yrdailytopentadTLn(dailyv(:,:,:,:))
            pentav!1 = "lev"
            pentav&lev = vwnd&level
            pentav!2 = "lat"
            pentav!3 = "lon"
            pentav&lat = dailyv&$dailyv!2$
            pentav&lon = dailyv&$dailyv!3$
        end if
        return pentav
    end if

    if(varname.eq."r")then
        df = addfile("/walker1/pgchiu/cmap/precip.pentad.mean.nc","r")
        precip = df->precip
        ptime = ut_calendar(precip&time,-2) ; as YYYYMMDD
        precip&time = ptime
        if (odims.eq."TZ")then
            pentar = dim_avg_Wrap(dim_avg_Wrap(precip({year*10000:(year+1)*10000-1},{5:20},{120:150})))
        end if
        if (odims.eq."TZYX")then
            pentar = precip({year*10000:(year+1)*10000-1},:,:)
            cmaplat= pentar&lat
            pentar!1 = "cmaplat"
            pentar&cmaplat = cmaplat
        end if
        pentar&time = fspan(1,73,73)
        return pentar
    end if
end

undef("read_jtwc_race_penta")
function read_jtwc_race_penta(year,varname)
begin
    if(isatt(varname,"odims"))then
        odims = varname@odims
    else
        odims = "TZ"
    end if
    if(varname.eq."race")then
        df = addfile("/walker2/pgchiu/data/ncl/daily_RACE.nc","r")
        dailyRACE = df->race
        if (odims.eq."TZ")then
            ;pentaRACE = dailytopentad(dim_avg_Wrap(dim_avg_Wrap(dailyRACE({year*10000:(year+1)*10000-1},{5:20},{120:150}))))
            pentaRACE = dailytopentad(dim_sum_Wrap(dim_sum_Wrap(dailyRACE({year*10000:(year+1)*10000-1},:,:))))
        end if
        if (odims.eq."TYX")then
            pentaRACE = yrdailytopentadTLn(dailyRACE({year*10000:(year+1)*10000-1},:,:))
        end if
        return pentaRACE
    end if
end

undef("plot_clm_onset")
function plot_clm_onset(filename)
begin
    ;; read u,v,p
    yb = 1979
    ye = 2007
    ny = ye-yb+1

    pentau = read_ncep_cmap_penta(yb,"u")
    pentav = read_ncep_cmap_penta(yb,"v")
    pentar = read_ncep_cmap_penta(yb,"r")
    race   = read_jtwc_race_penta(yb,"race")
    do year = yb+1,ye
        pentau = pentau + read_ncep_cmap_penta(year,"u")
        pentav = pentav + read_ncep_cmap_penta(year,"v")
        pentar = pentar + read_ncep_cmap_penta(year,"r")
        race   = race   + read_jtwc_race_penta(year,"race")
    end do
    pentau = pentau /ny
    pentav = pentav /ny
    pentar = pentar /ny
    race   = race/ny

    res = True
    res@title = "WNP(120E-150E,5N-20N) Wind barb & CMAP cli"
    res@filename = filename
    res@refline  =  6
    res@fftrain  = 12
    res@race = race
    res@onset = True
    p = plot_onset(pentau,pentav,pentar,res)
    return True
end

undef("plot_compset_onset")
function plot_compset_onset(years,filename)
begin
    ;; read u,v,p
    ny = dimsizes(years)
    if(isatt(years,"title"))then
        title = years@title
    else
        title = ""
    end if
    if(isatt(filename,"title"))then
        title = filename@title
    end if

    pentau = read_ncep_cmap_penta(years(0),"u")
    pentav = read_ncep_cmap_penta(years(0),"v")
    pentar = read_ncep_cmap_penta(years(0),"r")
    race   = read_jtwc_race_penta(years(0),"race")
    do y = 1, ny-1
        pentau = pentau + read_ncep_cmap_penta(years(y),"u")
        pentav = pentav + read_ncep_cmap_penta(years(y),"v")
        pentar = pentar + read_ncep_cmap_penta(years(y),"r")
        race   = race   + read_jtwc_race_penta(years(y),"race")
    end do
    pentau = pentau /ny
    pentav = pentav /ny
    pentar = pentar /ny
    race   = race/ny

    res = True
    res@title = "WNP(120E-150E,5N-20N) Wind barb & CMAP "+title
    res@filename = filename
    res@refline  =  6
    res@fftrain  = 12
    res@race = race
    res@onset = True
    p = plot_onset(pentau,pentav,pentar,res)
    return True
end
