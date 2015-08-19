load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;load "func_read_era_all.ncl"
load "func_read_reanalysis.ncl"
load "func_plot_trend.ncl"
load "func_plot_time_series.ncl"

;;undef("annual_sum_TLL") 
function annual_sum_TLL(monts[*][*][*],mon[*]) ;; monthly -> yearly
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
        yrts = month_to_annual(monts,0)
        nm = 12
    end if
    yrts@nmonth = nm
    return yrts
end

;;undef("clm_avg_TLL") 
function clm_avg_TLL(monts[*][*][*],mon[*]) ;; monthly -> one time(ex.JJA only)
begin
    dims = dimsizes(monts)
    clm = monts(0,:,:)
    clm = clm@_FillValue
    yrly = annual_sum_TLL(monts,mon)
    clm = dim_avg_n_Wrap(yrly,0)
    clm = clm/clm@nmonth
    delete(clm@nmonth)
    return clm
end

;;undef("plot_mflux") ;; moisture flux (obslated)
        ;;function plot_mflux(years[*],mons[*],title,filename)
        ;;begin
        ;;    if(isatt(years,"lev"))then
        ;;        lev = years@lev
        ;;    else
        ;;        lev = 850
        ;;    end if
        ;;
        ;;    ;; read u,v
        ;;    vn = "U"
        ;;    vn@lev = lev
        ;;    u = read_reanalysis_monclm(years,vn)
        ;;    vn = "V"
        ;;    v = read_reanalysis_monclm(years,vn)
        ;;
        ;;    ;; read specific humidity
        ;;    vn = "Q"
        ;;    q = read_reanalysis_monclm(years,vn)
        ;;
        ;;    ;; gradient q
        ;;    dxq = q
        ;;    dyq = q
        ;;    gradsf(q,dxq,dyq)
        ;;
        ;;    ;; read temperature
        ;;    vn = "T"
        ;;    temp = read_reanalysis_monclm(years,vn)
        ;;
        ;;    ;; density = P/RT
        ;;    R = 287.
        ;;    density = temp
        ;;    density = lev/R/temp
        ;;    
        ;;    ;; moisture flux = density*specific humidity gradent*V
        ;;    mu = u
        ;;    mu = density*dxq*u*1000 ;; let q units = g/kg
        ;;    mv = v
        ;;    mv = density*dyq*v*1000
        ;;    pmu = dim_avg_n_Wrap(mu({mons},:,:),0)
        ;;    pmv = dim_avg_n_Wrap(mv({mons},:,:),0)
        ;;    ;; moisture divergence
        ;;    pmflux = pmu
        ;;    pmflux = uv2dv_cfd(pmu,pmv,pmu&lat,pmu&lon,1)
        ;;    scale = 1.e07
        ;;    scale = 1.e10
        ;;    pmflux = pmflux *scale
        ;;
        ;;    ;; read vorticity 
        ;;    vn = "VOR"
        ;;    vn@lev = lev
        ;;    vor = read_reanalysis_monclm(years,vn)
        ;;
        ;;    scale = 1.e05
        ;;    vor = vor*scale
        ;;    pvor = dim_avg_n_Wrap(vor({mons},:,:),0)
        ;;
        ;;    ;; set plot res
        ;;    res = True
        ;;    res@gsnScalarContour     = True               ; contours desired
        ;;
        ;;    res@tiMainString = title
        ;;    res@gsnLeftString = ""
        ;;    res@gsnRightString = ""
        ;;    res@mpCenterLonF         = 180.               ; center plot at 180
        ;;    res@mpMinLonF            = 110.               ; select a subregion
        ;;    res@mpMaxLonF            = 180.
        ;;    res@mpMinLatF            =  00. 
        ;;    res@mpMaxLatF            =  40.
        ;;    res@mpFillOn = False
        ;;    res@cnLevelSelectionMode = "ManualLevels"
        ;;    res@cnMaxLevelValF    =  2. 
        ;;    res@cnMinLevelValF    = -2.
        ;;    res@cnLevelSpacingF   = .5
        ;;    ;res@cnMaxLevelValF    =  2. 
        ;;    ;res@cnMinLevelValF    = -2.
        ;;    ;res@cnLevelSpacingF   = .5
        ;;    res@gsnSpreadColors = True
        ;;    res@cnLinesOn = False
        ;;    res@cnFillOn             = True               ; turn on color for contours
        ;;    ;res@lbOrientation        = "Vertical"         ; vertical label bar
        ;;    ;res@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer
        ;;
        ;;    res@vcRefLengthF              = 0.055           ; define length of vec ref
        ;;    res@vcRefMagnitudeF = 1.
        ;;    res@vcRefAnnoArrowLineColor = "black"
        ;;    res@vcRefAnnoArrowUseVecColor = False
        ;;    res@vcGlyphStyle            = "CurlyVector"     ; turn on curley vectors
        ;;    res@vcLineArrowColor        = "white"           ; change vector color
        ;;    res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
        ;;
        ;;
        ;;      ;; anomaly
        ;;        ;if(isatt(years,"clm").and.(.not.years@clm))then
        ;;        ;end if
        ;;    ;; plot
        ;;    wks = gsn_open_wks("ps",filename) 
        ;;    gsn_define_colormap(wks,"BlWhRe")
        ;;
        ;;    ;plot = gsn_csm_vector_scalar_map_ce(wks,pmu,pmv,pvor,res)
        ;;    plot = gsn_csm_vector_scalar_map_ce(wks,pmu,pmv,pmflux,res)
        ;;    print("plot: "+filename)
        ;;
        ;;    return True
        ;;end

;;undef("read_mflux_monly")
function read_mflux_monly(years[*])
begin
    if(isatt(years,"lev"))then
        lev = years@lev
    else
        lev = 850
    end if
    only = "none"
    if(isatt(years,"only").and.years@only .eq."u")then
        only = "u"
    end if
    if(isatt(years,"only").and.years@only .eq."v")then
        only = "v"
    end if
    if(isatt(years,"only").and.years@only .eq."q")then
        ;; read specific humidity
        vn = "Q"
        vn@lev = lev
        q = read_reanalysis_monly(years,vn) ;; q units = kg/kg
        if(any(q@units.eq.(/"g/kg","g kg**-1"/)))then
            q = q/1000
        end if
        ;; moisture = specific humidity
        q = q*2500080 ;; Lwater = 2260000 by wikipedia
                      ;;  J/kg    2.50008*10^6 by book
        q@units = "J kg**-1"
        return q
    end if

    ;; read u,v for moist fluxergence
    if(only.ne."v")then
        vn = "U"
        vn@lev = lev
        u = read_reanalysis_monly(years,vn)
    end if
    if(only.ne."u")then
        vn = "V"
        vn@lev = lev
        v = read_reanalysis_monly(years,vn)
    end if

    ;; read specific humidity
    vn = "Q"
    q = read_reanalysis_monly(years,vn) ;; q units = kg/kg
    if(any(q@units.eq.(/"g/kg","g kg**-1"/)))then
        q = q/1000
    end if
    ;; moisture = specific humidity
    q = q*2500080 ;; Lwater = 2260000 by wikipedia
                  ;;  J/kg    2.50008*10^6 by book

    ;; moisture flux = \bigtriangledown \cdot V*q
    if(only.ne."v")then
        mu = u
        mu = q*u   
        mudx = u
        mudx = mudx@_FillValue
        mudy = u
        mudy = mudy@_FillValue
        gradsf(mu,mudx,mudy)
        mudx@units = "J kg**-1 s**-1"
        return mudx
    end if
    
    if(only.ne."u")then
        mv = v
        mv = q*v 
        mvdx = v
        mvdx = mvdx@_FillValue
        mvdy = v
        mvdy = mvdy@_FillValue
        gradsf(mv,mvdx,mvdy)
        mvdy@units = "J kg**-1 s**-1"
        return mvdy
    end if

    if(only.eq."none")then
        divVq = mudx + mvdy
        
        copy_VarCoords(mu,divVq)
        divVq@units = "J kg**-1 s**-1"
        return divVq
    end if
end

;;undef("read_mflux_monclm")
function read_mflux_monclm(years[*],mons[*])
begin
    mflux_monly = read_mflux_monly(years)
    mfluxclm = clmMonTLL(mflux_monly)

    pmflux = dim_avg_n_Wrap(mfluxclm({mons-1},:,:),0) ;; clmMonTLL() make month from 0 to 11
    return pmflux
end

;;undef("read_mflux_intP") ;; moisture fluxergence vertical integration
function read_mflux_intP(iyears[*],mons[*])
begin
    load "res_years.ncl"
    levs = (/ 100,150,200,250,300,400,500,600,700,775,850,925,1000 /) ;; from era lv_ISBL0
    levs   = levs(::-1)  ;; from bottom to top
    levspa = levs*100    ;; hPa to Pa
    nlev = dimsizes(levs)

    years = iyears
    years@lev = levs(0)
    mflux1 = read_mflux_monly(years)
    dims = dimsizes(mflux1)

    mflux_levs = new((/dims(0),nlev,dims(1),dims(2)/),typeof(mflux1)) ;; TLLL
    mflux_levs(:,0,:,:) = mflux1
    delete(mflux1)

    do l = 1, nlev-1
        years@lev = levs(l)
        mflux1 = read_mflux_monly(years)
        mflux_levs(:,l,:,:) = mflux1
        delete(mflux1)
    end do
    mflux_levs!1 = "level"
    mflux_levs&level = levspa

    vn = "P"
    Psfc = read_reanalysis_monly(years,vn)
    if(max(Psfc).lt.2000)then  ;; unit in hPa
        Psfc = Psfc*100   
        Psfc@units = "Pa"
    end if

    linlog = 2 ;; 1=linear 2=log
    pbot = 1000.  *100.   ;hPa to Pa
    ptop =  100.  *100.

    ;mflux_col = vibeta(levs,mflux_levs(time|:,lat|1:71,lon|:,level|:),linlog ,Psfc(:,1:71,:) ,pbot,ptop)
    mflux_col = vibeta(levspa,mflux_levs(time|:,lat|:,lon|:,level|:),linlog ,Psfc ,pbot,ptop)

    mflux_col = mflux_col/9.8

    copy_VarCoords(mflux_levs(time|:,lat|:,lon|:,level|0),mflux_col)

    if(isatt(years,"monly").and.years@monly)then
        return mflux_col
    end if

    mflux_col_ann = clm_avg_TLL(mflux_col,mons)

    mflux_col_ann@units = "W m**-2"

    return mflux_col_ann(:,:)
end
undef("read_mflux_intP_clm12")
function read_mflux_intP_clm12(only)
begin
    allyears = ispan(1965,2008,1)
    if(any(only.eq.(/"u","v"/)))then
        cachefile = "mflux_"+only+"_intP_clm12.nc"
    else
        cachefile = "mflux_intP_clm12.nc"
    end if
    if(isfilepresent(cachefile))then
        df = addfile(cachefile,"r")
        clm = df->mflux_intP_clm12
        return clm
    else
        print("creating cache file: "+cachefile)
        allyears@monly = True
        allyears@only = only
        monly = read_mflux_intP(allyears,ispan(1,12,1))
        clm = clmMonTLL(monly)
        clm&month = ispan(1,12,1)
        df = addfile(cachefile,"c")
        df->mflux_intP_clm12 = clm
        return clm
    end if
end

undef("read_E_P_watt")
function read_E_P_watt(years,EorP)
begin
    if(EorP.eq."E")then
        E = read_reanalysis_monly(years,"E")
        E = E*4000   ;; m/day to mm/day
        E = E*28     ;; to W/m^2
        E@units = "W/m^2"
        return E
    end if
    if(EorP.eq."P")then
        crain = read_reanalysis_monly(years,"crain")
        srain = read_reanalysis_monly(years,"srain")
        P = crain
        P = (crain+srain)
        P = P*4000   ;; m to mm/day
        P = P*28     ;; to W/m^2
        P@units = "W/m^2"
        return P
    end if
    print("var "+EorP+" not support")
    exit
end


;;undef("do_mflux_plot")
function do_mflux_plot(var_cn,var_sh,title,filename) ;; var for contour and shading
begin
    if(isatt(var_sh,"cmax"))then
        cmax = var_sh@cmax
    else
        cmmi = nice_mnmxintvl(min(var_sh),max(var_sh),10,True)
        cmax = max(abs(cmmi))
    end if
    if(isatt(var_sh,"cmin"))then
        cmin = var_sh@cmin
    else
        cmin = -1*cmax
    end if
    if(isatt(var_sh,"cint"))then
        cint = var_sh@cint
    else
        cint = cmax/10.
    end if

        ;;print("cmin/cmax/cint"+cmin+"/"+cmax+"/"+cint)
    if(isatt(var_cn,"cint"))then
        ccint = var_sh@cint*2
    else 
        ccint = cint*2
    end if
    ;; set plot res
    res = True   ;; shading
        ;res@tiMainString = title
        ;res@tiMainJust = "CenterLeft"
        ;res@tiMainPosition      = "Left"
        res@gsnLeftString = title
        res@gsnRightString = ""
        res@gsnDraw      = False
        res@gsnFrame     = False
        res@mpCenterLonF         = 180.               ; center plot at 180
        res@mpMinLonF            = 110.               ; select a subregion
        res@mpMaxLonF            = 180.
        res@mpMinLatF            =  00. 
        res@mpMaxLatF            =  40.
        res@cnLevelSelectionMode = "ManualLevels"
        res@cnMaxLevelValF    =  cmax
        res@cnMinLevelValF    =  cmin
        res@cnLevelSpacingF   =  cint
        if(isatt(var_sh,"ano"))then
            res@cnLevelSpacingF   =   0.1
            res@cnMaxLevelValF    =   0.5 
            res@cnMinLevelValF    =  -0.5
        end if
        res@gsnSpreadColors = True
        res@cnSmoothingOn = True
        res@cnSmoothingDistanceF = 0.005
        res@lbLabelStride  =10
        res@lbBoxLinesOn   = False

    res2 = True
        ;res2@cnLevelSelectionMode = "ManualLevels"
        ;res2@cnMaxLevelValF    = 100. 
        ;res2@cnMinLevelValF    =   0.
        res2@cnLevelSpacingF   =   ccint
        if(isatt(var_cn,"ano"))then
            res2@cnLevelSpacingF   =   .2
        end if
        res2@gsnContourZeroLineThicknessF = 5.0
        res2@gsnContourNegLineDashPattern = 1
        res2@cnLineThicknessF= 3.0
        res2@cnSmoothingOn = True
        res2@cnSmoothingDistanceF = 0.001
    
        res2@cnInfoLabelOn = False

    ;; plot
    wks = gsn_open_wks("ps",filename) 
    gsn_define_colormap(wks,"BlWhRe")
    gsn_reverse_colormap(wks)

    var_sh@long_name= ""      ;; prevent subtitle
    var_sh@units    = ""
    plot = gsn_csm_contour_map_overlay(wks,var_sh,var_cn,res,res2)
    draw(plot)
    frame(wks)
    print("plot: "+filename)

    return True
end

;;undef("plot_mflux") ;; moisture divergence & relative vorticity
function plot_mflux(years[*],mons[*],title,filename)
begin
    load "res_years.ncl"
    if(isatt(years,"lev"))then
        lev = years@lev
    else
        lev = 850
    end if
    if(isatt(years,"ano"))then
        ano = years@ano
    else
        ano = False
    end if

    pmflux = read_mflux_monclm(years,mons)
    scale = 1.e13
    pmflux = pmflux *scale
    if(ano)then
        pmflux = pmflux - scale*read_mflux_monclm(allyears,mons)
    end if

    ;; read vorticity
    vn = "VOR"
    vn@lev = lev
    vor = read_reanalysis_monclm(years,vn)
    pvor = dim_avg_n_Wrap(vor({mons},:,:),0)

    scale = 1.e05
    pvor = pvor*scale
    if(ano)then
        vorclm = read_reanalysis_monclm(allyears,vn)
        pvorclm = dim_avg_n_Wrap(vorclm({mons},:,:),0)
        ;print(max(pvor))
        ;print(max(pvorclm*scale))
        pvor = pvor - (scale*pvorclm)
    end if
    a = do_mflux_plot(pvor,pmflux,title,filename)
    return pmflux
end


;;undef("plot_mflux_trend")
function plot_mflux_trend(years[*],mons[*],title,filename)
begin
    mflux_monly = read_mflux_monly(years)
    mflux_yrly  = annual_sum_TLL(mflux_monly,mons)
    mflux_yrly  = mflux_yrly/mflux_yrly@nmonth
    mflux_trend = dtrendR2(mflux_yrly)
    mflux_trend@sscale = 10^15
    mflux_trend@cmin   = -3.
    mflux_trend@cmax   =  3.
    mflux_trend@cint   =   .1
    ;;a = do_mflux_plot(pmflux,pmflux,title,filename)
    a = do_trend_plot(mflux_trend,title,filename)
    return a
end

;;undef("plot_mflux_intP")
function plot_mflux_intP(years[*],mons[*],title,filename)
begin
    mflux = read_mflux_intP(years,mons)
    mflux = mflux
    print("mflux max: "+max(mflux))
    print("mflux max: "+min(mflux))
    print(mflux@units)
    mflux = mflux
    mflux@cmax = 500.
    mflux@cint =  20.
    a = do_mflux_plot(mflux,mflux,title,filename) ;; var for contour and shading
    return a
end

undef("plot_mflux_intP_trend")
function plot_mflux_intP_trend(iyears[*],mons[*],title,filename)
begin
    years = iyears
    years@monly = True  ;; read monthly data
    mflux_col = read_mflux_intP(years,mons)
    mflux_col_yrly = annual_sum_TLL(mflux_col,mons) ;; monthly -> yearly
    mflux_col_yrly = mflux_col_yrly/mflux_col_yrly@nmonth
    mflux_col_trend = dtrendR2(mflux_col_yrly)

    mflux_col_trend@cmin   = - 5.
    mflux_col_trend@cmax   =   5.
    mflux_col_trend@cint   =    .5

    a = do_trend_plot(mflux_col_trend,title,filename)
    return a
end

;;undef("plot_E_P") ;; Evaporation - Preciption
function plot_E_P(iyears[*],mons[*],title,filename)
begin
    E = read_E_P_watt(iyears,"E")
    P = read_E_P_watt(iyears,"P")

    Eclm = clm_avg_TLL(E,mons)
    Pclm = clm_avg_TLL(P,mons)
    EP = Eclm
    EP = (-Eclm-Pclm)
    print("Eclm max = "+max((Eclm({0:40},{110:180}))))
    print("Eclm min = "+min(Eclm({0:40},{110:180})))
    print("Pclm max = "+max((Pclm({0:40},{110:180}))))
    print("Pclm min = "+min(Pclm({0:40},{110:180})))
    print("E-P max  = "+max(abs(EP({0:40},{110:180}))))
    Pclm@cint = 10.
    Eclm@cmax = 500.
    Eclm@cint =  20
    EP@cmax = 500.
    EP@cint =  20.
    a = do_mflux_plot(EP,    EP,title,filename) 
    ;a = do_mflux_plot(Pclm,    Eclm,title,filename) 
                    ;;contour  shading
    return a
end

;;undef("plot_E_P_trend")
function plot_E_P_trend(years[*],mons[*],title,filename)
begin
    E = read_E_P_watt(years,"E")
    P = read_E_P_watt(years,"P")
    E_monly = E

    E_monly = -E - P   ;; E<0 => evaporation  E>0 => condisation

    E_yrly  = annual_sum_TLL(E_monly,mons)
    E_yrly  = E_yrly/E_yrly@nmonth
    E_trend = dtrendR2(E_yrly)
    E_trend@cmin   = - 5.
    E_trend@cmax   =   5.
    E_trend@cint   =   0.5
    ;;a = do_mflux_plot(pmflux,pmflux,title,filename)
    a = do_trend_plot(E_trend,title,filename)
    return a
end

undef("plot_E_P_ts")
function plot_E_P_ts(years,mons,title,filename)
begin
    E = read_E_P_watt(years,"E")
    P = read_E_P_watt(years,"P")
    E_monly = E
    E_monly = -E - P

    yrly  = annual_sum_TLL(E_monly,mons) 
    yrly  = yrly/yrly@nmonth
    yrts = wgt_areaave(yrly(:,{10:25},{110:150}),1.,1.,0)
    yrts!0 = "year"
    yrts&$yrts!0$ = years
    yrts@yref = 0. ;avg(yrts)
    yrts@ymax = 200 ; max(abs(yrts))
    yrts@ymin = 0 ;yrts@ymax*(-1)
    yrts = -yrts ;- avg(yrts)
    ;print(max(yrts))
    ;print(min(yrts))
    p = plot_time_series(yrts,title,filename)
    return p
end

;;undef("plot_mflux_SONdJJA") ;; moisture divergence & relative vorticity
function plot_mflux_SONdJJA(years[*],title,filename)
begin
    load "res_years.ncl"
    if(isatt(years,"lev"))then
        lev = years@lev
    else
        lev = 850
    end if

    pmfluxJJA = read_mflux_monclm(years,JJA)
    pmfluxSON = read_mflux_monclm(years,SON)
    pmflux    = pmfluxSON
    pmflux    = pmfluxSON - pmfluxJJA
    scale = 1.e07
    pmflux = pmflux *scale

    ;; read vorticity
    vn = "VOR"
    vn@lev = lev
    vor = read_reanalysis_monclm(years,vn)
    pvorJJA = dim_avg_n_Wrap(vor({JJA},:,:),0)
    pvorSON = dim_avg_n_Wrap(vor({SON},:,:),0)
    pvor = pvorSON
    pvor = pvorSON - pvorJJA

    scale = 1.e05
    pvor = pvor*scale

    ;; set plot res
    res = True
    res@tiMainString = title
    res@tiMainJust = "CenterLeft"
    res@tiMainPosition      = "Left"
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@gsnDraw      = False
    res@gsnFrame     = False
    res@mpCenterLonF         = 180.               ; center plot at 180
    res@mpMinLonF            = 110.               ; select a subregion
    res@mpMaxLonF            = 180.
    res@mpMinLatF            =  00. 
    res@mpMaxLatF            =  40.
    res@cnLevelSelectionMode = "ManualLevels"
        res@cnLevelSpacingF   =   0.1
        res@cnMaxLevelValF    =   0.5 
        res@cnMinLevelValF    =  -0.5
    res@gsnSpreadColors = True
    res@cnSmoothingOn = True
    res@cnSmoothingDistanceF = 0.005

    res2 = True
    ;res2@cnLevelSelectionMode = "ManualLevels"
    ;res2@cnMaxLevelValF    = 100. 
    ;res2@cnMinLevelValF    =   0.
        res2@cnLevelSpacingF   =   .5
    res2@gsnContourZeroLineThicknessF = 5.0
    res2@gsnContourNegLineDashPattern = 1
    res2@cnLineThicknessF= 3.0
    res2@cnSmoothingOn = True
    res2@cnSmoothingDistanceF = 0.001

    res2@cnInfoLabelOn = False

    ;; plot
    wks = gsn_open_wks("ps",filename) 
    gsn_define_colormap(wks,"BlWhRe")

    pvor@long_name= ""
    pvor@units    = ""
    ;prh@long_name = ""
    ;prh@units    = ""
    plot = gsn_csm_contour_map_overlay(wks,pvor,pmflux,res,res2)
    draw(plot)
    frame(wks)
    print("plot: "+filename)

    return True
end

undef("plot_mflux_intP_ts")
function plot_mflux_intP_ts(iyears,mons,title,filename)
begin
    years = iyears
    years@monly = True
    monly = read_mflux_intP(years,mons)
    yrly  = annual_sum_TLL(monly,mons) 
    yrly  = yrly/yrly@nmonth
    yrts = wgt_areaave(yrly(:,{10:25},{110:150}),1.,1.,0)
    yrts!0 = "year"
    yrts&$yrts!0$ = years
    yrts = -yrts ;- avg(yrts)
    yrts@yref = 0. ;avg(yrts)
    yrts@ymax = 200 ; max(abs(yrts))
    yrts@ymin = 0 ;yrts@ymax*(-1)
    ;print(max(yrts))
    ;print(min(yrts))
    p = plot_time_series(yrts,title,filename)
    return p
end

undef("plot_MoistBudget_ABCDE")
function plot_MoistBudget_ABCDE(term,iyears,mons,ano[1]:logical,title,filename)
begin
    ;; term list
    ;; A = \frac{\partial <q>}{\partial t}
    ;; B = -\frac{\partial}{\partial x} <uq>   years@only = "u"
    ;; C = -\frac{\partial}{\partial y} <vq>
    ;; D =  E (evaporation)
    ;; E = -P (precipitation)
    allyears = ispan(1965,2008,1)
    years = iyears
    if(ano)then
        cmax =  500.
        cint =   25.
    else
        cmax = 1000.
        cint =   20.
    end if
    if(isatt(iyears,"trend"))then
        trend = iyears@trend
    else
        trend = False
    end if

    if(term.eq."A")then
        ny = dimsizes(years)
        if(ny.ne.1)then
            ;print("plot_MoistBudget_ABCDE(): dq/dt should be 1 year")
            ;exit
        end if
        years@only = "q"
        years@monly = True
        ymp1 = years
        ymp1 = years-1
        q_ym1_12 = read_mflux_intP(ymp1,ispan(1,12,1))
        ymp1 = years+1
        q_yp1_12 = read_mflux_intP(ymp1,ispan(1,12,1))
        dq12    = (q_yp1_12-q_ym1_12)/(2*365*86400*ny)  ;; dif. between 2 years in month resolution
        if(trend)then
            dqyrly = annual_sum_TLL(dq12,mons)
            dqyrly = dqyrly/dqyrly@nmonth
            dqtrend = dtrendR2(dqyrly)
            copy_VarCoords(q_ym1_12(0:ny-1,:,:),dqtrend)
            dqtrend&$dqtrend!0$ = years
            ;dqtrend@cmin = 
            ;dqtrend@cmax = 
            ;dqtrend@cint = 
            ;dqtrend@@sscale = 10^15
            a = do_trend_plot(dqtrend,title,filename)
            return a
        end if
        dq      = clm_avg_TLL(dq12,mons)
        copy_VarCoords(q_ym1_12(0,:,:),dq)
        dq@cmax = cmax
        dq@cint = cint
        a = do_mflux_plot(dq, dq, title,filename)
        return dq
    end if

    if(term.eq."B")then
        years@monly = True
        years@only  = "u"
        mondmudx = read_mflux_intP(years,ispan(1,12,1))
        if(trend)then
            ny = dimsizes(years)
            yrly = annual_sum_TLL(mondmudx,mons)
            yrly = yrly/yrly@nmonth
            yrtrend = dtrendR2(yrly)
            copy_VarCoords(mondmudx(0:ny-1,:,:),yrtrend)
            a = do_trend_plot(yrtrend,title,filename)
            return a
        end if
        if(ano)then
            allyears@monly = years@monly
            allyears@only  = years@only
            clmdmudx = read_mflux_intP_clm12(allyears@only)
            dmudxP   = calcMonAnomTLL(mondmudx, clmdmudx)
            copy_VarCoords(mondmudx,dmudxP)
        else
            dmudxP    = mondmudx
        end if
        dmudx = clm_avg_TLL(dmudxP,mons)
        dmudx@cmax = cmax
        dmudx@cint = cint
        dmudx = -dmudx
        a = do_mflux_plot(dmudx, dmudx, title,filename)
        return dmudx
    end if

    if(term.eq."C")then
        years@monly = True
        years@only  = "v"
        mondmvdy = read_mflux_intP(years,ispan(1,12,1))
        if(trend)then
            ny = dimsizes(years)
            yrly = annual_sum_TLL(mondmvdy,mons)
            yrly = yrly/yrly@nmonth
            yrtrend = dtrendR2(yrly)
            copy_VarCoords(mondmvdy(0:ny-1,:,:),yrtrend)
            a = do_trend_plot(yrtrend,title,filename)
            return a
        end if
        if(ano)then
            allyears@monly = years@monly
            allyears@only  = years@only
            clmdmvdy = read_mflux_intP_clm12(allyears@only)
            dmvdyP   = calcMonAnomTLL(mondmvdy, clmdmvdy)
            copy_VarCoords(mondmvdy,dmvdyP)
        else
            dmvdyP    = mondmvdy
        end if
        dmvdy = clm_avg_TLL(dmvdyP,mons)
        dmvdy@cmax = cmax
        dmvdy@cint = cint
        dmvdy = -dmvdy
        a = do_mflux_plot(dmvdy, dmvdy, title,filename)
        return dmvdy
    end if

    if(term.eq."D")then
        Emonly = read_E_P_watt(years,"E")
        Emonly = -Emonly   ;; eva. <0 in data
        if(trend)then
            ny = dimsizes(years)
            yrly = annual_sum_TLL(Emonly,mons)
            yrly = yrly/yrly@nmonth
            yrtrend = dtrendR2(yrly)
            copy_VarCoords(Emonly(0:ny-1,:,:),yrtrend)
            a = do_trend_plot(yrtrend,title,filename)
            return a
        end if
        if(ano)then
            Eall = read_E_P_watt(allyears,"E")
            Eclm = clmMonTLL(Eall)
            Eclm&month = ispan(1,12,1)
            Eano = calcMonAnomTLL(Emonly, Eclm)
            copy_VarCoords(Emonly,Eano)
            E = clm_avg_TLL(Eano,mons)
        else
            E = clm_avg_TLL(Emonly,mons)
        end if
        E@cmax = cmax
        E@cint = cint
        a = do_mflux_plot(E, E, title,filename)
        return E
    end if

    if(term.eq."E")then
        Pmonly = read_E_P_watt(years,"P")
        Pmonly = -Pmonly    
        if(trend)then
            ny = dimsizes(years)
            yrly = annual_sum_TLL(Pmonly,mons)
            yrly = yrly/yrly@nmonth
            yrtrend = dtrendR2(yrly)
            copy_VarCoords(Pmonly(0:ny-1,:,:),yrtrend)
            a = do_trend_plot(yrtrend,title,filename)
            return a
        end if
        if(ano)then
            Pall = read_E_P_watt(allyears,"P")
            Pclm = clmMonTLL(Pall)
            Pclm&month = ispan(1,12,1)
            Pano = calcMonAnomTLL(Pmonly, Pclm)
            copy_VarCoords(Pmonly,Pano)
            P = clm_avg_TLL(Pano,mons)
        else
            P = clm_avg_TLL(Pmonly,mons)
        end if
        P@cmax = cmax
        P@cint = cint
        a = do_mflux_plot(P, P, title,filename)
        return P
    end if
end


;; test
;;a = plot_MoistBudget_ABCDE("A",(/1999/),(/6,7,8/),False,"test","test")
