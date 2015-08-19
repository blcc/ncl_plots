load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "func_read_era_all.ncl"
load "func_read_reanalysis.ncl"
load "func_dailytopentad.ncl"

undef("plot_vpsf")
function  plot_vpsf(u[*][*],v[*][*],opt)
begin
    if(isatt(opt,"filename"))then
        filename = opt@filename
    else
        filename = "test_sfvp"
    end if
    if(isatt(opt,"plottype"))then
        plottype = opt@plottype
    else
        plottype = "vp" ; or "sf"
    end if
    if(isatt(opt,"title"))then
        title = opt@title
    else
        title = filename
    end if
    
    sfvp = uv2sfvpF(u,v)
    sfvp!1 = u!0
    sfvp!2 = u!1
    sfvp&$sfvp!1$ = u&$u!0$
    sfvp&$sfvp!2$ = u&$u!1$
    sf   = sfvp(0,:,:)
    vp   = sfvp(1,:,:)
    zsf  = sf
    zvp  = vp
    zsf  = 0
    zvp  = 0
    divu = zvp
    divv = zvp
    rotu = zsf
    rotv = zsf
    
    sfvp2uvf(sf,zvp,rotu,rotv)
    sfvp2uvf(zsf,vp,divu,divv)

    print("ploting: "+filename)
    wks = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"testcmap")

    res = True

    res@tiMainString         = title

    res@cnFillOn  = True
    res@cnLineLabelsOn       = False              ; turn off contour line labels
    res@cnLinesOn            = False
    res@gsnScalarContour     = True               ; contours desired
    res@gsnSpreadColors      = True               ; use full color map
    res@gsnSpreadColorStart  = 17                 ; start at color 17
    res@gsnSpreadColorEnd    = 200                ; end at color 200
    ;res@lbLabelBarOn = False

    res@mpMinLonF            =  90.               ; select a subregion
    res@mpMaxLonF            = 300.
    res@mpMinLatF            = -40. 
    res@mpMaxLatF            =  60.
    res@mpCenterLonF         = 180.
    res@mpGridAndLimbOn      = True
    res@mpGridLatSpacingF    = 60.
    res@mpGridLonSpacingF    = 180.
    res@mpGridLineColor      = "blue"
    res@mpFillOn             = False

    ;res@vcRefMagnitudeF           = 4.0             ; define vector ref mag
    res@vcRefLengthF              = 0.050           ; define length of vec ref
    res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
    res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
    res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
    ;res@vcGlyphStyle  = "FillArrow"
    ;res@vcGlyphStyle = "WindBarb"
    ;res@vcFillArrowEdgeColor  = -1
    res@vcLineArrowThicknessF  = 3.0
    res@vcLineArrowHeadMaxSizeF = 0.04

    ;plot = gsn_csm_vector_scalar_map_ce(wks,rotu,rotv,sf,res)
    if (plottype.eq."vp")then
        pu = divu
        pv = divv
        pfield = vp
    else
        pu = rotu
        pv = rotv
        pfield = sf
    end if
    velo    = sqrt(pu*pu+pv*pv)
    maxwind = max(velo)
    minwind = 0
    if(maxwind .lt. 10) then
        nicewind = 5
    else
        nicewind = 10
    end if
    if(maxwind .lt.  3) then
        nicewind = 2
    end if
    if(maxwind .lt.  2) then
        nicewind = 1
    end if
    res@vcRefMagnitudeF           = nicewind             ; define vector ref mag
    pfield = pfield/10000000
    plot = gsn_csm_vector_scalar_map_ce(wks,pu(::4,::4),pv(::4,::4),pfield,res)
    ;filename = filename+"_uv"
    ;print("ploting: "+filename)
    ;wks = gsn_open_wks("ps",filename)
    ;gsn_define_colormap(wks,"testcmap")
    ;plot = gsn_csm_vector_map_ce(wks,u,v,res)
    return True
end

undef("plot_uv_JJASON_vpsf")
function plot_uv_JJASON_vpsf(u[*][*][*],v[*][*][*],opt)
begin
    dims = dimsizes(u)
    nt   = dims(0)

    if(nt.eq.365)then
            nt1 = 151      ; for June 1 = 152th day of year
            nt2 = 343      ; for Nov.30 = 344th day of year
    end if
    if(nt.eq.366)then
            nt1 = 152      ; for June 1 = 153th day of year
            nt2 = 344      ; for Nov.30 = 345th day of year
    end if
    if(nt.eq.73)then
            nt1 = 30       ; where penta 31 = May31-Jun04
            nt2 = 67       ; where penta 68 = Nov27-Dec02 
    end if
    if(nt.eq.12)then
            nt1 = 5        ; June
            nt2 = 10       ; November
    end if

    if(.not.isvar("nt1"))then  ;; if not any of them
        print("unknow time step "+nt)
        exit
    end if

    u1 = dim_avg_n_Wrap(u(nt1:nt2,:,:),0)
    v1 = dim_avg_n_Wrap(v(nt1:nt2,:,:),0)

    a = plot_vpsf(u1,v1,opt)
    return a
end

undef("plot_JJASON_vpsf")
function plot_JJASON_vpsf(years[*],opt)
begin
    if (isatt(opt,"plottype"))then
        plottype = opt@plottype
    else
        plottype = "sf"
    end if
    if (isatt(opt,"lev"))then
        lev = opt@lev
    else
        lev = 850
    end if
    if (isatt(opt,"compset"))then
        compset = opt@compset
    else
        compset = "test_compset"
    end if
    if (isatt(opt,"filename"))then
        filename = opt@filename
    else
        filename = ""+plottype+"_"+compset+"_JJASON_"+lev+"hPa"
    end if
    if (isatt(opt,"title"))then
        title = opt@title
    else
        title = filename
    end if

    ny = dimsizes(years)
    uvar = "U"
    uvar@lev = lev
    vvar = uvar
    vvar = "V"
    compu = read_reanalysis_monclm(years,uvar)
    compv = read_reanalysis_monclm(years,vvar)

    opt@filename = filename
    opt@title    = title
    dimu = dimsizes(compu)
    a =plot_uv_JJASON_vpsf(compu,compv,opt)
    compuv = new((/2,dimu(0),dimu(1),dimu(2)/),"float")
    compuv!0 = "u_v"
    compuv(0,:,:,:) = compu
    compuv(1,:,:,:) = compv
    delete(compu)
    delete(compv)
    delete(filename)
    delete(opt@filename)
    return compuv
end


undef("plot_clm_vpsf")
function plot_clm_vpsf(years[*],mons[*],plottype,title,filename)
begin
    load "res_years.ncl"
    if (isatt(years,"lev"))then
        lev = years@lev
    else
        lev = 850
    end if
    if (isatt(years,"ano"))then
        ano = years@ano
    else
        ano = False
    end if
    ryears = years
    ryears@clm = True

    ny = dimsizes(years)
    uvar = "U"
    uvar@lev = lev
    vvar = uvar
    vvar = "V"
    compu = read_reanalysis_monclm(ryears,uvar)
    compv = read_reanalysis_monclm(ryears,vvar)
    u1 = dim_avg_n_Wrap(compu({mons},:,:),0)
    v1 = dim_avg_n_Wrap(compv({mons},:,:),0)
    if (ano)then
        clmu = read_reanalysis_monclm(allyears,uvar)
        clmv = read_reanalysis_monclm(allyears,vvar)
        clmus = dim_avg_n_Wrap(clmu({mons},:,:),0)
        clmvs = dim_avg_n_Wrap(clmv({mons},:,:),0)
        u1 = u1-clmus
        v1 = v1-clmvs
        delete(clmu)
        delete(clmv)
        delete(clmus)
        delete(clmvs)
    end if

    opt  = True
    opt@filename = filename
    opt@title    = title
    opt@plottype = plottype
    dimu = dimsizes(compu)
    a =plot_vpsf(u1,v1,opt)
    compuv = new((/2,dimu(0),dimu(1),dimu(2)/),"float")
    compuv!0 = "u_v"
    compuv(0,:,:,:) = compu
    compuv(1,:,:,:) = compv
    delete(compu)
    delete(compv)
    return compuv
end


undef("plot_clm_SONdJJA_vpsf")
function plot_clm_SONdJJA_vpsf(years[*],plottype,title,filename)
begin
    load "res_years.ncl"
    if (isatt(years,"lev"))then
        lev = years@lev
    else
        lev = 850
    end if
    ryears = years
    ryears@clm = True

    ny = dimsizes(years)
    uvar = "U"
    uvar@lev = lev
    vvar = uvar
    vvar = "V"
    compu = read_reanalysis_monclm(ryears,uvar)
    compv = read_reanalysis_monclm(ryears,vvar)
    u1 = dim_avg_n_Wrap(compu({SON},:,:),0)
    u1 = u1 - dim_avg_n_Wrap(compu({JJA},:,:),0)
    v1 = dim_avg_n_Wrap(compv({SON},:,:),0)
    v1 = v1 - dim_avg_n_Wrap(compv({JJA},:,:),0)

    opt  = True
    opt@filename = filename
    opt@title    = title
    opt@plottype = plottype
    dimu = dimsizes(compu)
    a =plot_vpsf(u1,v1,opt)
    compuv = new((/2,dimu(0),dimu(1),dimu(2)/),"float")
    compuv!0 = "u_v"
    compuv(0,:,:,:) = compu
    compuv(1,:,:,:) = compv
    delete(compu)
    delete(compv)
    return compuv
end


