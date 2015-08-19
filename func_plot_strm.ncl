load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plot_strm_addpoint")
function plot_strm_addpoint(u[*][*],v[*][*],filename,xpts[*],ypts[*],ores)
begin
    res = ores
    res@gsnDraw    = False                         ; don't draw
    res@gsnFrame   = False                         ; don't advance frame
    res@mpFillOn = True
    res@mpMinLonF =   90
    res@mpMinLonF =  110
    res@mpMaxLonF =  180
    res@mpCenterLonF = 180.
    res@mpMaxLatF  =  40.                           ; select subregion
    res@mpMinLatF  =  00.
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@vcGlyphStyle = "CurlyVector"
    res@vcRefMagnitudeF = 2.
    res@vcRefLengthF     =  .05 
    res@tiMainString = filename

    pres = True
    pres@gsMarkerIndex = 16          ; Use filled dots for markers.
    pres@gsMarkerThicknessF = 3.
    pres@tfPolyDrawOrder = "PostDraw"
    print("Regenv: "+filename)

    wks  = gsn_open_wks("ps",filename)
    plot = gsn_csm_vector_map_ce(wks,u,v,res)
    ;plot = gsn_csm_streamline_map_ce(wks,u,v,res)

    a = gsn_add_polymarker(wks,plot,xpts,ypts,pres)
    draw(plot)
    frame(wks)

    return True
end

undef("plot_strm_z_addpoint")
function plot_strm_z_addpoint(u[*][*],v[*][*],z[*][*],shadlev,filename,xpts[*],ypts[*],ores)
begin
    res = ores
    res@gsnDraw    = False                         ; don't draw
    res@gsnFrame   = False                         ; don't advance frame
    res@mpFillOn = True
    res@mpMinLonF =  110
    res@mpMaxLonF =  180
    res@mpCenterLonF = 180.
    res@mpMaxLatF  =  40.                           ; select subregion
    res@mpMinLatF  =   0.
    res@mpFillDrawOrder = "PreDraw"
    res@mpOutlineOn            = True

    res@gsnLeftString = ""
    res@gsnRightString = ""
    if (.not.isatt(ores,"vcRefMagnitudeF"))then
        res@vcRefMagnitudeF = 15.
    end if
    res@vcRefLengthF     =  .015
    ;res@vcRefAnnoOn  = False
    res@tiMainString = "(A) "+filename
    res@tiMainJust = "CenterLeft"
    res@tiMainPosition = "Left"
    res@vcGlyphStyle = "LineArrow"
    res@vcGlyphStyle = "CurlyVector"
    res@vcRefLengthF     =  .075
    ;res@vcGlyphStyle = "Windbarb"
    ;u = u*1.943844   ; m/s to knots
    ;v = v*1.943844

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels     = shadlev
    res@cnFillColors = (/"transparent","yellow","yellowgreen"/)
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnScalarContour     = True  
    res@cnExplicitLegendLabelsOn = True
    ;res@lbOrientation = "Vertical"  
    res@pmLabelBarOrthogonalPosF  = -0.2
    res@pmLabelBarParallelPosF  = 0.4

    pres = True
    pres@gsMarkerIndex = 16          ; Use filled dots for markers.
    pres@gsMarkerThicknessF = 3.
    pres@tfPolyDrawOrder = "PostDraw"
    print("Regenv: "+filename)

    wks  = gsn_open_wks("ps",filename)
    delete(res@lev)
    if(ores@lev.eq.200)then
        ;plot = gsn_csm_vector_map_ce(wks,u,v,res)
        plot = gsn_csm_vector_map_ce(wks,u,v,res)
    else
        plot = gsn_csm_vector_scalar_map_ce(wks,u,v,z,res)
    end if

    a = gsn_add_polymarker(wks,plot,xpts,ypts,pres)
    ;; add map coast line
    draw(plot)
    frame(wks)

    return True
end

undef("plot_strm_rain")
function plot_strm_rain(u[*][*],v[*][*],r[*][*],filename,title)
begin
    res = True
    res@gsnDraw    = False                         ; don't draw
    res@gsnFrame   = False                         ; don't advance frame
    res@mpFillOn = False
    res@mpMinLonF =   90.
    res@mpMinLonF =  110.
    res@mpMaxLonF =  180.
    res@mpCenterLonF = 180.
    res@mpMaxLatF  =  40.                           ; select subregion
    res@mpMinLatF  =  00.
    res@mpGeophysicalLineThicknessF = 2.0
    ;res@mpFillDrawOrder = "PreDraw"
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@tiMainString = title
    res@tiMainJust = "CenterLeft"
    res@tiMainPosition      = "Left"


    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLevelSelectionMode        = "ManualLevels"       ; manual levels
    ;res@lbOrientation = "Vertical"  
    res@stLineStartStride = 2
    res@stArrowLengthF = .01
    print("Regenv: "+filename)

    ;; to mm
    rd = r
    rd = r*1000*4
    wks  = gsn_open_wks("ps",filename)
    if(isatt(u,"clm") .and. u@clm)then
        res@cnMaxLevelValF = 12.
        res@cnMinLevelValF =  6.
        res@cnLevelSpacingF             = 1.
        gsn_define_colormap(wks,"precip2_17lev")
    else
        res@cnMaxLevelValF =   5.
        res@cnMinLevelValF =  -5.
        res@cnLevelSpacingF  = 1.
        res@gsnSpreadColors  = True
        gsn_define_colormap(wks,"precip4_diff_19lev")
    end if
    if(isatt(u,"lev") .and. u@lev.eq.850)then
        plot = gsn_csm_streamline_contour_map_ce(wks,u,v,rd,res)
    else
        vv = u
        vv = sqrt(u*u+v*v)
        res@cnFillOn = False
        res@cnLinesOn = True
        plot = gsn_csm_streamline_contour_map_ce(wks,u,v,vv,res)
        ;plot = gsn_csm_streamline_map_ce(wks,u,v,res)
    end if

    draw(plot)
    frame(wks)

    return True
end

undef("plot_strm")
function plot_strm(years[*],mons[*],filename,title,opt)
begin
    load "res_years.ncl"
    nm = dimsizes(mons)
    ny = dimsizes(years)
    if(isatt(opt,"lev"))then
        lev = opt@lev
    else
        lev = 850
    end if

    vn = "U"
    vn@lev = lev
    u = read_reanalysis_monly(years(0),vn)
    vn = "V"
    v = read_reanalysis_monly(years(0),vn)
    vn@lev = lev
    vn = "crain"
    cr = read_reanalysis_monly(years(0),vn)
    vn = "srain"
    sr = read_reanalysis_monly(years(0),vn)
    rain = cr
    rain = cr + sr
    do y = 1,ny-1
        vn = "U"
        vn@lev = lev
        u = u + read_reanalysis_monly(years(y),vn)
        vn = "V"
        v = v + read_reanalysis_monly(years(y),vn)
        vn@lev = 850
        vn = "crain"
        cr = read_reanalysis_monly(years(y),vn)
        vn = "srain"
        sr = read_reanalysis_monly(years(y),vn)
        rain = rain + cr + sr
    end do
    u = u/ny
    v = v/ny
    rain = rain/ny
    u!1 = "lat"
    u!2 = "lon"
    v!1 = "lat"
    v!2 = "lon"
    rain!1 = "lat"
    rain!2 = "lon"
    ;; compset for mons
    if (nm.eq.1)then
        pu   = u(lat|:,lon|:,{time|mons})
        pv   = v(lat|:,lon|:,{time|mons})
        pr   = rain(lat|:,lon|:,{time|mons})
    else
        pu   = dim_avg_Wrap(u(lat|:,lon|:,{time|mons}))
        pv   = dim_avg_Wrap(v(lat|:,lon|:,{time|mons}))
        pr   = dim_avg_Wrap(rain(lat|:,lon|:,{time|mons}))
    end if
    if(isatt(years,"clm") .and. years@clm)then
        pu@clm = True
    else
        vn = "crain"
        clmcr = clmMonTLL(read_reanalysis_monly(allyears,vn))
        vn = "srain"
        clmsr = clmMonTLL(read_reanalysis_monly(allyears,vn))
        clmr = clmcr
        clmr = clmcr + clmsr
        pr = pr - dim_avg_Wrap(clmr(lat|:,lon|:,{month|mons-1}))
    end if
    pu@lev = lev
    a = plot_strm_rain(pu,pv,pr,filename,title)
    return True
end

undef("plot_strm_SONdJJA")
function plot_strm_SONdJJA(years[*],filename,title,opt)
begin
    load "res_years.ncl"
    ny = dimsizes(years)
    if(isatt(opt,"lev"))then
        lev = opt@lev
    else
        lev = 850
    end if

    vn = "U"
    vn@lev = lev
    u = read_reanalysis_monclm(years,vn)
    vn = "V"
    v = read_reanalysis_monclm(years,vn)
    vn@lev = lev
    vn = "crain"
    cr = read_reanalysis_monclm(years,vn)
    vn = "srain"
    sr = read_reanalysis_monclm(years,vn)
    rain = cr
    rain = cr + sr

    pu   = u(0,:,:)
    pv   = v(0,:,:)
    pr   = rain(0,:,:)
    pu   = (/dim_avg_n_Wrap(u({SON},:,:) - u({JJA},:,:),0)/)
    pv   = (/dim_avg_n_Wrap(v({SON},:,:) - v({JJA},:,:),0)/)
    pr   = (/dim_avg_n_Wrap(rain({SON},:,:) - rain({JJA},:,:),0)/)

    pu@lev = lev
    a = plot_strm_rain(pu,pv,pr,filename,title)
    return True
end

