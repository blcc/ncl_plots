load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "func_read_era_all.ncl"

undef("toTycData")
function toTycData(data,clat,clon)
begin
    plotRadius  = 25     ;; lat/lon degrees
    plotResolu  = 2.5    ;; lat/lon degrees
    plotGrids   = floattoint((plotRadius*2/plotResolu) +1)
    latB        = clat-plotRadius
    latE        = clat+plotRadius
    lonB        = clon-plotRadius
    lonE        = clon+plotRadius
    plotLatList = fspan(latB,latE,plotGrids)
    plotLonList = fspan(lonB,lonE,plotGrids)
    plotLatList@units = "degree_north"
    plotLonList@units = "degree_east"
    plotLats    = new((/plotGrids,plotGrids/),"float")
    plotLons    = new((/plotGrids,plotGrids/),"float")
    ;; rgrid2rcm
    do i = 0, plotGrids-1
        plotLats(i,:) = plotLatList
        plotLons(:,i) = plotLonList
    end do
    odata    = rgrid2rcm(data&lat,data&lon, data, \
                            plotLats    , plotLons   , 0)
    return odata
end

undef("plot_tycplot")
function  plot_tycplot(u[*][*],v[*][*],z[*][*],lev,clat,clon,filename,title)
begin
    ;; 以颱風為中心作風場圖
    plotRadius  = 25     ;; lat/lon degrees
    plotResolu  = 2.5    ;; lat/lon degrees
    plotGrids   = floattoint((plotRadius*2/plotResolu) +1)
    latB        = clat-plotRadius
    latE        = clat+plotRadius
    lonB        = clon-plotRadius
    lonE        = clon+plotRadius

    resTyC                  = True
    resTyC@gsnDraw          = False
    resTyC@gsnFrame         = False
    resTyC@gsnAddCyclic     = False            ; regional data
    resTyC@vcRefLengthF     =  .03
    resTyC@tiMainString     = title
    resTyC@gsnLeftString    = ""
    resTyC@gsnRightString   = ""
    resTyC@vcGlyphStyle = "WindBarb"
        u = u/0.514   ;; m/s to knots
        v = v/0.514 

    pres = True
    pres@gsMarkerColor  = "magenta"
    pres@gsMarkerIndex = 16
    pres@gsMarkerSizeF = 0.010
    pres@gsFillDotSizeF = 0.5
    pres@gsLineColor   = "red"
    pres@gsLineThicknessF = 3.0

    if(lev.eq.200)then
        resTyC@vcRefMagnitudeF = 50.
        resTyC@cnLevelSpacingF = 30.
    end if
    if(lev.eq.500)then
        resTyC@vcRefMagnitudeF = 30.
    end if
    if(lev.eq.850)then
        resTyC@vcRefMagnitudeF = 15.
    end if
    if(lev.eq.700)then
        resTyC@vcRefMagnitudeF = 20.
    end if

    wks2 = gsn_open_wks("ps",filename)
    if(clat.eq.0.and.clon.eq.0.)then
        plot     = gsn_csm_vector(wks2,u,v,resTyC)
        plotp    = gsn_add_polymarker(wks2,plot,0.,0.,pres)
    else
        resTyC@mpMaxLatF        =  latE
        resTyC@mpMinLatF        =  latB
        resTyC@mpMaxLonF        =  lonE
        resTyC@mpMinLonF        =  lonB
        resTyC@mpFillOn         = True
        resTyC@mpLandFillColor  = "gray"           ; set land to be
        plot     = gsn_csm_vector_map_ce(wks2,u,v,resTyC)
        plotp    = gsn_add_polymarker(wks2,plot,clon,clat,pres)
    end if

    draw(plot)
    frame(wks2)
    delete(plot)
    delete(wks2)

    return True
end


undef("plot_thickness_addpoint")
function plot_thickness_addpoint(dz[*][*],xpts[*],ypts[*],filename,ores)
begin
    res = ores
    load "res_map.ncl"
    res@mpFillOn = True
    res@mpMinLonF =   90
    res@mpMaxLonF =  210
    res@gsnLeftString = ""
    res@gsnRightString = ""

    pres = True
    pres@gsMarkerIndex = 16          ; Use filled dots for markers.
    pres@gsMarkerThicknessF = 3.
    pres@tfPolyDrawOrder = "PostDraw"
    print("Regenv: "+filename)

    wks  = gsn_open_wks("ps",filename)
    plot = gsn_csm_contour_map_ce(wks,dz,res)

    a = gsn_add_polymarker(wks,plot,xpts,ypts,pres)
    draw(plot)
    frame(wks)
    return True
end


undef("alignarray")
function alignarray(dd[*])
begin
    n = num(.not.ismissing(dd))
    od = new(n,typeof(dd))
    j = 0
    do i = 0,n-1
        if(.not.ismissing(dd(i)))then
            od(j) = dd(i)
        end if
    end do
    return od
end

undef("plot_tycenv")
function plot_tycenv(iyear[*]:integer,imon[*]:integer,iday[*]:integer,ihr[*]:integer,iclat[*],iclon[*],lev,res)
begin
    year = alignarray(iyear)
    mon  = alignarray(imon)
    day  = alignarray(iday)
    hr   = alignarray(ihr)
    clat = alignarray(iclat)
    clon = alignarray(iclon)
    nt = dimsizes(iyear)
    if(isatt(res,"filename"))then
        filename = res@filename
        else
        filename = "tyc_env_"+year(0)+"-"+sprinti("%2.2i",mon(0))+"-"+sprinti("%2.2i",day(0))+"-"+sprinti("%2.2i",hr(0))
        end if
    if(isatt(res,"title"))then
        title =  res@title
        else
        title = filename
        end if

    vn = "U"
    vn@lev = lev
    if(nt.eq.1)then
        u = read_ERA6hr(year,mon,day,hr,vn)
        vn = "V"
        v = read_ERA6hr(year,mon,day,hr,vn)
        vn = "Z"
        z = read_ERA6hr(year,mon,day,hr,vn)
        a = plot_tycplot(u,v,z,lev,clat,clon,filename,title)
        return True
    end if
    
    u = toTycData(read_ERA6hr(year(0),mon(0),day(0),hr(0),vn),clat(0),clon(0))
    vn = "V"
    v = toTycData(read_ERA6hr(year(0),mon(0),day(0),hr(0),vn),clat(0),clon(0))
    vn = "Z"
    z = toTycData(read_ERA6hr(year(0),mon(0),day(0),hr(0),vn),clat(0),clon(0))

    do i = 1, nt -1
        vn = "U"
        u  = u + toTycData(read_ERA6hr(year(i),mon(i),day(i),hr(i),vn),clat(i),clon(i))
        vn = "V"
        v  = v + toTycData(read_ERA6hr(year(i),mon(i),day(i),hr(i),vn),clat(i),clon(i))
        vn = "Z"
        z  = z + toTycData(read_ERA6hr(year(i),mon(i),day(i),hr(i),vn),clat(i),clon(i))
    end do
    u = u/nt
    v = v/nt
    z = z/nt
    a = plot_tycplot(u,v,z,lev,0.,0.,filename,title)
    delete(u)
    delete(v)
    delete(z)
    return True
end

undef("plot_tyenv")
function plot_tyenv(iyear[*]:integer,imon[*]:integer,iday[*]:integer,ihr[*]:integer,iclat[*],iclon[*],lev,res)
begin
    year = alignarray(iyear)
    mon  = alignarray(imon)
    day  = alignarray(iday)
    hr   = alignarray(ihr)
    clat = alignarray(iclat)
    clon = alignarray(iclon)
    nt = dimsizes(iyear)
    if(isatt(res,"filename"))then
        filename = res@filename
        else
        filename = "env_"+year(0)+"-"+sprinti("%2.2i",mon(0))+"-"+sprinti("%2.2i",day(0))+"-"+sprinti("%2.2i",hr(0))
        end if
    if(isatt(res,"title"))then
        title =  res@title
        else
        title = filename
        end if

    vn = "U"
    vn@lev = lev
    if(nt.eq.1)then
        u = read_ERA6hr(year,mon,day,hr,vn)
        vn = "V"
        v = read_ERA6hr(year,mon,day,hr,vn)
        vn = "Z"
        z = read_ERA6hr(year,mon,day,hr,vn)
        a = plot_tycplot(u,v,z,lev,clat,clon,filename,title)
        return True
    end if
    
    u = toTycData(read_ERA6hr(year(0),mon(0),day(0),hr(0),vn),clat(0),clon(0))
    vn = "V"
    v = toTycData(read_ERA6hr(year(0),mon(0),day(0),hr(0),vn),clat(0),clon(0))
    vn = "Z"
    z = toTycData(read_ERA6hr(year(0),mon(0),day(0),hr(0),vn),clat(0),clon(0))

    do i = 1, nt -1
        vn = "U"
        u  = u + toTycData(read_ERA6hr(year(i),mon(i),day(i),hr(i),vn),clat(i),clon(i))
        vn = "V"
        v  = v + toTycData(read_ERA6hr(year(i),mon(i),day(i),hr(i),vn),clat(i),clon(i))
        vn = "Z"
        z  = z + toTycData(read_ERA6hr(year(i),mon(i),day(i),hr(i),vn),clat(i),clon(i))
    end do
    u = u/nt
    v = v/nt
    z = z/nt
    a = plot_tycplot(u,v,z,lev,0.,0.,filename,title)
    delete(u)
    delete(v)
    delete(z)
    return True
end


undef("plot_tyenv_thickness")
function plot_tyenv_thickness(iyear[*]:integer,imon[*]:integer,iday[*]:integer,ihr[*]:integer,iclat[*],iclon[*],lev1,lev2,res)
begin
    year = alignarray(iyear)
    mon  = alignarray(imon)
    day  = alignarray(iday)
    hr   = alignarray(ihr)
    clat = alignarray(iclat)
    clon = alignarray(iclon)
    nt = dimsizes(iyear)
    if(isatt(res,"filename"))then
        filename = res@filename
        else
        filename = "env_dz_"+year(0)+"-"+sprinti("%2.2i",mon(0))+"-"+sprinti("%2.2i",day(0))+"-"+sprinti("%2.2i",hr(0))
        end if
    if(isatt(res,"title"))then
        title =  res@title
        else
        title = filename
        end if

    vn = "Z"
    if(nt.eq.1)then
        vn@lev = min((/lev1,lev2/))
        zupper = read_ERA6hr(year,mon,day,hr,vn)
        vn@lev = max((/lev1,lev2/))
        zlower = read_ERA6hr(year,mon,day,hr,vn)
        dz = zlower
        dz = zlower - zupper
        a = plot_thickness_addpoint(dz,clon,clat,filename,res)
        return True
    end if
    
    do i = 1, nt -1
        vn@lev = min((/lev1,lev2/))
        zupper  = zupper + read_ERA6hr(year(i),mon(i),day(i),hr(i),vn)
        vn@lev = max((/lev1,lev2/))
        zlower  = zlower + read_ERA6hr(year(i),mon(i),day(i),hr(i),vn)
    end do
    zupper = zupper/nt
    zlower = zlower/nt
    dz = zlower
    dz = zlower - zupper
    a = plot_thickness_addpoint(dz,clon,clat,filename,res)
    return True
end


