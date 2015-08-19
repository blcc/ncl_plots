load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

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
