load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("plot_latlon")
function plot_latlon(title,filename,var[*][*],ires)
begin
    ; input figtitle, fig filename and 2d var
    if(all(ismissing(var)))then
        print("plot plot_latlon(): input var is all missing.")
        return False
    end if
    print("plot plot_latlon(): "+filename)

    res = True
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnSpreadColors = True
    res@mpFillOn              =  False        
    res@mpCenterLonF  = 180.
    res@tiMainString = ""
    res@gsnStringFont = 22
    res@gsnLeftString = title
    if(isatt(var,"units"))then
        res@gsnRightString = var@units
    end if
    res@lbLabelAutoStride = True
    
    if(min(var)*max(var).ge.0)then
        
    else
        symMinMaxPlt(var,21,True,res)
    end if

    wks = gsn_open_wks("ps",filename)
    colormap = "BlWhRe"
    copy_VarAtts(ires,res)
    if(isatt(res,"colormap"))then
        cc = stringtochar(res@colormap)
        lcc = dimsizes(cc)
        if(cc(lcc-5:lcc-1).eq.".rgb")then
            if(fileexists(res@colormap))then
                colormap := RGBtoCmap(res@colormap)
            end if
        else
            colormap = res@colormap
        end if
        delete_VarAtts(res, "colormap")
    end if
        ;;print("plot_latlon(): set color map to "+res@colormap)
    gsn_define_colormap(wks, colormap)

    if(isatt(res,"reversecolor"))then
        delete_VarAtts(res, "reversecolor")
        gsn_reverse_colormap(wks)
    end if
    if(isatt(res,"plotaddcng"))then  ;; will plot green contour
        addcng = res@plotaddcng   ;; should be (lat,lon)
        copy_VarCoords(var,addcng)
        delete_VarAtts(res,"plotaddcng")
        res@gsnFrame = False
        res@gsnDraw  = False
    end if
    if(isatt(res,"plotaddcnr"))then  ;; will plot red contour
        addcnr = res@plotaddcnr   ;; should be (lat,lon)
        copy_VarCoords(var,addcnr)
        delete_VarAtts(res,"plotaddcnr")
        res@gsnFrame = False
        res@gsnDraw  = False
    end if
    if(isatt(res,"plotaddvcu"))then  ;; will plot vector
        addvcu = res@plotaddvcu   ;; should be (lat,lon)
        addvcv = res@plotaddvcv   ;; should be (lat,lon)
        assignFillValue(var,addvcu)
        assignFillValue(var,addvcv)
        if(isatt(res,"plotaddvcc"))then   ;; vector color
            ;;; vector color is limited by colormap in wks
            ;;; still need improve to usable
            addvcc = res@plotaddvcc   ;; should be (lat,lon)
            copy_VarCoords(var,addvcc)
            delete_VarAtts(res,"plotaddvcc")
        end if
        ;;;not done
        copy_VarCoords(var,addvcu)
        copy_VarCoords(var,addvcv)
        delete_VarAtts(res,"plotaddvcu")
        delete_VarAtts(res,"plotaddvcv")
        res@gsnFrame = False
        res@gsnDraw  = False
    end if
    if(isatt(res,"plotline"))then
        lines = res@plotline   ;; should be (nline,npoint,lat/lon)
        delete_VarAtts(res,"plotline")
        res@gsnFrame = False
        res@gsnDraw  = False
        if(isatt(res,"plotlinecolors"))then
            linecolors = res@plotlinecolors
            delete_VarAtts(res,"plotlinecolors")
        else
            linecolors = "red"
        end if
        nlcolor = dimsizes(linecolor)
    end if

    if(isatt(res,"nomap") .and. res@nomap)then
        delete(res@nomap)
        plot = gsn_csm_contour(wks,var,res)
    else
        plot = gsn_csm_contour_map_ce(wks,var,res)
    end if
    if(isvar("lines"))then
        ldims = dimsizes(lines)
        ;;print("ploting lines: "+ldims(0))
        lineres = True
        lineres@gsLineColor = "mediumaquamarine"
        lineres@gsLineThicknessF  = 3
        do l = 0,ldims(0)-1
            lineres@gsLineColor = linecolors(where(l.gt.nlcolor,nlcolor,l))
            plotstr = "plot"+l
            ldims@$plotstr$ = gsn_add_polyline(wks,plot,lines(l,:,1),lines(l,:,0),lineres)
            plotstr = "plotm"+l
            maxp = min(ind(ismissing(lines(l,:,1))))-1
            lineres@gsMarkerIndex  = 1
            if(ismissing(maxp))then
                maxp = 0
                ;print(lines(l,:,0)+","+lines(l,:,1))
                ;exit
            end if
            ldims@$plotstr$ = gsn_add_polymarker(wks,plot,lines(l,maxp,1),lines(l,maxp,0),lineres)
        end do
    end if
    if(isvar("addcng"))then  ;; plot red contour
        gres = True
        gres@cnLineColor = "green"
        gres@cnLineLabelsOn = False
        gres@cnLineThicknessF = 2.
        gres@gsnRightString   = ""   
        gres@gsnLeftString    = ""
        gres@gsnCenterString  = ""
        gres@gsnDraw = False
        gres@gsnFrame = False
        gres@cnInfoLabelOn  = False
        plotr = gsn_csm_contour(wks,addcng,gres)
        overlay(plot,plotr)
    end if
    if(isvar("addcnr"))then  ;; plot red contour
        rres = True
        rres@cnLineLabelsOn = False
        rres@cnLineColor = "red"
        rres@cnLineThicknessF = 2.
        rres@gsnRightString   = ""   
        rres@gsnLeftString    = ""
        rres@gsnCenterString  = ""
        rres@gsnDraw = False
        rres@gsnFrame = False
        rres@cnInfoLabelOn  = False
        plotr = gsn_csm_contour(wks,addcnr,rres)
        overlay(plot,plotr)
    end if
    if(isvar("addvcu").and.(.not.all(ismissing(addvcu))))then  ;; plot vector
        ;;print("ploting vector")
        vres = True
            vres@vcRefAnnoOn = False
            vres@vcGlyphStyle            = "CurlyVector" 
            vres@vcMinDistanceF          = 0.017  
            vres@gsnDraw = False
            vres@gsnFrame = False
            vres@vcMinMagnitudeF          = 0.001
            vres@vcRefLengthF             = 0.045
            vres@vcRefMagnitudeF          = 20.0
            vres@vcLineArrowColor = "green"
        if(isvar("addvcc"))then
            plotv = gsn_csm_vector_scalar(wks,addvcu,addvcv,addvcc,vres)
        else
            plotv = gsn_csm_vector(wks,addvcu,addvcv,vres)
        end if
        ;;print("ploting vector before overlay")
        overlay(plot,plotv)
    end if
    if(isatt(res,"gsnDraw").and. .not.res@gsnDraw)then
        draw(plot)
    end if
    if(isatt(res,"gsnFrame").and. .not.res@gsnFrame)then
        frame(wks)
    end if

    return plot
end

undef("plot_cn_PL")
function plot_cn_PL(title,filename,var[*][*],ires)
begin
    ; assume var(p,latorlon)
    print("plot plot_cn_PL(): "+filename)
    if(isatt(var,"_FillValue_original"))then
        var@_FillValue = var@_FillValue_original
    end if

    res = True
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@tiMainString = ""
    res@gsnStringFont = 22
    res@gsnLeftString = title
    if(isatt(var,"units"))then
        res@gsnRightString = var@units
    end if
    res@lbLabelAutoStride = True
    
    if(min(var)*max(var).ge.0)then
        mmi = nice_mnmxintvl(min(var),max(var),21,True)
        res@cnLevelSelectionMode = "ManualLevels"
        res@cnMinLevelValF = mmi(0)
        res@cnMaxLevelValF = mmi(1)
        res@cnLevelSpacingF = mmi(2)
    else
        symMinMaxPlt(var,21,True,res)
    end if

    wks = gsn_open_wks("ps",filename)
    colormap = "BlWhRe"
    copy_VarAtts(ires,res)

    if(isatt(res,"colormap"))then
        colormap = res@colormap
        delete_VarAtts(res, "colormap")
        if(isatt(res,"clmcolormap"))then
            delete_VarAtts(res, "clmcolormap")
        end if
        if(isatt(res,"difcolormap"))then
            delete_VarAtts(res, "difcolormap")
        end if
    end if
    gsn_define_colormap(wks, colormap)

    if(isatt(res,"reversecolor"))then
        delete_VarAtts(res, "reversecolor")
        gsn_reverse_colormap(wks)
    end if

    if(isatt(res,"plotline"))then
        lines = res@plotline   ;; should be (nline,npoint,plev/x)
        delete_VarAtts(res,"plotline")
        res@gsnFrame = False
        res@gsnDraw  = False
        if(isatt(res,"plotlinecolors"))then
            linecolors = res@plotlinecolors
            delete_VarAtts(res,"plotlinecolors")
        else
            linecolors = "red"
        end if
        nlcolor = dimsizes(linecolors)
    end if

    plot = gsn_csm_pres_hgt(wks,var,res)

    if(isvar("lines"))then
        ldims = dimsizes(lines)
        ;;print("ploting lines: "+ldims(0))
        lineres = True
        ;;lineres@gsLineColor = "mediumaquamarine"
        lineres@gsLineThicknessF  = 3
        do l = 0,ldims(0)-1
            if(all(ismissing(lines(l,:,:))))then
                continue
            end if
            lineres@gsLineColor = linecolors(where(l.gt.nlcolor,nlcolor,l))
            plotstr = "plot"+l
            ldims@$plotstr$ = gsn_add_polyline(wks,plot,lines(l,:,1),lines(l,:,0),lineres)
            plotstr = "plotm"+l
            maxp = min(ind(ismissing(lines(l,:,1))))-1
            lineres@gsMarkerIndex  = 1
            if(ismissing(maxp))then
                maxp = 0
                ;print(lines(l,:,0)+","+lines(l,:,1))
                ;exit
            end if
            ldims@$plotstr$ = gsn_add_polymarker(wks,plot,lines(l,maxp,1),lines(l,maxp,0),lineres)
        end do
    end if

    if(isatt(res,"gsnDraw").and. .not.res@gsnDraw)then
        draw(plot)
    end if
    if(isatt(res,"gsnFrame").and. .not.res@gsnFrame)then
        frame(wks)
    end if

    return plot
end
undef("plot_vc_cn_PL")
function plot_vc_cn_PL(title,filename,var[*][*],xvar[*][*],zvar[*][*],ires)
begin
    ; assume var(p,latorlon)
    print("plot "+filename)
    if(isatt(var,"_FillValue_original"))then
        var@_FillValue = var@_FillValue_original
    end if

    res = True
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnSpreadColors = True
    res@vcRefLengthF = .03
    res@vcRefMagnitudeF = 30.
    res@vcGlyphStyle  = "CurlyVector"
    res@tiMainString = ""
    res@gsnStringFont = 22
    res@gsnLeftString = title
    res@vcMinDistanceF = 0.01
    vcxmaxpoint = 36
    dims = dimsizes(xvar)
    xstride = max((/1,(dims(1)/vcxmaxpoint)/))
    
    scale = 1.
    if(isatt(var,"units"))then
        res@gsnRightString = var@units+"*"+scale
    end if
    res@lbLabelAutoStride = True
    
    if(min(var)*max(var).ge.0)then
        mmi = nice_mnmxintvl(min(var),max(var),21,True)
        res@cnLevelSelectionMode = "ManualLevels"
        res@cnMinLevelValF = mmi(0)
        res@cnMaxLevelValF = mmi(1)
        res@cnLevelSpacingF = mmi(2)
    else
        symMinMaxPlt(var,21,True,res)
    end if

    wks = gsn_open_wks("ps",filename)
    colormap = "BlWhRe"
    copy_VarAtts(ires,res)

    if(isatt(res,"colormap"))then
        colormap = res@colormap
        delete_VarAtts(res, "colormap")
    end if
    gsn_define_colormap(wks, colormap)

    if(isatt(res,"reversecolor"))then
        delete_VarAtts(res, "reversecolor")
        gsn_reverse_colormap(wks)
    end if

    pzvar = zvar
    pzvar = zvar*scale
    pxvar = xvar
    plot = gsn_csm_pres_hgt_vector(wks,var,pxvar(:,::xstride),pzvar(:,::xstride),res)

    return plot
end
undef("plot_Pprofiles2")
function plot_Pprofiles2(title,filename,varname1,var1[*],varname2,var2[*],ires)
begin
    ; assume 2 vars
    res = ires
    print("plot plot_Pprofiles2(): "+filename)

    linecolors = (/"Black","Red"/)
    linedashpattern = (/0,1/)

    if(isatt(var1,"_FillValue_original"))then
        var1@_FillValue = var1@_FillValue_original
    end if
    if(isatt(var2,"_FillValue_original"))then
        var2@_FillValue = var2@_FillValue_original
    end if

    res = True
        res@tiMainString = ""
        res@gsnStringFont = 22
        res@gsnLeftString = title
        res@gsnFrame = False
        res@trYReverse        = True
        
        ;res@pmLegendDisplayMode    = "Always"
        ;res@pmLegendParallelPosF   = .90
        ;res@pmLegendOrthogonalPosF = -0.8 
        ;res@pmLegendWidthF         = 0.52
        ;res@pmLegendHeightF        = 0.25
    
        res@xyYStyle          = "Log"
        res@tmYLMode          = "Explicit"

        res@xyLineColor       = linecolors(0)
        res@xyDashPattern    = linedashpattern(0)
    
        res@tmYLValues        = (/ 1000, 700,500,400,300,200,100,50,30,10,5,1/)
        res@tmYLLabels        = ""+res@tmYLValues
    res1Leg = True
        res1Leg@gsLineColor = linecolors(0)
        res1Leg@gsLineDashPattern = linedashpattern(0)
    res2 = True
        res2@gsLineColor       = linecolors(1)
        res2@gsLineDashPattern    =  linedashpattern(1)

    resText = True
        resText@txFontHeightF        = 0.020

    wks = gsn_open_wks("ps",filename)
    colormap = "BlWhRe"
    copy_VarAtts(ires,res)

    if(isatt(res,"colormap"))then
        colormap = res@colormap
        delete_VarAtts(res, "colormap")
    end if
    gsn_define_colormap(wks, colormap)

    if(isatt(res,"reversecolor"))then
        delete_VarAtts(res, "reversecolor")
        gsn_reverse_colormap(wks)
    end if

    plot = gsn_csm_xy(wks,var1,var1&$var1!0$,res)
    gsn_polyline(wks,plot,var2,var2&$var2!0$,res2)

    gsn_text_ndc(wks,varname1,0.7,0.75 ,resText)
    gsn_text_ndc(wks,varname2,0.7,0.70,resText)
    gsn_polyline_ndc(wks,(/0.60,0.65/),(/0.75,0.75/),res1Leg)
    gsn_polyline_ndc(wks,(/0.60,0.65/),(/0.70,0.70/),res2)

    frame(wks)

    return plot
end
undef("plot_Pprofiles2_dif")
function plot_Pprofiles2_dif(title,filename,var1[*],var2[*],ires)
begin
    ; assume 2 vars
    res = ires
    print("plot "+filename)



    linecolors = (/"Black"/)
    linedashpattern = (/0/)

    if(isatt(var1,"_FillValue_original"))then
        var1@_FillValue = var1@_FillValue_original
    end if
    if(isatt(var2,"_FillValue_original"))then
        var2@_FillValue = var2@_FillValue_original
    end if

    nvar1 = int2p_Wrap(var1&$var1!0$,var1,var2&$var2!0$,2)
    vardif = nvar1
    vardif = nvar1 - var2

    res = True
        res@tiMainString = ""
        res@gsnStringFont = 22
        res@gsnLeftString = title
        res@gsnFrame = False
        res@trYReverse        = True
        
        res@xyYStyle          = "Log"
        res@tmYLMode          = "Explicit"

        res@xyLineColor       = linecolors(0)
        res@xyDashPattern    = linedashpattern(0)
    
        res@tmYLValues        = (/ 1000, 700,500,400,300,200,100,50,30,10,5,1/)
        res@tmYLLabels        = ""+res@tmYLValues

    resref = True
        resref@gsLineColor = "green"

    wks = gsn_open_wks("ps",filename)
    colormap = "BlAqGrYeOrReVi200"
    copy_VarAtts(ires,res)

    if(isatt(res,"colormap"))then
        colormap = res@colormap
        delete_VarAtts(res, "colormap")
    end if
    gsn_define_colormap(wks, colormap)

    if(isatt(res,"reversecolor"))then
        delete_VarAtts(res, "reversecolor")
        gsn_reverse_colormap(wks)
    end if

    plot = gsn_csm_xy(wks,vardif,vardif&$vardif!0$,res)
    gsn_polyline(wks,plot,(/0,0/),(/1000,100/),resref)

    frame(wks)

    return plot
end

undef("plot_cn")
function plot_cn(title,filename,var[*][*],ires)
begin
    ;; detet var dimname to plot lon-lat or lon-p or lat-p 
    dimnames = getvardims(var)
    if(isatt(var,"valid_range"))then
        rag = var@valid_range
        var = where((var.lt.rag(0)).or.(var.gt.rag(1)),var@_FillValue,var)
    end if
    lat = True ;str_index_of_substr(dimnames(0),"lat",-1).gt.-1
    lon = True ;str_index_of_substr(dimnames(1),"lon",-1).gt.-1
    if(lat.and.lon)then
        a = plot_latlon(title,filename,var,ires)
    else
        a = plot_cn_PL(title,filename,var,ires)
    end if
    return a
end

undef("plot_vc_cn")
function plot_vc_cn(title,filename,u[*][*],v[*][*],cn[*][*],ires)
begin
    print("plot vector-countour : "+filename)
    res = True
    res@gsnSpreadColors = True
    res@cnFillOn = True
    res@cnLinesOn = False
    res@tiMainString = ""
    res@gsnStringFont = 22
    res@gsnLeftString = title
    if(isatt(cn,"units"))then
        res@gsnRightString = cn@units
    end if
    res@lbLabelAutoStride = True
    res@pmLabelBarOrthogonalPosF = -0.15          ; move label bar closer
    res@mpFillOn              =  False        
    res@mpCenterLonF  = 180.

    res@gsnScalarContour     = True 
    res@vcLineArrowColor        = "darkgreen"           ; change vector color
    res@vcGlyphStyle            = "CurlyVector"     ; turn on curley vectors
    res@vcMinDistanceF = 0.02
    vref = nice_mnmxintvl (0.,max(abs(u)),15,True)
    res@vcRefMagnitudeF           = vref(1)             ; define vector ref mag
    res@vcRefLengthF              = 0.045           ; define length of vec ref
    res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color

    wks = gsn_open_wks("ps",filename)

    copy_VarAtts(ires,res)

    colormap = "BlWhRe"
    if(isatt(res,"colormap"))then
        colormap = res@colormap
        delete_VarAtts(res, "colormap")
    end if

    if(isatt(res,"reversecolor"))then
        delete_VarAtts(res, "reversecolor")
        gsn_reverse_colormap(wks)
    end if

    gsn_define_colormap(wks, colormap)


    plot = gsn_csm_vector_scalar_map_ce(wks,u,v,cn,res)

    return plot
end

undef("diff_var")
function diff_var(var1,var2)
begin
    dims1 = dimsizes(var1)
    dims2 = dimsizes(var2)
    
    nvar1 = linint2(var1&$var1!1$,var1&$var1!0$,var1,False,var2&$var2!1$,var2&$var2!0$,0)
    vardif = nvar1 - var2
    copy_VarCoords(var2,vardif)
    return vardif
end

undef("diff_var_PL")
function diff_var_PL(var1,var2)
begin
    dims1 = dimsizes(var1)
    dims2 = dimsizes(var2)
    
    npvar1 = int2p_n_Wrap(var1&$var1!0$,var1,var2&$var2!0$,2,0)
    nvar1 = linint1_n_Wrap(npvar1&$npvar1!1$,npvar1,False,var2&$var2!1$,0,1)

    vardif = nvar1 - var2
    copy_VarCoords(var2,vardif)
    return vardif
end

undef("plot_diff_cn")
function plot_diff_cn(title[1]:string,filename[1]:string,var1,var2,res[1]:logical)
begin
    ;; plot contour of va1-var2
    vardif = diff_var(var1,var2)
    a = plot_cn(title,filename,vardif,res)
    return vardif
end

undef("plot_diff_cn_PL")
function plot_diff_cn_PL(title[1]:string,filename[1]:string,var1,var2,res[1]:logical)
begin
    ;; plot contour of va1-var2
    vardif = diff_var_PL(var1,var2)
    a = plot_cn_PL(title,filename,vardif,res)
    ;print(num(ismissing(vardif)))
    ;print(num(.not.ismissing(vardif)))
    return vardif
end

undef("plot_diff_vc_cn")
function plot_diff_vc_cn(title[1]:string,filename[1]:string,u1,v1,c1,u2,v2,c2,res[1]:logical)
begin
    ;; plot vector and contour of va1-var2
    udif = diff_var(u1,u2)
    vdif = diff_var(v1,v2)
    cdif = diff_var(c1,c2)
    a = plot_vc_cn(title,filename,udif,vdif,cdif,res)
    return a
end
;; test
    ;;load "$HOME/HiRAM_ncl/func_read_c360.ncl"
    
    ;a = read_1month("ta",1999,8)
    ;printVarSummary(a)
    ;res = True
    ;b = plot_cn("lon-lat map","test1",a(0,:,:),res)
    ;c = plot_cn("lat-p profile","test1",a(:,:,0),res)
    ;u = read_1month("c384","ucomp",1979,1)
    ;v = read_1month("c384","vcomp",1979,1)
    ;t = read_1month("c384","temp" ,1979,1)
    ;a = plot_vc_cn("test","test",u(5,::20,::20),v(5,::20,::20),t(5,::20,::20),True)


