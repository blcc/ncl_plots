load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
undef("plotT")
function plotT(dfile,ores)
begin
    print("ploting T...")
    res = ores
    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res@cnMinLevelValF       =  200.               ; set min contour level
    res@cnMaxLevelValF       =  300.               ; set max contour level
    res@cnLevelSpacingF      =     .5              ; set contour spacing
    res@lbBoxLinesOn         = False
    res@lbLabelStride        = 50

 
    var = dfile->T
    dims = dimsizes(var)

    wks = gsn_open_wks("ps","T")                  ; open a ps file
    gsn_define_colormap(wks,"3gauss")        ; choose colormap
    do i = 0,dims(0) -1
        varPlot = var(i,:,:)
        res@gsnCenterString  = sprintf("%4.0f",varPlot@lev)+" hPa"   ; plot center string
        plot = gsn_csm_contour_map_ce(wks,varPlot,res)  ; create a default plot
        draw(plot)
        frame(wks)
        delete(varPlot)
    end do
    return True
end
undef("plotTS")
function plotTS(dfile,ores)
begin
    print("ploting TS...")
    res = ores
    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res@cnMinLevelValF       =  260.               ; set min contour level
    res@cnMaxLevelValF       =  310.               ; set max contour level
    res@cnLevelSpacingF      =    .5               ; set contour spacing
    res@lbLabelStride        = 20
    var = dfile->TS(:,:)
    wks = gsn_open_wks("ps","TS")                  ; open a ps file
    gsn_define_colormap(wks,"3gauss")        ; choose colormap
    plot = gsn_csm_contour_map_ce(wks,var,res)  ; create a default plot
    draw(plot)
    frame(wks)
    return True
end
undef("plotUVT")
function plotUVT(dfile,ores)
begin
    print("ploting UVT...")
    res = ores
    res@gsnScalarContour     = True               ; contours desired

    res@vcGlyphStyle            = "CurlyVector"     ; turn on curley vectors
    res@vcRefMagnitudeF         = 3.0     ; define vector ref mag
    res@vcMinDistanceF  = 0.03               ; thin out vectors
    res@vcRefLengthF    = 0.045              ; define length of vec ref
    res@vcLineArrowColor        = "black"           ; change vector color
    res@vcLineArrowThicknessF   = 1.4               ; change vector thickness
    res@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
    res@mpFillOn                = False
    res@mpFillDrawOrder         = "Draw"
    res@lbLabelStride           = 50
    res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
    res@cnMinLevelValF       =  200.               ; set min contour level
    res@cnMaxLevelValF       =  300.               ; set max contour level
    res@cnLevelSpacingF      =    .5               ; set contour spacing

    wks = gsn_open_wks("ps","UVT")                  ; open a ps file
    gsn_define_colormap(wks,"3gauss")        ; choose colormap
    U   = dfile->U
    V   = dfile->V
    T   = dfile->T
    dims= dimsizes(T)
    do i = 0,dims(0) -1
        if (i.eq.26)then
        end if
        UPlot = U(i,:,:)
        VPlot = V(i,:,:)
        TPlot = T(i,:,:)
        res@gsnCenterString  = sprintf("%4.0f",TPlot@lev)+" hPa"   ; plot center string
        plot  = gsn_csm_vector_scalar_map_ce(wks,UPlot,VPlot,TPlot,res)
        draw(plot)
        frame(wks)
        delete(UPlot)
        delete(VPlot)
        delete(TPlot)
    end do
    
    return True
end

undef("plotZ")
function plotZ(dfile,ores)
begin
    print("ploting Z...")
    res = ores
    res@cnLevelSelectionMode = "AutomaticLevels"     ; set manual contour levels
    res@cnFillOn = False
    res@cnLinesOn = True
    res@cnLevelSpacingF      =    5.               ; set contour spacing
    var = dfile->Z3
    wks = gsn_open_wks("ps","Z")                  ; open a ps file
    dims = dimsizes(var)
    do i = 0,dims(0)-1
        varPlot = var(i,:,:)
        res@gsnCenterString  = sprintf("%4.0f",varPlot@lev)+" hPa"   ; plot center string
        plot = gsn_csm_contour_map_ce(wks,varPlot,res)  ; create a default plot
        delete(varPlot)
        draw(plot)
        frame(wks)
    end do
    return True
end
