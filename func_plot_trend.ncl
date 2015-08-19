undef("dtrendR2")
function dtrendR2(yrts[*][*][*]) ;; assume time,lat,lon
begin
    yrts@_FillValue = -99.99
    dims = dimsizes(yrts)
    dtrendyr = dtrend_msg_n(ispan(1,dims(0),1),yrts,True,True,0)
    yin  = new(dims(1:),typeof(dtrendyr))
    slope= new(dims(1:),typeof(dtrendyr)) 
    yin  = onedtond(dtrendyr@y_intercept,dims(1:))
    slope= onedtond(dtrendyr@slope,dims(1:))
    copy_VarCoords(yrts(0,:,:),yin)
    copy_VarCoords(yrts(0,:,:),slope)
    copy_VarCoords(yrts,dtrendyr)

    SST  = new(dims(1:),typeof(dtrendyr))
    SSE  = new(dims(1:),typeof(dtrendyr))
    SSR  = new(dims(1:),typeof(dtrendyr))

    SSE = dim_sum_n(dtrendyr*dtrendyr,0)
    ano = dim_rmvmean_n(yrts,0)
    SST = dim_sum_n(ano*ano,0)
    SST = where(SST.eq.0,SST@_FillValue,SST)
    SSR = SST-SSE
    R2  = SSR/SST

    dtrendyr@R2 = R2

    return dtrendyr
end
undef("do_trend_plot")
function do_trend_plot(dtrendyr[*][*][*],title,filename) ;; preprocess by dtrendR2()
begin
    dims = dimsizes(dtrendyr)
    yin  = new(dims(1:),typeof(dtrendyr))
    slope= new(dims(1:),typeof(dtrendyr)) 
    R2   = new(dims(1:),typeof(dtrendyr)) 
    yin  = onedtond(dtrendyr@y_intercept,dims(1:))
    slope= onedtond(dtrendyr@slope,dims(1:))
    R2   = onedtond(dtrendyr@R2,dims(1:))
    if(isatt(dtrendyr,"sscale"))then ;; slope scale
        slope = slope*dtrendyr@sscale
    end if
    copy_VarCoords(dtrendyr(0,:,:),yin)
    copy_VarCoords(dtrendyr(0,:,:),slope)
    copy_VarCoords(dtrendyr(0,:,:),R2)
    slope@_FillValue = -99.99   ;; default by dtrend_msg_n()
    slope = where(slope .lt. -90,slope@_FillValue,slope) ;; well...

    wks = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
    res            = True
        ;;res@tiMainString = title
        res@gsnLeftString = title
        res@gsnDraw    = False                         ; don't draw
        res@gsnFrame   = False                         ; don't advance frame
        res@vpWidthF   = 0.80                          ; make map bigger
        res@vpHeightF  = 0.80
        res@mpMaxLatF  =  40.                           ; select subregion
        res@mpMinLatF  =  00.
        res@mpMinLonF  = 110.
        res@mpMaxLonF  = 180.
        res@mpCenterLonF = 180.
        res@mpFillDrawOrder = "PreDraw"
        res@mpFillOn = False
        res@cnFillOn = True
        res@cnLinesOn             = False   
        ;;res@gsnSpreadColors  = True
        res@cnFillDrawOrder = "Draw"
        res@lbLabelAutoStride  = True
        res@cnFillColors = (/"white","yellow","green"/)
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels  = (/0.25,0.49/)

    res2 = True
        res2@cnFillOn = False
        res2@cnInfoLabelOn = False
        res2@gsnContourNegLineDashPattern = 1
        if(isatt(dtrendyr,"cmax"))then
            res2@cnLevelSelectionMode = "ManualLevels"
            res2@cnMinLevelValF = dtrendyr@cmin
            res2@cnMaxLevelValF = dtrendyr@cmax
            res2@cnLevelSpacingF = dtrendyr@cint ; manually sets the contour levels.
        end if


    print("trend plot: "+filename)
    print(max(abs(slope({0:30},{110:180}))))
    ;plot = gsn_csm_contour_map_ce(wks,slope,res)
    plot = gsn_csm_contour_map_overlay(wks,R2,slope,res,res2)

    ;; plot MDR square, bottom left corner start clockwise
    linex = (/110.,110.,160.,160.,110./)
    liney = (/ 10., 30., 30., 10., 10./)
    ;; plot MGR square, bottom left corner start clockwise
    linex = (/110.,110.,160.,160.,110./)
    liney = (/ 10., 25., 25., 10., 10./)
    mdrres = True
    mdrres@gsLineColor  ="red" 
    mdrres@gsLineThicknessF = 5.0    
    c = gsn_add_polyline(wks,plot,linex,liney,mdrres)

    draw(plot)
    frame(wks)

    return True
end
