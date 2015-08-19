load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
undef("plot_eof")
function plot_eof(race[*][*][*]) ; race(time,lat,lon)
begin
    runmean = 13
    neof = 10
    eofopt = True
    eofopt@pcrit = 100.
    yb = floattoint(min(race&time)/100)
    ye = floattoint(max(race&time)/100)
    ny = ye-yb+1

    eofrace   = eofunc_Wrap(race(lat|:,lon|:,time|:),neof,eofopt)
    eofracets = eofunc_ts_Wrap(race(lat|:,lon|:,time|:),eofrace,eofopt)
    ;; normalize 
    do i = 0, neof-1
        eofrace(i,:,:) = eofrace(i,:,:)/max(abs(eofrace(i,:,:)))
        eofracets(i,:) = eofracets(i,:)/max(abs(eofracets(i,:)))
        if (max(abs(eofracets(i,:))) .ne. max(eofracets(i,:)))then
            eofrace(i,:,:) = eofrace(i,:,:) * -1.
            eofracets(i,:) = eofracets(i,:) * -1.
        end if
    end do

    wks = gsn_open_wks("ps","eof")
    gsn_define_colormap(wks,"BlWhRe") 
    resEOF = True
        resEOF@gsnAddCyclic = False
        resEOF@gsnSpreadColors      = True         ; spread out color table
        resEOF@cnFillOn             = True         ; turn on color fill
        resEOF@cnLinesOn            = False    
        resEOF@cnFillMode           = "RasterFill"
        resEOF@cnLevelSelectionMode = "ManualLevels"
        resEOF@cnLevelSpacingF      =  0.1
        resEOF@cnMaxLevelValF       =  1.0
        resEOF@cnMinLevelValF       = -1.0
        resEOF@lbLabelStride        = 2
        resEOF@gsnStringFont  = 21
        resEOF@gsnStringFontHeightF     = 0.02
        resEOF@lbOrientation        = "Vertical"   ; vertical label bar
        resEOF@lbLeftMarginF        = -0.2
        resEOF@lbRightMarginF       =  0.3
        resEOF@gsnFrame = False
        resEOF@gsnDraw  = False

    rests= resEOF
    resEOF@mpCenterLonF         = 180.
    resEOF@mpMinLatF            = min(race&lat)         ; zoom in on map
    resEOF@mpMaxLatF            = max(race&lat)
    resEOF@mpMinLonF            = min(race&lon)
    resEOF@mpMaxLonF            = max(race&lon)

        rests@cnFillMode            = "AreaFill"
        rests@tmXBMode              = "Explicit"
        rests@tmXBValues            = ispan(0,11,1)
        rests@tmXBLabels            = ispan(1,12,1)
        rests@tmYLMode              = "Explicit"
        rests@tmYLValues            = ispan(0,ny-1,5)
        rests@tmYLLabels            = ispan(yb,ye,5)
        rests@tmYLLabelFont         = 12
        rests@lbLabelStride         = 5


    rts     = True
        rts@vpHeightF = 0.20        ; Changes the aspect ratio
        rts@vpWidthF  = 0.75
        rts@vpXF      = 0.10        ; change start locations
        rts@vpYF      = 0.85        ; the plot
        rts@gsnYRefLine           = 0.              ; reference line   
        rts@gsnXYBarChart         = True            ; create bar chart 
        rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
        rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
        rts@gsnStringFontHeightF = resEOF@gsnStringFontHeightF
        rts@gsnStringFont        = resEOF@gsnStringFont
        rts@gsnFrame = False
        rts@gsnDraw  = False
        rts@trYMinF               = -1.
        rts@trYMaxF               =  1.
        rts@xyLineColors            = (/-1,-1/)
        rts@tmXBMode = "Explicit"
        rts@tmXBValues = ispan(0,(ny*12)-1,60)
        rts@tmXBLabels = ispan(yb,ye,5)
        rts@tmYLLabelFontHeightF = 0.012
    plres = True
    plres@gsLineColor    = "black"
    plres@gsLineThicknessF = 5.0

    mons = fspan(0,(ny*12)-1,ny*12)

    delete(rts@tiMainString) ; = ""
    resEOF@gsnLeftString     = "(A)"
    rts@gsnLeftString     = "(B)"
    pres = True
    pres@txFont = 21
    do i = 0,neof-1
        pres@txString     = race@vartitle+" EOF "+(i+1)
        resEOF@gsnRightString   = sprintf("%4.2f",eofrace@pcvar(i))+"%"
        ;rts@gsnRightString      = eofrace@pcvar(i)+"%"
        rests@gsnRightString    = eofrace@pcvar(i)+"%"
        plot1 = gsn_csm_contour_map(wks,eofrace(evn|i,lat|:,lon|:),resEOF)
        plot2 = gsn_csm_xy (wks,mons,eofracets(i,:),rts)
        if (runmean .gt. 1)then
            ponly = where(eofracets(i,:).ge.0,eofracets(i,:),0.)
            rmean = runave_Wrap(ponly,6,0)
            yearly = rmean(7::12)
            rmyearly = runave_Wrap(yearly,runmean,0)
            plot0 = gsn_add_polyline(wks, plot2, mons(7::12), yearly, plres)
        end if
        pres@gsnPanelYWhiteSpacePercent = 0.
        gsn_panel(wks,(/plot1,plot2/),(/2,1/),pres)
        ;draw(plot2)
        ;frame(wks)
        ym   = onedtond(eofracets(i,:),(/ny,12/))
        plot = gsn_csm_contour(wks,ym,rests)
    end do 
    return True
end

undef("plot_cal_eof")
function plot_cal_eof(var[*][*][*],title,filename) ; race(time,lat,lon)
begin
    runmean = 0
    neof = 10
    eofopt = True
    eofopt@pcrit =  100.
    yb = floattoint(min(var&time)/100)
    ye = floattoint(max(var&time)/100)
    ny = ye-yb+1

    eofvar   = eofunc_Wrap(var(lat|:,lon|:,time|:),neof,eofopt)
    eofvarts = eofunc_ts_Wrap(var(lat|:,lon|:,time|:),eofvar,eofopt)
    ;; normalize 
    do i = 0, neof-1
        eofvar(i,:,:) = eofvar(i,:,:)/max(abs(eofvar(i,:,:)))
        eofvarts(i,:) = eofvarts(i,:)/max(abs(eofvarts(i,:)))
        if (max(abs(eofvarts(i,:))) .ne. max(eofvarts(i,:)))then
            eofvar(i,:,:) = eofvar(i,:,:) * -1.
            eofvarts(i,:) = eofvarts(i,:) * -1.
        end if
    end do

    print("plot "+filename)
    wks = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"BlWhRe") 
    resEOF = True
        resEOF@gsnAddCyclic = False
        resEOF@gsnSpreadColors      = True         ; spread out color table
        resEOF@cnFillOn             = True         ; turn on color fill
        resEOF@cnLinesOn            = False    
        resEOF@cnFillMode           = "RasterFill"
        resEOF@cnLevelSelectionMode = "ManualLevels"
        resEOF@cnLevelSpacingF      =  0.1
        resEOF@cnMaxLevelValF       =  1.0
        resEOF@cnMinLevelValF       = -1.0
        resEOF@lbLabelStride        = 2
        resEOF@gsnStringFont  = 21
        resEOF@gsnStringFontHeightF     = 0.02
        resEOF@lbOrientation        = "Vertical"   ; vertical label bar
        resEOF@lbLeftMarginF        = -0.2
        resEOF@lbRightMarginF       =  0.3
        resEOF@gsnFrame = False
        resEOF@gsnDraw  = False

    rests= resEOF
        resEOF@mpCenterLonF         = 180.
        resEOF@mpMinLatF            = min(var&lat)         ; zoom in on map
        resEOF@mpMaxLatF            = max(var&lat)
        resEOF@mpMinLonF            = min(var&lon)
        resEOF@mpMaxLonF            = max(var&lon)

        rests@cnFillMode            = "AreaFill"
        rests@tmXBMode              = "Explicit"
        rests@tmXBValues            = ispan(0,11,1)
        rests@tmXBLabels            = ispan(1,12,1)
        rests@tmYLMode              = "Explicit"
        rests@tmYLValues            = ispan(0,ny-1,5)
        rests@tmYLLabels            = ispan(yb,ye,5)
        rests@tmYLLabelFont         = 12
        rests@lbLabelStride         = 5


    rts     = True
        rts@vpHeightF = 0.20        ; Changes the aspect ratio
        rts@vpWidthF  = 0.75
        rts@vpXF      = 0.10        ; change start locations
        rts@vpYF      = 0.85        ; the plot
        rts@gsnYRefLine           = 0.              ; reference line   
        rts@gsnXYBarChart         = True            ; create bar chart 
        rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
        rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
        rts@gsnStringFontHeightF = resEOF@gsnStringFontHeightF
        rts@gsnStringFont        = resEOF@gsnStringFont
        rts@gsnFrame = False
        rts@gsnDraw  = False
        rts@trYMinF               = -1.
        rts@trYMaxF               =  1.
        rts@xyLineColors            = (/-1,-1/)
        rts@tmXBMode = "Explicit"
        rts@tmXBValues = ispan(0,(ny*12)-1,60)
        rts@tmXBLabels = ispan(yb,ye,5)
        rts@tmYLLabelFontHeightF = 0.012
    plres = True
        plres@gsLineColor    = "black"
        plres@gsLineThicknessF = 5.0

    mons = fspan(0,(ny*12)-1,ny*12)

    ;delete(rts@tiMainString) ; = ""
    rts@tiMainString = ""
    resEOF@gsnLeftString     = "(A)"+title
    rts@gsnLeftString     = "(B) Principal component"
    pres = True
    pres@txFont = 21
    eofvarts@long_name = ""
    do i = 0,neof-1
        ;;pres@txString     = var@vartitle+" EOF "+(i+1)
        ;;pres@txString     = title+" EOF "+(i+1)
        resEOF@gsnRightString   = sprintf("%4.2f",eofvar@pcvar(i))+"%"
        ;rts@gsnRightString      = eofvar@pcvar(i)+"%"
        rests@gsnRightString    = eofvar@pcvar(i)+"%"
        plot1 = gsn_csm_contour_map(wks,eofvar(evn|i,lat|:,lon|:),resEOF)
        plot2 = gsn_csm_xy (wks,mons,eofvarts(i,:),rts)
        if (runmean .gt. 1)then
            ponly = where(eofvarts(i,:).ge.0,eofvarts(i,:),0.)
            rmean = runave_Wrap(ponly,6,0)
            yearly = rmean(7::12)
            rmyearly = runave_Wrap(yearly,runmean,0)
            plot0 = gsn_add_polyline(wks, plot2, mons(7::12), yearly, plres)
        end if
        pres@gsnPanelYWhiteSpacePercent = 0.
        gsn_panel(wks,(/plot1,plot2/),(/2,1/),pres)
        ;draw(plot2)
        ;frame(wks)
        ym   = onedtond(eofvarts(i,:),(/ny,12/))
        plot = gsn_csm_contour(wks,ym,rests)
    end do 
    return True
end
