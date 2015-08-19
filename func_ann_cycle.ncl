load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  ;;cmap = RGBtoCmap (rgb_file)
  ;; gsn_define_colormap(wks,cmap)
function plot_annual_cycle(title1:string,title2:string,filename:string,A[*][*],B[*][*]) 
begin
    plot_type = "ps"
    color_type = "COLOR" ;"MONO" or "COLOR"
    time_stamp = "True"
    version = "blc"
    lat1  = A&lat
    lat2  = B&lat
    nlat1 = dimsizes(A(0,:))
    nlat2 = dimsizes(B(0,:))
    
;----------------------- Open files for plots ----------------------
    wks = gsn_open_wks(plot_type,filename)
    ;gsn_define_colormap(wks,"BlWhRe")
    cmap = RGBtoCmap ("rgb/amwg.rgb")
    gsn_define_colormap(wks,cmap)
;------------------------------------------------------------------------

; case contour plots of time(x) vs. latitude(y)
    min1 = min(A)
    max1 = max(A)
    min2 = min(B)
    max2 = max(B)
    res = True

     plot = new(3,graphic)
    res                        = True
    res@gsnDraw                = False
    res@gsnFrame               = False
    res@txFontHeightF          = 0.025
    res@sfXArray               = ispan(0,12,1)
    res@tiMainFontHeightF      = 0.03
    res@tmXBMode               = "Explicit"
    res@tmXBValues             = ispan(0,12,1)
    res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S",\
                    "O","N","D","J"/)
    res@cnFillOn  = True
    res@cnLinesOn = False
    res@lbTitleOn = True
    res@lbLabelFontHeightF = 0.018
    res@lbTitleFontHeightF = 0.02
    res@lbBoxMinorExtentF = 0.18
    res@gsnSpreadColors  = True
    ;res@cnLevelSelectionMode="ExplicitLevels"
    res@cnLevelSelectionMode="AutomaticLevels"

    minab = min((/min1,min2/))
    maxab = max((/max1,max2/))
    mnmxint = nice_mnmxintvl(minab,maxab,21,False)
    res@cnMinLevelValF = mnmxint(0)
    res@cnMaxLevelValF = mnmxint(1)
    res@cnLevelSpacingF = mnmxint(2)

    res@gsnLeftString = A@long_name
    res@tiMainString = title1
    res@sfYArray = A&lat
    plot(0) = gsn_csm_lat_time(wks,A(lat|:,time|:),res)
    plot(0) = ZeroNegDashLineContour (plot(0))

    delete (res@lbTitleString)
    delete(res@sfYArray)
    delete(res@tiMainString)

    res@tiMainString = title2

    res@sfYArray = B&lat
    if (color_type .eq. "COLOR") then
        res@lbTitleString = "MIN = "+sprintf("%6.2f",min2)+ \
                           " MAX = "+sprintf("%6.2f",max2)
    end if
    plot(1) = gsn_csm_lat_time(wks,B(lat|:,time|:),res)
    if (color_type .eq. "MONO") then
        plot(1) = ZeroNegDashLineContour (plot(1))
    else
        delete (res@cnLevels)
        delete (res@lbTitleString)
    end if

    delete (res@sfYArray)
    delete (res@gsnLeftString)

;----------------------------------------------------------------------
; difference plot of time(x) vs. latitude(y)
; check for different number of latitudes and then
; lineary interpolate to the smaller number of latitudes

    if (nlat1 .ne. nlat2) then                 ; lat grids different
        if (nlat1 .gt. nlat2) then
            tmp1 = linint1 (lat1,A,False,lat2,0)   ; a(time,lat)
            C = B                                  ; copy dims,coords
            C = (/tmp1-B/)                         ; get diff values
            delete (tmp1)
            res@sfYArray = lat2
        else
            tmp2 = linint1 (lat2,B,False,lat1,0)   ; b(time,lat)
            C = A                                  ; copy dims,coords
            C = (/A-tmp2/)                         ; get diff values
            delete (tmp2)
            res@sfYArray = lat1
        end if
    else                                       ; same grid latitudes
        C = A                                    ; copy dims,coords
        C = (/A-B/)                              ; get diff values
        res@sfYArray = lat1
    end if

    mind = min(C)
    maxd = max(C)

    res@tiMainString = "Difference"
    if (color_type .eq. "COLOR") then
        ;res@cnLevels = dcntr(i,:) 
        res@lbLabelStride = 1
        res@lbOrientation = "Vertical"
        res@lbTitleString = "MIN = "+sprintf("%6.2f",mind)+ \
                            " MAX = "+sprintf("%6.2f",maxd)
    else
        mnmxint = nice_mnmxintvl(mind,maxd,dcnlvls,False)
        res@cnMinLevelValF = mnmxint(0)
        res@cnMaxLevelValF = mnmxint(1)
        res@cnLevelSpacingF = mnmxint(2)
    end if

    plot(2) = gsn_csm_lat_time(wks,C(lat|:,time|:),res)
    plot(2) = ZeroNegDashLineContour (plot(2))
    if (color_type .eq. "COLOR") then
        delete (res@lbTitleString)
        delete (res@cnLevels)
    end if
    delete (mind)
    delete (maxd)

    pan = True
    pan@gsnMaximize = True
    pan@gsnFrame = False
    pan@gsnPaperOrientation = "portrait"
    if (time_stamp .eq. "True") then
        pan@gsnPanelBottom = 0.05
        gsn_panel(wks,plot,(/2,2/),pan)
        infoTimeStamp(wks, 0.011, "DIAG Version: "+version)   
    else
        gsn_panel(wks,plot,(/2,2/),pan)
    end if
    frame (wks)

    return True     
end

load "$HOME/HiRAM_ncl/func_read_c360.ncl"
t = read_1month("ts",1999,4)
tA = t(lon| 0:11,lat|:)
tA!0 = "time"
tB = t(lon|12:23,lat|:)
tB!0 = "time"
a = plot_annual_cycle("title1","title2","diffplot",tA,tB) 
