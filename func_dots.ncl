load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("add_dots")
function add_dots(wks,plot,dots[*][*],ores)
begin
    res = ores
    monocolor = False
    if(isatt(dots,"monocolor"))then
        monocolor = dots@monocolor
    end if
    dims = dimsizes(dots)
    if(isatt(ores,"overmark"))then
        overmark = ores@overmark
    else
        overmark = 0.
    end if
    if(isatt(ores,"datamax"))then
        datamax = ores@datamax
        ;print("datamax: "+datamax)
    else
        datamax = 0.
    end if
    if(isatt(ores,"datamin"))then
        datamin = ores@datamin
    else
        datamin = 0.
    end if

    mint= 5
    if (datamax .eq.0)then
        mmi = nice_mnmxintvl(datamin,max(abs(dots)),mint,False)
        datamax = mmi(1)
    else
        mmi = nice_mnmxintvl(datamin,datamax,mint,False)
    end if
    datamin = max((/mmi(0),datamin/))
    dataint = (datamax-datamin)/(mint-1)
    reso = 5.
    nx = dims(1)
    ny = dims(0)
    maxdot = 45.

    gsres               = True
    gsres@gsMarkerIndex = 16 
    gsres@gsMarkerSizeF  = maxdot
    gsres@gsMarkerColor   = "magenta"
    if(monocolor)then
        gsres@gsMarkerColor = "gray40"
    end if

    mgsres = True
    mgsres@gsMarkerIndex  = 3
    mgsres@gsMarkerColor  = "yellow"
    mgsres@gsMarkerSizeF = maxdot/3


    txres               = True
    txres@txFontHeightF = 0.015
    txres@txFont        = 21
    txres@txAngleF      = 15

    do j = 0, ny-1
    do i = 0, nx-1
        gsres@gsMarkerThicknessF  = 1.0
        if(.not.ismissing(dots(j,i)))then
            gsres@gsMarkerSizeF = maxdot * sqrt(abs(dots(j,i)-datamin)/(datamax-datamin))
            if (dots(j,i).lt.0)then
                if(monocolor)then
                    gsres@gsMarkerIndex =  4
                    gsres@gsMarkerThicknessF  = 4.0
                    gsres@gsMarkerColor   = "gray40"
                else
                    gsres@gsMarkerColor   = "blue"
                end if
            else
                if(monocolor)then
                    gsres@gsMarkerIndex = 16 
                    gsres@gsMarkerColor   = "gray40"
                else
                    gsres@gsMarkerColor   = "magenta"
                end if
            end if
            if(abs(dots(j,i)).gt. (dataint/2) ) then
                str = unique_string("aa")
                plot@$str$ = gsn_add_polymarker(wks,plot,dots&lon(i),dots&lat(j),gsres)
                label = ""+dots(j,i)
                ;gsn_text (wks,plot,label,dots&lon(i),dots&lat(j),txres)
                if(overmark.ne.0 .and. abs(dots(j,i)).ge. overmark .and. (abs(dots(j,i)).lt.datamax*100)) then
                    str = unique_string("aa")
                    if(monocolor.and.dots(j,i).gt. 0.)then
                        mgsres@gsMarkerColor  = "white"
                    end if
                    if(monocolor.and.dots(j,i).lt. 0.)then
                        mgsres@gsMarkerColor  = "black"
                    end if
                    ;print("overmark: "+dots(j,i)+" "+overmark)
                    plot@$str$ = gsn_add_polymarker(wks,plot,dots&lon(i),dots&lat(j),mgsres)
                end if
            end if
        else
            str = unique_string("aa")
            plot@$str$ = gsn_add_polymarker(wks,plot,dots&lon(i),dots&lat(j),False)
            
        end if
    end do
    end do

    txres               = True
    txres@txFontHeightF = 0.015
    txres@txFont        = 21
    txres@txAngleF      = 0

    gsres@gsMarkerColor   = "magenta"
    if(monocolor)then
        gsres@gsMarkerColor   = "gray40"
    end if
    dx = .5 / (mint)
    do i = 0,mint-1
        val = datamin+dataint*i
        label = ""+val
        gsres@gsMarkerSizeF = maxdot*sqrt(abs(val)/datamax)
        ;print("val : "+val+" size: "+gsres@gsMarkerSizeF)
        if(val.gt.0 .and. gsres@gsMarkerSizeF.le.maxdot)then
            gsn_polymarker_ndc(wks,0.5+dx*(i-mint/2-.5),.21,gsres)
            gsn_text_ndc      (wks,label,0.5+dx*(i-mint/2-.5),.18,txres)
        end if
        if(gsres@gsMarkerSizeF.gt.maxdot)then
            exit
        end if
    end do

    return plot
end

undef("plot_dot")
function plot_dot(dots[*][*],ores)
begin
    if(isatt(ores,"title"))then
        title = ores@title
    else
        title = "dot_test_title"
    end if
    if(isatt(ores,"filename"))then
        filename = ores@filename
    else
        filename = "dot"
    end if
    if(isatt(ores,"datamin"))then
        datamin = ores@datamin
    else
        datamin = -999.
    end if
    if(isatt(ores,"datamax"))then
        datamax = ores@datamax
    else
        datamax = 0.
    end if
    if(isatt(ores,"cnlev"))then
        cnlev = ores@cnlev
        cndata= ores@cndata
        cndata!0 = "lat"
        cndata!1 = "lon"
        cnlat    = ores@cndatalat 
        cnlon    = ores@cndatalon
        cnlat@units = "degrees_north"
        cnlon@units = "degrees_east"
        cndata&lat = cnlat 
        cndata&lon = cnlon
    else
        cnlev = 0.
    end if

    res            = True
    res@gsnDraw    = False                         ; don't draw
    res@gsnFrame   = False                         ; don't advance frame
    ;res@gsnAddCyclic = False
    res@vpWidthF   = 0.80                          ; make map bigger
    res@vpHeightF  = 0.80
    res@mpMaxLatF  =  40.                           ; select subregion
    res@mpMinLatF  =   0.
    res@mpMinLonF  = 110.
    res@mpMaxLonF  = 180.
    res@mpCenterLonF = 180.
    res@mpFillDrawOrder = "PreDraw"
    res@mpFillOn = False
    if(isatt(ores,"plotarea").and.ores@plotarea.eq."Pac")then
        res@mpMinLatF  =  -40.
        res@mpMaxLonF  = 300.
    end if
    res@tiMainString = title
    res@tiMainJust = "CenterLeft"
    res@tiMainPosition      = "Left"
    if(isatt(ores,"gsnRightString"))then
        res@gsnRightString = ores@gsnRightString
        res@tiMainString = ""
        res@gsnLeftString  = title
        res@gsnLeftStringFontHeightF = 0.018
        res@gsnRightStringFontHeightF = 0.018
        res@gsnStringFont  = 21
    end if
    if(cnlev .gt. 0.)then
        ;res@cnLevelSelectionMode = "ExplicitLevels"
        ;res@cnLevels = (/5820., 5880., 5920, 5950/)  ; for z 500 hPa
        res@cnLevelSpacingF = 5.
        res@gsnContourZeroLineThicknessF = 3.0
        res@gsnContourNegLineDashPattern = 1
    end if
    ;if(isatt(ores,"tiMainPosition"))then
    ;    res@tiMainPosition  = ores@tiMainPosition 
    ;    res@tiMainJust = "CenterLeft"
    ;end if

    reso = 5.
    wks = gsn_open_wks("ps",filename)
    print("plot dot: "+filename)
    if(cnlev .gt. 0.)then
        ;print("cn on dots")
        res@cnLevelSpacingF = 0.5
        ;print(max(cndata))
        ;print(min(cndata))
        plot = gsn_csm_contour_map(wks,cndata,res)
        ares = True
        ares@cnLevelSelectionMode = "ExplicitLevels"  ;; add 5880 thickness
        ares@cnLevels = 0.  ; for z 500 hPa
        ares@cnLineThicknessF = 3.0
        ares@cnLineLabelsOn = False
        ares@cnInfoLabelOn  = False
        ares@gsnDraw    = False                         ; don't draw
        ares@gsnFrame   = False                         ; don't advance frame
        ;plot5880 = gsn_csm_contour(wks,cndata,ares)
        ;overlay(plot,plot5880)
    else
        plot = gsn_csm_map(wks,res)
    end if
    plot1 = add_dots(wks,plot,dots,ores)
    if(isatt(dots,"rectangle"))then
        ;; bottom left corner start clockwise
        if(dots@rectangle .eq."MDR")then
            linex = (/110.,110.,150.,150.,110./)
            liney = (/ 10., 30., 30., 10., 10./)
        end if
        if(dots@rectangle .eq."MGR")then
            linex = (/110.,110.,160.,160.,110./)
            liney = (/ 10., 25., 25., 10., 10./)
            totalwnp = sum(dots({ 0:40},{110:180}))
            totalmgr = sum(dots({10:25},{110:160}))
            mgrratio = totalmgr/totalwnp
            print("mgr/wnp = "+totalmgr+"/"+totalwnp+" = "+mgrratio)
        end if
        mdrres = True
        mdrres@gsLineColor  ="red" 
        mdrres@gsLineThicknessF = 5.0    
        gsn_polyline(wks,plot,linex,liney,mdrres)
    end if
    draw(plot1)
    frame(wks)
    return True
end



undef("plot_ssta_wpsh_race_wind")
function plot_ssta_wpsh_race_wind(ssta,z[*][*],race[*][*],u[*][*],v[*][*],ores)
begin
    if(isatt(ores,"title"))then
        title = ores@title
    else
        title = "dot_test_title"
    end if
    if(isatt(ores,"filename"))then
        filename = ores@filename
    else
        filename = "testHybridFig"
    end if
    if(isatt(ores,"dotmax"))then
        dotmax = ores@dotmax
    else
        dotmax = 0.
    end if
    
    wks = gsn_open_wks("ps",filename)
    print("ploting... "+filename)

    ;; base map with ssta, but ssta is not imply yet.
    res1 = True
        res1@gsnDraw    = False
        res1@gsnFrame   = False
        res1@mpMaxLatF  =  40.                           ; select subregion
        res1@mpMinLatF  =  00.
        res1@mpMinLonF  = 110.
        res1@mpMaxLonF  = 180.
        res1@mpCenterLonF = 180.
        res1@mpFillOn   = False
        res1@tfPolyDrawOrder  = "Draw"
        res1@gsnLeftString = ""
        res1@gsnRightString = ""
        res1@tiMainString  = title
        res1@tiMainJust = "CenterLeft"
        res1@tiMainPosition      = "Left"


    ;;plot1 = gsn_csm_contour_map_ce(wks,ssta,res1) 
    plot1 = gsn_csm_map_ce(wks,res1) 
    
    ;; WPSH for overlay
    res2 = True
        res2@gsnDraw    = False
        res2@gsnFrame   = False
        res2@cnLevelSelectionMode = "ExplicitLevels"
        res2@cnLevels = (/1520,1540/)
        res2@cnFillOn = True
        res2@cnFillDrawOrder = "PreDraw"
        res2@cnLinesOn = False
        res2@cnFillColors = (/"transparent","yellow","yellowgreen"/)
        res2@lbOrientation = "Vertical"
        res2@lbLabelFontHeightF =  0.01
        res2@lbLeftMarginF  = 0.
        res2@lbRightMarginF  = 0.3
        res2@gsnLeftString = ""
        res2@gsnRightString = ""
    
    plot2 = gsn_csm_contour(wks,z,res2)
    
    ;; wind vector for overlay
    res3 = True
        res3@gsnDraw    = False
        res3@gsnFrame   = False
        res3@gsnLeftString = ""
        res3@gsnRightString = ""
        res3@vcRefMagnitudeF  = 15.
        ;res3@vcRefLengthF     =   .10
        res3@vcRefLengthF     =   .02
        res3@vcGlyphStyle = "WindBarb"
        res3@vcRefAnnoOn  = False
        u = u/0.514  ; m/s to knots
        v = v/0.514
    
    plot3 = gsn_csm_vector(wks,u,v,res3)
    
    ;; overlay and add RACE dots
    res1@datamax = dotmax
    plot4 = add_dots(wks,plot1,race,res1)
    overlay(plot1,plot2)
    overlay(plot1,plot3)
    
    draw(plot1)
    frame(wks)
    return True
end

