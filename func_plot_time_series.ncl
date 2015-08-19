load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" ;; for drawNDCGrid()
undef("plot_mcol_bar")
function plot_mcol_bar(wks,dd[*][*],res) ;plot multi color bar chart
begin
    dims    = dimsizes(dd)
    ncolor  = dims(0)+1
    nt      = dims(1)
    lplot   = new(ncolor,"graphic")
    lres    = res
    lres@tiYAxisString = ""
    lres@gsnDraw = False
    lres@gsnFrame= False
    lres@gsnXYBarChartColors = 9
    lres@gsnYRefLine   = 0
    lres@xyCurveDrawOrder = "PreDraw"
    if (isatt(dd,"nfft"))then
        nfft = dd@nfft
    else
        nfft = 0
    end if
    if (isatt(dd,"rave"))then
        rave = dd@rave
    else
        rave = 0
    end if
    if(isatt(res,"yavgref"))then
        delete(lres@yavgref)
    end if
    ldd     = dd
    xcord   = ldd!1
    x       = ldd&$xcord$
    ;; sum, background bars must higher than front
    do i = 0, ncolor-2
    do j = i+1, ncolor-2
        ldd(i,:) = ldd(i,:) + ldd(j,:)
    end do
    end do
    ;lres@gsnRightString = "avg = "+sprintf("%4.2f",avg(ldd(0,:)))
    if (isatt(dd,"LeftString"))then
        lres@gsnLeftString = dd@LeftString
    end if
    if (isatt(dd,"CenterString"))then
        lres@gsnCenterString = dd@CenterString
    end if
    if (max(ldd).le.101. .and. max(ldd).ge.99.)then
        print(max(ldd))
        lres@trYMaxF = 100.
        lres@gsnRightString = ""
    end if

    lplot(0) = gsn_csm_xy(wks,x,ldd(0,:),lres)
    do i = 1, ncolor-2
        lres@gsnXYBarChartColors = lres@gsnXYBarChartColors+3
        lplot(i) = gsn_csm_xy(wks,x,ldd(i,:),lres)
        overlay(lplot(0),lplot(i))
    end do
        addp = ldd(0,:)  ;; for some weird behavior: addline cause last bar outline missing.
        addp = 0.
        lplot(ncolor-1) = gsn_csm_xy(wks,x,addp,lres)
        overlay(lplot(0),lplot(ncolor-1))
    yavg   = avg(ldd(0,:))
    yref   = avg(ldd(0,:))
    reflines = True
    reflines@tfPolyDrawOrder = "PostDraw"
    if(isatt(res,"yavgref").and.res@yavgref)then
        plot2 = gsn_add_polyline(wks,lplot(ncolor-1),(/0,max(ldd&$ldd!1$)+1/),(/yref,yref/),reflines)
    end if

    if (isatt(dd,"stddevscale"))then
        stddevscale = dd@stddevscale
    else
        stddevscale = 1.0
    end if
    if (isatt(dd,"nostddev").and.dd@nostddev)then
        
    else 
        stdlevel = 1.0 * stddevscale
        ydev   = stdlevel*stddev(ldd(0,:))
        stdval = (ldd(0,:)-avg(ldd(0,:)))/ydev
            print("stddev line: "+stdlevel+"*stddev")
            ;print(ldd&$ldd!1$+" "+stdval)
        yref   = avg(ldd(0,:)) - ydev
        plot3 = gsn_add_polyline(wks,lplot(ncolor-1),(/0,max(ldd&$ldd!1$)+1/),(/yref,yref/),reflines)
        yref   = avg(ldd(0,:)) + ydev
        plot4 = gsn_add_polyline(wks,lplot(ncolor-1),(/0,max(ldd&$ldd!1$)+1/),(/yref,yref/),reflines)
    end if
    if (isatt(dd,"addline"))then
        addline = dd@addline
        if(isatt(dd,"addlineyb"))then
            addlineyb = dd@addlineyb
            addlineye = addlineyb + dimsizes(addline)
        else
            addlineyb = min(x)
            addlineye = max(x)
        end if
        alavg   = avg(addline)
        alstd   = 0.5 ; stddev(addline) ; 0.5 for CPC's ONI
        ;; standardlize addline time series and fit to main plot
        ;addline = ((addline-alavg)/alstd * ydev)  +yavg
        print("line avg/std: "+alavg+"/"+alstd)

        sstline = True
        sstline@gsLineColor             = "red"
        sstline@gsLineThicknessF    = 3.0
        sstline@tfPolyDrawOrder = "PostDraw"
        ;plot5 = gsn_add_polyline(wks,lplot(ncolor-1),ispan(addlineyb,addlineye,1),addline,sstline)
        plot5 = gsn_add_polyline(wks,lplot(ncolor-1),x,addline,sstline)
        ;cor = escorc(addline,ldd(0,:))
        ;print("corr: "+cor)
    end if
    if(nfft .gt. 0)then
        dtsfft = ezfftf(ldd(0,:))
        dtsfft(:,nfft:) = 0
        fftres = ezfftb(dtsfft,dtsfft@xbar)
        plotfft = gsn_add_polyline(wks,lplot(ncolor-1),x,fftres,False)
    end if
    if(rave .gt. 0)then
        raves = runave(ldd(0,:),rave,0)
        avelineres =  True
        avelineres@gsLineThicknessF = 2.0
        plotrave = gsn_add_polyline(wks,lplot(ncolor-1),x,raves,avelineres)
    end if
    draw(lplot(0))
    frame(wks)
    return True
end


undef("plot_time_series")
function plot_time_series(dts,title,filename)
begin
    ;; dts(basin,time)
    if(isatt(title,"RightMean"))then
        RightMean = title@RightMean
    else 
        RightMean = False
    end if
    if(isatt(title,"RightSum"))then
        RightSum = title@RightSum
    else
        RightSum = False
    end if
    if(isatt(title,"noRight"))then
        noRight = title@noRight
    else
        noRight= False
    end if
    dims = dimsizes(dts)
    if (isatt(dts,"nfft"))then
        nfft = dts@nfft
    else
        nfft = 0
    end if
    if (isatt(dts,"ymin"))then
        ymin = dts@ymin
    else
        ymin = min(dts)
    end if
    if(isatt(title,"plotline"))then
        plotline = title@plotline
    else
        plotline = False
    end if

    res = True
    ;;res@tiMainString = title
    res@tiMainString = ""
    res@gsnLeftString = title

    if(isatt(title,"tiMainPosition"))then
        res@tiMainPosition = title@tiMainPosition
        res@tiMainJust = title@tiMainJust
    else
        res@tiMainJust = "TopLeft"
        res@tiMainPosition = "Left"
    end if

    if(plotline)then
        res@gsnXYBarChart   = False
    else
        res@gsnXYBarChart   = True
    end if

    res@vpWidthF   = 0.8
    res@vpHeightF   = 0.3
    res@txFont = 21
    res@tiMainFont      = 21
    res@tiMainFontHeightF = 0.018
    res@gsnLeftStringFontHeightF = res@tiMainFontHeightF
    res@tiXAxisFont      = 21
    res@tiXAxisFontHeightF = 0.018
    if(isatt(dts,"ymin"))then
        res@trYMinF   = dts@ymin
    end if
    if(isatt(dts,"ymax"))then
        res@trYMaxF   = dts@ymax
    end if
    if(RightMean)then
        dtsmean = avg(dim_sum_n(dts,0))
        res@gsnRightString = "Mean="+sprintf("%4.2f",dtsmean)
        if(dtsmean.lt.0.00001)then
            res@gsnRightString = ""
        end if
        res@gsnRightStringFontHeightF = res@tiMainFontHeightF
    end if
    if(RightSum)then
        dtssum = sum(dim_sum_n(dts,0))
        res@gsnRightString = "Sum="+sprintf("%4.2f",dtssum)
        res@gsnRightStringFontHeightF = res@tiMainFontHeightF
    end if
    if(isatt(dts,"tiXAxisString"))then
        res@tiXAxisString = dts@tiXAxisString
    end if
    if(noRight)then
        res@gsnRightString = ""
    end if


    ;res@tmYLFormat =  "0@*+^sg" 

    print("plot: "+filename)
    wks  = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"nrl_sirkes")
    if (dimsizes(dims).eq.1.or. plotline)then
        if(plotline .or. (dims .ge.15 .and. dims.ne.73))then 
            stdlevel  = 1.
            if(isatt(dts,"yref"))then
                res@gsnYRefLine   = dts@yref
            else
                res@gsnYRefLine   = avg(dts)
            end if
            if(RightMean .or. RightSum.or.noRight)then
            else
                print("stddev line: "+stdlevel+"*stddev")
                yref1 = avg(dts) + stdlevel*stddev(dts)
                yref2 = avg(dts) - stdlevel*stddev(dts)
                ;res@gsnRightString = "Mean="+sprintf("%4.2f",avg(dts))+" stddev = "+sprintf("%4.2f",stddev(dts))
                res@gsnRightString = " stddev = "+sprintf("%4.2f",stddev(dts))
            end if
            res@gsnAboveYRefLineColor = 15
            res@gsnBelowYRefLineColor = 15
            res@gsnDraw = False
            res@gsnFrame= False
        else
            res@gsnYRefLine   = 0
            res@gsnXYBarChartColors = 9
            res@tmLabelAutoStride = False
            res@tmXBMaxTicks = 13
            res@tmXBMinorOn = False
            res@tmXBPrecision = 2
        end if
        res@tiXAxisString   = ""
        if(dimsizes(dims).eq.1)then
            plot1 = gsn_csm_xy(wks,dts&$dts!0$,dts,res)        ; create histogram 
        else
            delete(res@gsnYRefLine)
            delete(res@gsnAboveYRefLineColor)
            delete(res@gsnBelowYRefLineColor)
            plot1 = gsn_csm_xy(wks,dts&$dts!1$,dts,res)        ; create line
        end if
        if(nfft .gt. 0)then
            dtsfft = ezfftf(dts)
            dtsfft(:,nfft:) = 0
            fftres = ezfftb(dtsfft,dtsfft@xbar)
            plotfft = gsn_add_polyline(wks,plot1,dts&$dts!0$,fftres,False)
        end if
        if(isdefined("yref1"))then
            plot2 = gsn_add_polyline(wks,plot1,(/0,max(dts&$dts!0$)+1/),(/yref1,yref1/),False)
            plot3 = gsn_add_polyline(wks,plot1,(/0,max(dts&$dts!0$)+1/),(/yref2,yref2/),False)
            plot4 = gsn_add_polyline(wks,plot1,(/0,max(dts&$dts!0$)+1/),(/avg(dts),avg(dts)/),False)
            delete(yref1)
        end if
        draw(plot1)
        frame(wks)
    else
        res@yavgref = True
        plot1 = plot_mcol_bar(wks,dts,res) ;plot multi color bar chart
        ;print("plot2")
    end if

    return True
end

undef("plot_allts_res")
function plot_allts_res(allyrts[*][*],title[1]:string,filename[1]:string,ires[1]:logical)
begin  ;; allyrts(kind,time) , res not imply yet(need to be done)
    res = ires
    if(isatt(allyrts,"plotmode"))then  ;; colors, lines(dash) or fft setting
        plotmode = allyrts@plotmode
    else
        plotmode = "none"
    end if
    if(isatt(res,"addbackcolorx"))then ;; (nblock,(/xstart,xend/))(n,2)  fill color @ background
        ;; backcolor example
        ;; res@backcolorx = (/(/x11,x12/),(/x21,x22/),(/x31,x32/)/)
        ;; res@backcolory = (/(/y11,y12/),(/y21,y22/),(/y31,y32/)/) optional

        backcolorx = res@addbackcolorx
        delete(res@addbackcolorx)
        if(isatt(res,"addbackcolory"))then ;; (nblock,(/ystart,yend/)) y range, optional
            backcolory = res@addbackcolory
            delete(res@addbackcolory)
        end if
    end if
    if(isatt(res,"addlinesx"))then ;; (nline,xpts),(nline,ypts)  add lines
        linex = res@addlinesx
        liney = res@addlinesy
        delete(res@addlinesx)
        delete(res@addlinesy)
    end if

    nts = dimsizes(allyrts&$allyrts!0$)
    ntsx = dimsizes(allyrts&$allyrts!1$)

    wks = gsn_open_wks("ps",filename)
    res = True
        res@tfPolyDrawOrder    = "Predraw" ;; if background color
        res@gsnFrame = False
        res@gsnDraw  = False
        res@vpWidthF   = 0.8
        res@vpHeightF   = 0.3
        res@tiMainFontHeightF = 0.018
        res@gsnLeftString = title
        res@gsnLeftStringFontHeightF = 21
    
        if(.not.isatt(res,"pmLegendDisplayMode"))then
            res@pmLegendDisplayMode    = "Always" 
        end if
        res@pmLegendParallelPosF   = 0.10                  ; move units right
        res@pmLegendOrthogonalPosF = -1.1 ;-0.5 ; more neg = up
        res@pmLegendWidthF         = 0.10                  ; Change width and
        res@pmLegendHeightF        = 0.08                  ; height of legend.
        res@lgPerimOn              = False                 ; no box around
        res@lgLabelFontHeightF     = .01                   ; change font height
        res@xyExplicitLegendLabels = allyrts&$allyrts!0$
    
        res@trXMinF = min(allyrts&$allyrts!1$)-0.1
        res@trXMaxF = max(allyrts&$allyrts!1$)+0.1
        if(isatt(allyrts,"ymin"))then
            res@trYMinF = allyrts@ymin
        end if
        if(isatt(allyrts,"ymax"))then
            res@trYMaxF = allyrts@ymax
        end if
        if(ntsx.eq.12)then
            res@tmXBMode  = "Explicit"
            res@tmXBValues = ispan(1,12,1)
            res@tmXBLabels = ispan(1,12,1)
        end if
        if(.not.isatt(res,"xyLineColors"))then
            res@xyLineColors      = (/"black","blue","green","red", "black", "RoyalBlue", "lightseagreen","Pink"/)   
            res@xyLineThicknesses = (/ 2.0   , 2.0  , 2.0 , 2.0    , 3.0   , 2.0        , 2.0            ,2.0/)
        end if
    if(plotmode.eq."fft")then
        cf = ezfftf(allyrts)
        cf(:,:,5:) = 0.     ;; about 9yr smoothing for 44 yr time series.
        fftyrts = ezfftb(cf,cf@xbar)

        delete(res@xyLineThicknesses)
        ;delete(res@xyLineColors)
        res@xyLineThicknessF  = 1.0
        res@xyDashPattern  = 1
    end if
    if(plotmode.eq."avgts")then
        avgts = dim_avg_n_Wrap(allyrts,0)
        delete(res@xyLineThicknesses)
        delete(res@xyLineColors)
        res@xyLineColor = "gray68"
        res@xyLineThicknessF  = 1.0
        res@xyDashPattern  = 0
        res@pmLegendDisplayMode    = "Never" 
    end if
    if(nts.eq.2.and.ntsx.gt.12.and.(.not.(isatt(allyrts,"nocor").and.allyrts@nocor)).and.(abs(rr).ne.1))then
        rr = escorc(allyrts(0,:),allyrts(1,:))
        res@gsnRightString = "r^2 = "+sprintf("%5.4g",rr)
    end if


    res@gsnYRefLine  = 0.

    print("plot xy lines via plot_allts_res(): "+filename)
    plot  = gsn_csm_xy (wks,allyrts&$allyrts!1$,allyrts,res)
    if(isvar("backcolorx"))then
        gsres = True
            gsres@gsFillColor       = "SlateBlue"    
        bdims = dimsizes(backcolorx)
        do i = 0, bdims(0)-1
            lefx = min(backcolorx(i,:))
            rigx = max(backcolorx(i,:))
            if(isvar("backcolory"))then
                topy = max(backcolory(i,:))
                boty = min(backcolory(i,:))
            else
                topy = max(allyrts)
                boty = min(allyrts)
            end if
            pgx = (/lefx,lefx,rigx,rigx/)
            pgy = (/boty,topy,topy,boty/)
            if(.not.any(ismissing(pgx)))then
                gsn_polygon(wks,plot,pgx,pgy,gsres)
            end if
        end do
    end if
    if(isvar("linex"))then
        linedim = dimsizes(linex)
        lineres = True
        do i = 0, linedim(0)-1
            if(.not.all(ismissing(linex(i,:))))then
                gsn_polyline(wks,plot,linex(i,:),liney(i,:),lineres)
            end if
        end do
    end if
    if(plotmode.eq."fft")then
        plotfft = True
        ares = True
        ares@gsLineDashPattern   = 0
        ares@gsLineThicknessF  = 3.0
        do n = 0,nts-1
            str = unique_string("polyline")
            ares@gsLineColor  = res@xyLineColors(n)
            plotfft@$str$ = gsn_add_polyline (wks,plot,allyrts&$allyrts!1$,fftyrts(n,:),ares)
        end do
    end if
    if(plotmode.eq."avgts")then
        ares = True
        ares@gsLineDashPattern   = 0
        ares@gsLineThicknessF  = 3.0
        plotavg = gsn_add_polyline (wks,plot,allyrts&$allyrts!1$,avgts,ares)

        modeltrend = True
        if(modeltrend)then
            tyb = 2002
            trend = dtrend(avgts({tyb:}),True)
            print("slope: "+trend@slope)
            endy  = (2100-tyb+1)*trend@slope
            ares@gsLineThicknessF  = 1.0
            ares@gsLineColor = "red"
            plottrend = gsn_add_polyline (wks,plot,(/tyb,2100/),(/0,endy/),ares)
        end if
        if(isatt(allyrts,"text"))then
            ;;drawNDCGrid(wks) ; Draw helpful grid lines showing
            tres       =  True
            tres@txFontHeightF = 0.015
            tres@txJust  = "CenterLeft"
            gsn_text_ndc(wks,allyrts@text,0.2,0.75,tres)
        end if
    end if
    ;;drawNDCGrid(wks) ; Draw helpful grid lines showing
    draw(plot)
    frame(wks)


    return allyrts
end

undef("plot_allts")
function plot_allts(allyrts[*][*],title[1]:string,filename[1]:string) ;; 
begin  ;; allyrts(kind,time)
    ;; transtional frontend interface for plot_all_res()
    res = True
    return plot_allts_res(allyrts,title,filename,res)
end
