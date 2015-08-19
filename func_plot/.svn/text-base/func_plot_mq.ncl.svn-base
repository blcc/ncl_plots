load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "func_read_ylchen_mmq.ncl"
load "func_read_mq.ncl"
load "func_plot_time_series.ncl"
load "func_read_rvke_sst.ncl"

undef("annual_sum3d")  ;; 12t,y,x
function annual_sum3d(monts[*][*][*],mon[*])
begin
    time = monts&$monts!0$
    years = toint(time(::12)/100)
    if(.not.any(mon.eq.0))then
        yrts = monts(::12,:,:)
        yrts = 0
        do m = 0,11
            mm = m+1
            if(any(mm.eq.mon))then
                yrts = yrts + monts(m::12,:,:)
            end if
        end do
        nm = dimsizes(mon)
    else
        yrts = month_to_annual(ts,0)
        nm = 12
    end if
    yrts@nmonth = nm
    yrts!0 = "year"
    yrts&year = years
    return yrts
end

function plotsquare(wks,plot)
begin
    ;; bottom left corner start clockwise
        ;;ts = wgt_areaave(gpic(:,{10:25},{110:160}),1.,1.,0) ;; MGR
    linex = (/110.,110.,160.,160.,110./)
    liney = (/ 10., 25., 25., 10., 10./)
    mdrres = True
    mdrres@gsLineColor  ="red" 
    mdrres@gsLineThicknessF = 5.0    
    a = gsn_add_polyline(wks,plot,linex,liney,mdrres)
    return a
end

undef("read_ylchen_mq_clm")
function read_ylchen_mq_clm(years,mons)
begin
    mqmonly = read_ylchen_mq(years)
    mqclm   = clmMonTLL(mqmonly)
    mq = dim_avg_n_Wrap(mqclm(mons-1,:,:),0)
    return mq
end

;undef("plot_mq_clm")
function plot_mq_clm(years[*],imons[*],title[1]:string,filename[1]:string)
begin
    if(any(imons.eq.0))then
        mons = ispan(1,12,1)
        print("mons content 0, set to 1-12")
    else
        mons = imons
    end if
    if(isatt(years,"region"))then
        region = years@region
    else
        region = "tropical"
    end if

    mmq = read_mmq_clm(years,mons)
    m  = mmq
    mq = mmq@mq
    copy_VarCoords(m,mq)
    print("m/mq max = "+max(m)+" "+max(mq))
    wks = gsn_open_wks("ps",filename)
    gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

    res                     = True
        res@gsnDraw = False
        res@gsnFrame = False
        res@vpHeightF            = 0.40           ; Changes the aspect ratio
        res@vpWidthF             = 0.85
        res@tmXTOn              = False
        res@tmYROn              = False
        ;if(region.eq."tropical")then  ;; default
            res@mpMinLatF            = -40.
            res@mpMaxLatF            =  40.
            res@mpMinLonF            =   0.
            res@mpMaxLonF            = 360.
        ;end if
        if(region .eq."WNP")then
            res@mpMinLatF            =   0.
            res@mpMaxLatF            =  40.
            res@mpMinLonF            = 110.
            res@mpMaxLonF            = 180.
        end if
    
        res@mpCenterLonF         = 180.
    
        res@cnLinesOn           = True
        res@cnFillOn            = True            ; turn on color fill
        res@tiMainString        = title
    
        res@tiMainFont     = 29
        res@txFont         = 29
        res@tmXBLabelFont  = 29
        res@tmYLLabelFont  = 29
    
        res@gsnMajorLatSpacing  = 10
        res@gsnMajorLonSpacing  = 60
    
        res@tiMainFontHeightF       = 0.015
        res@gsnStringFontHeightF    = 0.012
        res@tmXBLabelFontHeightF    = 0.010
        res@tmYLLabelFontHeightF    = 0.010
        res@tmXBLabelFontHeightF    = 0.010
        res@tmYLLabelFontHeightF    = 0.010
    
        res@lbLabelBarOn             = True        ; turn on individual lb's
        res@lbLabelFont              = 21
        res@lbLabelFontHeightF       = 0.010
        res@pmLabelBarOrthogonalPosF = 0.2
        res@pmLabelBarWidthF         = 0.85
        res@pmLabelBarHeightF        = 0.04
    
        res@cnLevelSelectionMode  = "ExplicitLevels"
        res@cnLevels              = ispan(800,2400,200)
        res@cnFillColors          = (/ 50,  62,  76,   83,  86,  92,  95,  97,  100,  102 ,102/)
        if(isatt(years,"sqRoot"))then
            delete(res@cnLevels)
            delete(res@cnFillColors)
            res@cnLevels          = ispan( 0,  60, 5)
            res@cnFillColors          = (/ 50,56,  62,69,76,80,83,84,86,  92,  95,  97,  100,  102 ,102/)
        end if

    plot = gsn_csm_contour_map_ce(wks,mq,res)
    plotsq = plotsquare(wks,plot)
    draw(plot)
    frame(wks)
    print("plot: "+filename)
    return mq
end

function plot_mq_ts(years[*],imons[*],title[1]:string,filename[1]:string)
begin
    yb = min(years)
    ye = max(years)
    if(any(imons.eq.0))then
        mons = ispan(1,12,1)
        print("mons content 0, set to 1-12")
    else
        mons = imons
    end if
    if(isatt(years,"range"))then
        range = years@range
    else
        range = "WNP"
    end if

    gen  = readmonTygen(yb,ye)
    genyr = annual_sum3d(gen,mons)
    genyr = genyr/genyr@nmonth

    mmon = read_mmq_monly(years)
    mqmon = mmon@mq
    copy_VarCoords(mmon,mqmon)

    mqyr = annual_sum3d(mqmon,mons)
    mqyr = mqyr/mqyr@nmonth

    if(range.eq."MGR")then
        genyrts = dim_avg_n_Wrap(genyr(:,{10:25},{110:160}),(/1,2/))  ;; MGR
        mqyrts = dim_avg_n_Wrap(mqyr(:,{10:25},{110:160}),(/1,2/))  ;; MGR
    end if
    if(range.eq."WNP")then
        genyrts = dim_avg_n_Wrap(genyr(:,{0:40},{110:180}),(/1,2/))  ;; WNP
        mqyrts = dim_avg_n_Wrap(mqyr(:,{0:40},{110:180}),(/1,2/))  ;; WNP
    end if

    ;myr = annual_sum3d(mmon,mons)
    ;myr = myr/myr@nmonth
    ;myrts = dim_avg_n_Wrap(myr(:,{10:25},{110:160}),(/1,2/)) ;; MGR

    dims = dimsizes(mqyrts)
    
    mmqyrts = new((/2,dims/),typeof(mqyrts))
    mmqyrts(0,:) = genyrts
    mmqyrts(1,:) = mqyrts
    mmqyrts!0 = "mmq"
    mmqyrts&mmq = (/"NTC","Mq"/)

    copy_VarAtts(years,mmqyrts)
    ;; standardize
    mmqyrts = dim_standardize_n_Wrap(mmqyrts,1,1)
    mmqyrts@ymin = -2.
    mmqyrts@ymax =  3.

    a = plot_allts(mmqyrts,title,filename)
    return a
end

    ;; test
    ;a = plot_mq_clm(ispan(1988,2001,1),0,"Gross Moisture Stratification (Mq)","Mq_clm_1988-2001")
    ;a = plot_mq_ts(ispan(1988,2001,1),0,"Gross Moisture Stratification (Mq)","Mq_clm_1988-2001")

