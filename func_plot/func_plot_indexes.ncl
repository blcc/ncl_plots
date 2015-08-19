load "func_global_rvke.ncl"
load "func_plot_time_series.ncl"
load "func_dots.ncl"
load "func_ttest.ncl"
load "func_read_rvke.ncl"
load "func_read_jma.ncl"
load "func_read_ibtracs.ncl"

undef("read_jtwc_data")
function read_jtwc_data()
begin
    jf = addfile("/walker2/pgchiu/data/tcBestTrackJTWC/ncl/JTWC_wp.nc","r")
    a  = jf->recTy
    return a
end

undef("read_besttrack_data")
function read_besttrack_data(dataset[1]:string)
begin
    ;; default(not jma) is jtwc
    print("reading dataset: "+dataset)
    if(dataset.eq."jma")then
        a = read_jma_textdata()
        return a
    end if
    if(dataset.eq."ibtracs")then
        ;a = read_ibtracs("/walker2/pgchiu/data/IBTrACS/Allstorms.ibtracs.v02r01.nc")
        ;a = read_ibtracs("/walker2/pgchiu/data/IBTrACS/Basin.WP.ibtracs_all.v03r02.nc")
        a = read_ibtracs("/walker2/pgchiu/data/IBTrACS/Basin.WP.ibtracs_wmo.v03r02.nc")
        return a
    end if
    a = read_jtwc_data()
    return a
end

undef("print_cats_stddev")
function print_cats_stddev(yrcats)
begin
    years = yrcats&$yrcats!1$
    cat01 = yrcats(0,:)
    cat23 = yrcats(1,:)
    cat45 = yrcats(2,:)
    allc  = cat01+cat23+cat45
    ;; to ratio
    cat01 = cat01/allc
    cat23 = cat23/allc
    cat45 = cat45/allc

    cat01 = (cat01 -avg(cat01))/stddev(cat01)
    cat23 = (cat23 -avg(cat23))/stddev(cat23)
    cat45 = (cat45 -avg(cat45))/stddev(cat45)
    print(years+" "+cat01+" "+cat45)
    print("corela: "+escorc(cat01,cat45))
    return True
end
undef("print_cats_percent")
function print_cats_percent(yrcats)
begin
    years = yrcats&$yrcats!1$
    cat01 = yrcats(0,:)
    cat23 = yrcats(1,:)
    cat45 = yrcats(2,:)
    allc  = cat01+cat23+cat45
    ;; to ratio
    rcat01 = sum(cat01)/sum(allc) *100
    rcat23 = sum(cat23)/sum(allc) *100
    rcat45 = sum(cat45)/sum(allc) *100

    print("percent of  cat01 cat23 cat45: "+rcat01+" "+rcat23+" "+rcat45)
    return True
end

undef("plot_index_ts")
function plot_index_ts(yb,ye,iindex,title,filename)
begin
    if(isatt(yb,"dataset"))then
        dataset = yb@dataset
    else
        dataset = "jtwc"
    end if
    a  = read_besttrack_data(dataset)
    if (isatt(yb,"emanuel2005Correcting").and.yb@emanuel2005Correcting)then
        print("emanuel2005Correcting")
        a = emanuelCorrecting(a) 
    end if
    a@todo = iindex
    index = global_MonRvKE(a)
    dims = dimsizes(index)
    if(False .and. iindex .eq."gen")then
        ymb = yb*100+1
        yme = ye*100+12
        gents = dim_sum_Wrap(dim_sum_Wrap(index({ymb:yme},:,:)))
        genyr = month_to_annual(gents,0)
        yrs   = gents&time(::12)/100
        genyr!0 = "year"
        genyr&year = yrs
        e = plot_time_series(genyr,title,filename)
        return True
    end if
    
    cats = new((/3,dims(0),dims(1),dims(2)/),"float")
    cats = 0.
    cats!0 = "cat"
    cats&cat = (/0,2,4/)
     
    ;; cat. 0+1, 2+3, 4+
    thresholds = (/ 34., 83.,114./)
    limits     = (/ 82.,113.,999./)
    if(isatt(yb,"RPDIpaper"))then
        thresholds = (/ 34., 64., 96./)
        limits     = (/ 63., 95.,999./)
    end if
    do i = 0,2
        a@threshold = thresholds(i)
        a@limit     = limits(i)
        cats(i,:,:,:) = global_MonRvKE(a)
    end do
    if (isatt(yb,"mon"))then
        mon = yb@mon
        do m = 0,11
            if(.not.any(mon.eq.(m+1)))then
                cats(:,m::12,:,:) = 0
            end if
        end do
    end if
    if (isatt(yb,"scs").and.yb@scs)then
        catsts = dim_sum_Wrap(dim_sum_Wrap(cats(:,{yb*100+1:ye*100+12},{0:60},{90:125}))) ;; scs only
    else
        ;;catsts = dim_sum_Wrap(dim_sum_Wrap(cats(:,{yb*100+1:ye*100+12},{0:60},{90:210}))) ;; wnp only
        catsts = dim_sum_Wrap(dim_sum_Wrap(cats(:,{yb*100+1:ye*100+12},:,:))) ;; all basins
    end if
    catstsat = month_to_annual(catsts(time|:,cat|:),0)
    catstsa  = catstsat(cat|:,year|:)
    catstsa&$catstsa!1$ = ispan(yb,ye,1)
    if (isatt(yb,"smooth121").and.yb@smooth121)then ;; from Emanuel2005
        catstsa = wgt_runave_Wrap(catstsa,(/1,2,1/),0)
        catstsa = wgt_runave_Wrap(catstsa,(/1,2,1/),0)
        scale = 1.09/sum(catstsa(cat|:,{year|1996}))
        scale = 2.8*10^-9
        scale = 1.
        print("scale: "+scale)
        catstsa = catstsa*scale
    end if
    ii = iindex
    if (     ii.eq."ke".or.ii.eq."ace".or.ii.eq."ACE" \
        .or. ii.eq."rvke".or.ii.eq."race".or.ii.eq."RACE")then  ; for ACE and RACE
        catstsa = catstsa /10000
    end if
    if (ii.eq."pdi".or.ii.eq."PDI".or.ii.eq."rpdi".or.ii.eq."RPDI")then  ; PDI and RPDI
        catstsa = catstsa /1000000000
    end if
    if (isatt(yb,"nostddev"))then
        catstsa@nostddev = yb@nostddev
    end if
    ;c = print_cats_stddev(catstsa)
    if(isatt(yb,"ratio") .and. yb@ratio)then
        do y = yb,ye
            catstsa(:,{y}) = catstsa(:,{y}) / sum(catstsa(:,{y})) *100
        end do
    end if
    if(isatt(yb,"stddevscale"))then
        catstsa@stddevscale = yb@stddevscale
    end if
    if(isatt(yb,"ymax"))then
        catstsa@ymax = yb@ymax
    end if
    e = plot_time_series(catstsa,title,filename)

    return True
end

undef("plot_index_monclm")
function plot_index_monclm(years,iindex,title,filename)
begin
    if(isatt(years,"dataset"))then
        dataset = years@dataset
    else
        dataset = "jtwc"
    end if
    a  = read_besttrack_data(dataset)
    if (isatt(years,"emanuel2005Correcting").and.years@emanuel2005Correcting)then
        print("emanuel2005Correcting")
        a = emanuelCorrecting(a) 
    end if
    
    if(iindex .eq."day")then
        a@todo = "hr"
    else
        a@todo = iindex
    end if
    index = global_MonRvKE(a)
    
    ny = dimsizes(years)
    yyyymm = new(ny*12,"integer")
    do y = 0, ny-1
        yyyymm(y*12:y*12+11) = yyyymm_time(years(y),years(y),"integer")
    end do
    clmmon = clmMonTLL(index({yyyymm},:,:))
    ts = dim_sum_Wrap(dim_sum_Wrap(clmmon))
    if(iindex.eq."ke".or.iindex.eq."rvke")then
        ts = ts/10000
    end if
    if(iindex.eq."pdi".or.iindex.eq."rpdi")then
        ts = ts/1000000000
    end if
    if(iindex .eq."day")then
        ts = ts/24
    end if
    if(isatt(years,"norm").and.years@norm)then
        title@gsnRightString ="Factor="+sprintf("%4.2f",max(ts))
        ts = ts/max(ts)
    end if
        title@tiMainPosition = "Left"
        title@tiMainJust = "CenterLeft"
    ts&$ts!0$ = ispan(1,12,1)

    p      = plot_time_series(ts,title,filename)
    return ts
end

undef("plot_index_pentaclm")
function plot_index_pentaclm(years,iindex,title,filename)
begin
    todo = iindex
    thresholds = (/ 34., 83.,114./)
    limits     = (/ 82.,113.,999./)
    if(isatt(years,"RPDIpaper"))then
        thresholds = (/ 34., 64., 96./)
        limits     = (/ 63., 95.,999./)
    end if
    cats  = new((/3,73/),"float")
        cats@clm = True
    do i =0, 2
        todo@threshold = thresholds(i)
        todo@limit = limits(i)
        penta = read_jtwc_penta_rvke(todo)
        cats(i,:) = dim_sum_Wrap(dim_sum_Wrap(dim_avg_n_Wrap(penta({years},:,:,:),0)))
    end do
    if(isatt(years,"norm").and.years@norm)then
        title@gsnRightString ="Factor="+sprintf("%4.2f",max(cats))
        allty = cats(0,:)
        allty = allty+cats(1,:)
        allty = allty+cats(2,:)
        do i = 0,2
            cats(i,:) = cats(i,:)/max(allty)
        end do
    end if
    ii = iindex
    if (     ii.eq."ke".or.ii.eq."ace".or.ii.eq."ACE" \
        .or. ii.eq."rvke".or.ii.eq."race".or.ii.eq."RACE")then  ; for ACE and RACE
        cats = cats/10000
    end if
    if (ii.eq."pdi".or.ii.eq."PDI".or.ii.eq."rpdi".or.ii.eq."RPDI")then  ; PDI and RPDI
        cats= cats/1000000000
    end if
        title@tiMainPosition = "Left"
        title@tiMainJust = "CenterLeft"
    cats@rave = 3
    cats@nostddev = True
    cats@tiXAxisString  = "Juilan Pentad"
    p      = plot_time_series(cats,title,filename)
    a = print_cats_percent(cats)
    return cats
end

undef("plot_index_dots")
function plot_index_dots(iyears,iindex,title,filename)
begin
    iyears!0 = "years"
    iyears&years  = iyears
    if(isatt(iyears,"dataset"))then
        dataset = iyears@dataset
    else
        dataset = "jtwc"
    end if
    if(any(dataset.eq.(/"jma"/)))then
        miny = 1977
        maxy = 2008
    else
        miny = 1965
        maxy = 2008
    end if
    years = iyears({miny:maxy})
    a  = read_besttrack_data(dataset)
    if (isatt(years,"emanuel2005Correcting").and.years@emanuel2005Correcting.and.dataset.eq."jtwc")then
        print("emanuel2005Correcting")
        a = emanuelCorrecting(a) 
    end if

    a@todo = iindex
    if(isatt(iindex,"threshold"))
        a@threshold = iindex@threshold
    end if
    if(isatt(iindex,"limit"))
        a@limit = iindex@limit
    end if
    index = global_MonRvKE(a)  ; read all data
    times = index&time
    yyyy  = ispan(min(times/100),max(times/100),1)
    dims = dimsizes(index)
    if (isatt(years,"norm"))then
        norm = years@norm
    else
        norm = False
    end if
    if (isatt(years,"mon"))then
        mon = years@mon
        do m = 0,11
            if(.not.any(mon.eq.(m+1)))then
                index(m::12,:,:) = 0
            end if
        end do
    else
        mon = 0  ; per year
    end if
    if (isatt(years,"sigtest"))then
        sigtest = years@sigtest
    else
        sigtest = False
    end if
    yrindex = month_to_annual(index,0) ; sum
    yrindex&year = yyyy
    yrindex = yrindex/dimsizes(mon)    ; per month if mon.ne.0
    print(years)
    indexareas = yrindex(lat|:,lon|:,{year|years}); filter needed years
    if (dimsizes(dimsizes(indexareas)).gt.2)then
        indexarea = dim_avg_Wrap(indexareas) 
    else
        indexarea = indexareas
    end if
    delete(indexareas)

        ;; cal ratio of area/total
        ;total = sum(indexarea)
        ;area  = sum(indexarea({lat|10:20},{lon|110:160}))
        ;ra    = area/total
        ;print("ratio: "+ra)
    res =  True
    res@title = title
    res@filename = filename
    if(iindex.eq."pdi")then
        res@datamax = 0.
    end if
    if(iindex.eq."ke")then
        res@datamax = 0.
    end if
    if(iindex.eq."gen")then
        res@datamax = 0.2 
        if(any(mon.eq.0))then
            res@datamax = 2.0 
        end if  
    end if
    if(iindex.eq."hr")then
        res@datamax = 00.
    end if
    if(iindex.eq."rvke")then
        res@datamax = 0.4
        if(any(mon.eq.0))then
            res@datamax = 0.8 
        end if  
    end if

    if(isatt(years,"clm").and.years@clm)then
        res@datamax = 0.0
    else
        clmindex  = dim_avg_Wrap(yrindex(lat|:,lon|:,{year|miny:maxy}))
        ;print("sum index = "+sum(indexarea*dimsizes(mon)/10000))
        ;print("clm sum index = "+sum(clmindex*dimsizes(mon)/10000))
        indexarea = indexarea - clmindex
        if(sigtest)then
            siglvl    = ttest_data(yrindex({year|1965:2008},lat|:,lon|:),yrindex({year|years},lat|:,lon|:))
            indexarea = where(siglvl.le.0.1,indexarea,0.) ;; only significant
            indexarea = where(ismissing(indexarea),0,indexarea) ;; set missing to 0
        end if
    end if
    ii = iindex 
    if (     ii.eq."ke".or.ii.eq."ace".or.ii.eq."ACE" \
        .or. ii.eq."rvke".or.ii.eq."race".or.ii.eq."RACE")then  ; for ACE and RACE
        indexarea = indexarea /10000
    end if
    if (ii.eq."pdi".or.ii.eq."PDI".or.ii.eq."rpdi".or.ii.eq."RPDI")then  ; PDI and RPDI
        indexarea = indexarea /1000000000
    end if
    if(norm)then
        normfactor = max(abs(indexarea))
        indexarea = indexarea/normfactor
        indexarea@normfactor = normfactor
        res@datamax = 1.
        ;res@tiMainPosition = "Left"
        ;res@tiMainJust = "CenterLeft"
        res@gsnRightString = "Factor="+sprintf("%4.2f",normfactor)
    end if
    if(isatt(years,"overmark"))then
        res@overmark = years@overmark
    else
        res@overmark = max(abs(indexarea))/2
    end if
    print(ii+" sum: "+sum(indexarea))
    if(isatt(years,"rectangle"))then
        indexarea@rectangle = years@rectangle
    end if
    if(isatt(years,"monocolor"))then
        indexarea@monocolor = years@monocolor
    end if

    c = plot_dot(indexarea,res)
    return indexarea
end

undef("index_pattern_corr")
function index_pattern_corr(years,iindex1,iindex2)
begin
    i1 = plot_index_dots(years,iindex1,"tmp","tmp")
    i2 = plot_index_dots(years,iindex2,"tmp","tmp")
    i1@_FillValue = -999.
    i2@_FillValue = -999.
    cor = escorc(ndtooned(i1({lat|5:30},{lon|110:170})),ndtooned(i2({lat|5:30},{lon|110:170})))
    ;;cor = escorc(ndtooned(i1({lat|:},{lon|:})),ndtooned(i2({lat|:},{lon|:})))
    return cor
end

undef("plot_race_ace_ratio_dots")
function plot_race_ace_ratio_dots(years,title,filename)
begin
    if(isatt(years,"dataset"))then
        dataset = yb@dataset
    else
        dataset = "jtwc"
    end if
    a  = read_besttrack_data(dataset)

    a@todo = "rvke"
    race = global_MonRvKE(a)
    a@todo = "ke"
    ace = global_MonRvKE(a)
    times = race&time
    yyyy  = ispan(min(times/100),max(times/100),1)
    dims = dimsizes(race)
    yrrace = month_to_annual(race,0)
    yrace  = month_to_annual( ace,0)
    yrrace&year = yyyy
    yrace&year  = yyyy
    racearea = dim_avg_Wrap(yrrace(lat|:,lon|:,{year|years}))
    acearea  = dim_avg_Wrap( yrace(lat|:,lon|:,{year|years}))
    acearea = where(acearea.eq.0,acearea@_FillValue,acearea)
    ratio   = racearea  ; trick for coordinate
    ratio   = ratio/acearea
    ratio   = where(ismissing(ratio),0,ratio)

    res =  True
    res@title = title
    res@filename = filename
    ;res@datamax = 1.

    if(isatt(years,"clm").and.years@clm)then
        res@datamax = 0.8
        res@datamin = 0.
    else
        clmrace  = dim_avg_Wrap(yrrace(lat|:,lon|:,{year|1965:2008}))
        clmace   = dim_avg_Wrap( yrace(lat|:,lon|:,{year|1965:2008}))
        ratio    = ratio - (clmrace/clmace)
    end if
    c = plot_dot(ratio,res)

    ;; try "RasterFill" 
    cres = True 
    cres@cnFillMode = "RasterFill"
    cres@cnFillOn         = True  
    cres@cnLinesOn        =  False 
    cres@cnLevelSelectionMode = "ManualLevels"
    cres@cnMinLevelValF = 0.2
    cres@cnMaxLevelValF = 0.55
    cres@cnLevelSpacingF= 0.05
    cres@vpWidthF   = 0.80                          ; make map bigger
    cres@vpHeightF  = 0.80
    cres@mpMaxLatF  =  40.                           ; select subregion
    cres@mpMinLatF  =  00.
    cres@mpMinLonF  = 110.
    cres@mpMaxLonF  = 180.
    cres@mpCenterLonF = 180.
    cres@tiMainString = title

    wks = gsn_open_wks("ps",filename+"_m")          ; open a ps file
    gsn_define_colormap(wks,"temp_diff_18lev")
    ratio = where(ratio.lt.0.1,ratio@_FillValue,ratio)
    plot = gsn_csm_contour_map_ce(wks,ratio, cres) 
    return ratio
end

undef("plot_hrratio_RACEratio")
function plot_hrratio_RACEratio(years,mons,title,filename)
begin
    thresholds = (/ 34., 83.,114./)
    limits     = (/ 82.,113.,999./)
    if(isatt(yb,"dataset"))then
        dataset = yb@dataset
    else
        dataset = "jtwc"
    end if
    a  = read_besttrack_data(dataset)

    a@todo = "rvke"
    race = global_MonRvKE(a) 
    a@todo = "ke"
    ace  = global_MonRvKE(a) 
    a@todo = "hr"
    hrall= global_MonRvKE(a)  ; read all data
    a@threshold  = thresholds(2)
    a@limit      = limits(2)
    hr4  = global_MonRvKE(a)  ; read all data
    a@threshold  = thresholds(0)
    a@limit      = limits(0)
    hr01 = global_MonRvKE(a)  ; read all data
    a@threshold  = thresholds(1)
    a@limit      = limits(1)
    hr23 = global_MonRvKE(a)  ; read all data

    if (isatt(years,"mon"))then
        mon = years@mon
        do m = 0,11
            if(.not.any(mon.eq.(m+1)))then
                race(m::12,:,:) = 0
                ace (m::12,:,:) = 0
                hr01(m::12,:,:) = 0
                hr4 (m::12,:,:) = 0
            end if
        end do
    else
        mon = 0  ; per year
    end if
    ;; all to yearly
    yrrace = month_to_annual(dim_sum(dim_sum(race(time|:,lat|:,lon|:))),0)
    yrace  = month_to_annual(dim_sum(dim_sum(ace (time|:,lat|:,lon|:))),0)
    yrhrall= month_to_annual(dim_sum(dim_sum(hrall(time|:,lat|:,lon|:))),0)
    yrhr01 = month_to_annual(dim_sum(dim_sum(hr01(time|:,lat|:,lon|:))),0)
    yrhr23 = month_to_annual(dim_sum(dim_sum(hr23(time|:,lat|:,lon|:))),0)
    yrhr4  = month_to_annual(dim_sum(dim_sum(hr4 (time|:,lat|:,lon|:))),0)
    yrfac  = yrrace
    yrfac  = yrrace/yrace
    yrhrr01= yrhr01
    yrhrr01= yrhr01/yrhrall
    yrhrr23= yrhr23
    yrhrr23= yrhr23/yrhrall
    yrhrr4 = yrhr4 
    yrhrr4 = yrhr4 /yrhrall
    ;; set year coordinate
    times = race&time
    yyyy  = ispan(min(times/100),max(times/100),1)
    yrfac&year   = yyyy
    yrhrr01&year = yyyy
    yrhrr23&year = yyyy
    yrhrr4&year  = yyyy

    ;; plot scatter
    wks   = gsn_open_wks ("ps",filename)      
    res                   = True   
    res@gsnDraw  = False
    res@gsnFrame = False
    res@tiMainString = title
    res@xyMarkLineMode    = "Markers"                ; choose to use markers
    res@xyMarkerColor     = "blue"           ; Marker color
    res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
    plot = gsn_csm_xy(wks,yrhrr01({year|years}),yrfac({year|years}),res)
    draw(plot)
    frame(wks)
    res@xyMarkerColor     = "ForestGreen"            ; Marker color
    plot = gsn_csm_xy(wks,yrhrr23({year|years}),yrfac({year|years}),res)
    draw(plot)
    frame(wks)
    res@xyMarkerColor     = "red"            ; Marker color
    plot = gsn_csm_xy(wks,yrhrr4 ({year|years}),yrfac({year|years}),res)
    draw(plot)
    frame(wks)
    return True

end

undef("plot_pdi_dots")
function plot_pdi_dots(years,title,filename)
begin
    a = plot_index_dots(years,"pdi",title,filename)
    return a
end

undef("plot_pdi_ts")
function plot_pdi_ts(yb,ye,title,filename)
begin
    a = plot_index_ts(yb,ye,"pdi",title,filename)
    return a
end

undef("plot_line_ts")  ;; yearly
function plot_line_ts(yb,ye,todos,title,filename)
begin
    if(isatt(yb,"dataset"))then
        dataset = yb@dataset
    else
        dataset = "jtwc"
    end if
    a  = read_besttrack_data(dataset)

    if (isatt(yb,"emanuel2005Correcting").and.yb@emanuel2005Correcting)then
        print("emanuel2005Correcting")
        a = emanuelCorrecting(a) 
    end if
    ny = ye-yb+1
    ntodo = dimsizes(todos)
    tss   = new((/ntodo,ny/),"float")
    tss!1 = "year"
    tss&year = ispan(yb,ye,1)
    tss!0 = "var"
    tss&var = todos
    scales = new(ntodo,"float")

    ymb = yb*100+1
    yme = ye*100+12
    do i = 0, ntodo-1
        a@todo = todos(i)
        index = global_MonRvKE(a)
        indexts = dim_sum_Wrap(dim_sum_Wrap(index({ymb:yme},:,:)))
        indextsyr = month_to_annual(indexts,0)
        tss(i,:) = (/indextsyr/)
        delete(index)
        delete(indexts)
        delete(indextsyr)
    end do

    if (isatt(yb,"smooth121").and.yb@smooth121)then ;; from Emanuel2005
        tss = wgt_runave_Wrap(tss,(/1,2,1/),0)
        tss = wgt_runave_Wrap(tss,(/1,2,1/),0)
    end if

    do i = 0, ntodo-1 ;; scale trick: tss(:,{1974}) = 1.   
        scales(i) = 1./sum(tss(var|i,{year|1973}))
        tss(i,:) = tss(i,:)*scales(i)
    end do

    ;print("tss: "+tss(0,:)+" "+tss(1,:)+" "+tss(2,:)+" "+tss(3,:)+" ")
    res = True
    res@tiMainString = title
    res@trYMinF  = 0.7
    res@pmLegendDisplayMode    = "Always"
    res@pmLegendSide = "Top"
    res@lgPerimOn = False 
    res@pmLegendParallelPosF   = .32                 ; move units right
    res@pmLegendOrthogonalPosF = -0.40                ; move units down
    res@xyExplicitLegendLabels = tss&var+" * "+scales
    res@pmLegendWidthF         = 0.25                ; Change width and
    res@pmLegendHeightF        = 0.15                ; height of legend.
    res@lgLabelFontHeightF     = .015                 ; change font height

    res@xyDashPattern = 0               ; Make curves all solid
    res@xyMarkLineMode = "MarkLines"    ; Markers *and* lines
    res@xyMarkers      = (/6,11,16,18,0,7/)    ; 3 different markers
    res@xyMarkerColors = (/"blue","red","darkgreen","black","tomato4","darkorange4"/)  ; 3 different colors
    res@xyLineColors = res@xyMarkerColors

    wks = gsn_open_wks("ps",filename) 

    plot=gsn_csm_xy(wks,tss&year,tss,res)

    return True
end
