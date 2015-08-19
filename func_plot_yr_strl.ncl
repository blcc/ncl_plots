load "func_plot_yr_onset.ncl"
undef("plot_strl")
function plot_strl(u1[*][*],v1[*][*],r1[*][*],u2[*][*],v2[*][*],r2[*][*],ires)
begin
    ;; plot 3 stream line figs at 1 file
    ;; fig 1 and 2 is preonset and postonset(in monsoon season)
    ;; fig 3 is the diff between 1 and 2
    ;; assume vars are lat*lon
    if(isatt(ires,"filename"))then
        filename = ires@filename
    else
        filename = "test_strl"
    end if
    if(isatt(ires,"lev"))then
        slev = ires@lev+" hPa + Prec."
    else
        slev = " + Prec."
    end if
    if(isatt(ires,"yeartitle"))then
        sytitle = " "+ires@sytitle+" "
    else
        sytitle = " "
    end if
    anu = u2
    anu = u2-u1
    anv = v2
    anv = v2-v1
    anr = r2
    anr = r2-r1

    res1 = True
    res1@gsnDraw             = False
    res1@gsnFrame            = False
    res1@stArrowLengthF      = 0.01
    res1@lbLabelBarOn        = False            ; turn off individual label bars
    res1@mpFillOn            = False            ; no need   
    res1@mpMinLonF           =  60.
    res1@mpMaxLonF           = 180.
    res1@mpMinLatF           = -20.
    res1@mpMaxLatF           =  50.
    res1@mpGridLineDashPattern= 2
    res1@gsnMajorLatSpacing  =  10
    res1@cnLevelSelectionMode= "ManualLevels"   ; manual set levels
    res1@cnMinLevelValF      =  6.0
    res1@cnMaxLevelValF      = 16.0
    res1@cnLevelSpacingF     =  1.           
    res1@cnFillOn            = True             ; color fill plot
    res1@cnLinesOn           = False
    res1@cnLineLabelsOn      = False
    res1@cnInfoLabelOn       = False
    res1@tiMainString        = ""
    res1@tiMainFont      = 21
    res1@gsnStringFont       = 21
    res1@gsnRightString       = ""
    res1@gsnLeftString        = ""
    res1@gsnCenterString        = ""

    res2 = res1
    res3 = res1

    ;;;gsn_define_colormap(wks,"ViBlGrWhYeOrRe")
    res1@cnFillColors = (/54,63,66,70,73,76,80,83,86,90,95,102/)
    res2@cnFillColors = (/54,63,66,70,73,76,80,83,86,90,95,102/)
    ;;;gsn_define_colormap(wks,"BlWhRe")
    res3@cnFillColors = (/ 2,10,31,36,45,49,55,63,69,74,89,99/)

    res1@gsnLeftString     = "(A)"+sytitle+"WNP preonset - Streamline "+slev
    res2@gsnLeftString     = "(B)"+sytitle+"WNP postonset - Streamline "+slev
    res3@gsnLeftString     = "(C)"+sytitle+"WNP onset diff - Streamline "+slev
    res3@cnMinLevelValF      = -5.0
    res3@cnMaxLevelValF      =  5.0

    panres1                   = True                   ; modify the panel plot
    panres1@gsnFrame          = False                  ; don't advance frame yet
    panres1@gsnPanelLabelBar  = True                   ; add common colorbar
    panres1@lbLabelAutoStride = True
    panres1@lbOrientation     = "vertical"
    panres1@pmLabelBarWidthF  = 0.070                 ; make thinner
    panres1@lbTitleFont         = 21
    panres1@lbTitleFontHeightF  = 0.010
    panres1@lbLabelFont         = 21
    panres1@lbJustification     = "TopLeft"
    panres1@lbTitleOffsetF      = -0.05
    ;panres1@lbRightMarginF      = 0.01
    panres1@lbTitleString       = "[mm/day]"

    panres2 = panres1

    panres1@pmLabelBarHeightF = 0.64  
    panres2@pmLabelBarHeightF = 0.30  
    panres1@pmLabelBarParallelPosF = 0.01

    panres1@gsnPanelTop       = 0.9 
    panres1@gsnPanelBottom    = 0.3 
    panres2@gsnPanelTop       = 0.3
    ;panres2@gsnPanelBottom    =  .1


    wks  = gsn_open_wks("ps",filename) 
    gsn_define_colormap(wks,"ViBlGrWhYeOrRe")

    plotA = gsn_csm_streamline_contour_map_ce(wks,u1,v1,r1,res1) 
    plotB = gsn_csm_streamline_contour_map_ce(wks,u2,v2,r2,res2) 
    gsn_define_colormap(wks,"BlWhRe")
    plotC = gsn_csm_streamline_contour_map_ce(wks,anu,anv,anr,res3) 

    gsn_panel(wks,(/plotA,plotB/),(/2,1/),panres1)
    gsn_panel(wks,plotC,(/1,1/),panres2) 
    frame(wks)
    return True
end

undef("plot_monsoon_onset_strl")
function plot_monsoon_onset_strl(year,onsetpenta)
begin
    cachefile = "strl_cache_"+year+".nc"
    if(isfilepresent(cachefile))then
        cf = addfile(cachefile,"r")
        pentau = cf->pentau 
        pentav = cf->pentav
        pentar = cf->pentar
    else
        varname = "u"
        varname@odims = "TZYX"
        pentau = read_ncep_cmap_penta(year,varname)
        varname = "v"
        pentav = read_ncep_cmap_penta(year,varname)
        varname = "r"
        pentar = read_ncep_cmap_penta(year,varname)
        cf = addfile(cachefile,"c")
        cf->pentau = pentau
        cf->pentav = pentav
        cf->pentar = pentar
    end if
    onset = onsetpenta
    opt = True
    a = plot_strl(pentau(onset-1,3,:,:),pentav(onset-1,3,:,:),pentar(onset-1,:,:),pentau(onset,3,:,:),pentav(onset,3,:,:),pentar(onset,:,:),opt)

    return a
end

undef("plot_monsoon_onset_strl_clm")
function plot_monsoon_onset_strl_clm(years[*])
begin
    ny = dimsizes(years)
    if(isatt(years,"title"))then
        title = years@title
    end if
    cachefile = "strl_cache_"+"clm_"+title+".nc"
    if(isfilepresent(cachefile))then
        cf = addfile(cachefile,"r")
        pentau = cf->pentau 
        pentav = cf->pentav
        pentar = cf->pentar
    else
        print("reading.."+years(0))
        ycfile = "strl_cache_"+years(0)+".nc"
        if(isfilepresent(ycfile))then
            cf = addfile(ycfile,"r")
            pentau = cf->pentau 
            pentav = cf->pentav
            pentar = cf->pentar
        else
            varname = "u"
            varname@odims = "TZYX"
            pentau = read_ncep_cmap_penta(years(0),varname)
            varname = "v"
            pentav = read_ncep_cmap_penta(years(0),varname)
            varname = "r"
            pentar = read_ncep_cmap_penta(years(0),varname)
            cf = addfile(ycfile,"c")
            cf->pentau = pentau
            cf->pentav = pentav
            cf->pentar = pentar
        end if
        do iy = 1 , ny-1
            print("reading.."+years(iy))
            ycfile = "strl_cache_"+years(0)+".nc"
            if(isfilepresent(ycfile))then
                cf = addfile(ycfile,"r")
                pentau = pentau + cf->pentau 
                pentav = pentav + cf->pentav
                pentar = pentar + cf->pentar
            else
                varname = "u"
                pentau = pentau + read_ncep_cmap_penta(years(iy),varname)
                varname = "v"
                pentav = pentav + read_ncep_cmap_penta(years(iy),varname)
                varname = "r"
                pentar = pentar + read_ncep_cmap_penta(years(iy),varname)
                cf = addfile(ycfile,"c")
                cf->pentau = pentau
                cf->pentav = pentav
                cf->pentar = pentar
            end if
        end do
        pentau = pentau/ny
        pentav = pentav/ny
        pentar = pentar/ny
        cf = addfile(cachefile,"c")
        cf->pentau = pentau
        cf->pentav = pentav
        cf->pentar = pentar
    end if
    ;; for area average
    pi = 4.*atan(1.0)
    rad = (pi/180.)
    coswgt = cos(rad*pentau&lat)
    coswgt!0 = "lat"
    coswgt&lat = pentau&lat

    ;;areau = dim_avg_Wrap(dim_avg_Wrap(pentau(:,:,{5:20},{120:150})))
    ;;areav = dim_avg_Wrap(dim_avg_Wrap(pentav(:,:,{5:20},{120:150})))
    ;;arear = dim_avg_Wrap(dim_avg_Wrap(pentar(:,{5:20},{120:150})))
    areau = wgt_areaave(pentau(:,:,{5:20},{120:150}),coswgt({5:20}),1.,0)
    areav = wgt_areaave(pentav(:,:,{5:20},{120:150}),coswgt({5:20}),1.,0)
    arear = wgt_areaave(pentar(:,{5:20},{120:150}),coswgt({5:20}),1.,0)
    onoffset = wnponset(arear,areau,areav)
    onset = onoffset(0)
    print("onset = "+onset)
    opt = True
    a = plot_strl(pentau(onset-1,3,:,:),pentav(onset-1,3,:,:),pentar(onset-1,:,:),pentau(onset,3,:,:),pentav(onset,3,:,:),pentar(onset,:,:),opt)

    return a
end
