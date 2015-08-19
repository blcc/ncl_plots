load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "/walker1/pgchiu/func_read/func_read_ersst.ncl"
;load "func_read_era_all.ncl"
load "func_read_reanalysis.ncl"

undef("read_wind_shear")
function read_wind_shear(years[*],mons[*])
begin
    vn = "U"
    vn@lev = 850
    u850 = read_reanalysis_monly(years,vn)
    vn = "V"
    v850 = read_reanalysis_monly(years,vn)
    vn = "U"
    vn@lev = 200
    u200 = read_reanalysis_monly(years,vn)
    vn = "V"
    v200 = read_reanalysis_monly(years,vn)
    ;; cal vertical wind shear in monthly data
    ushear = u200 - u850
    vshear = v200 - v850
    shearmons = u850
    shearmons = sqrt(ushear*ushear + vshear*vshear)
    shearclm  = clmMonTLL(shearmons)
    shearneed = dim_avg_n_Wrap(shearclm({mons-1},:,:),0)
    return shearneed
end
undef("plot_sst_vshear")
function plot_sst_vshear(years[*],mons[*],filename,title)
begin
    load "res_years.ncl"
    nm = dimsizes(mons)
    ;; read sst
    sst = read_ersst(years)
    sstclm = clmMonTLL(sst)
    sstneed = dim_avg_n_Wrap(sstclm({mons-1},:,:),0) ;; mon = 0-11
    delete(sst)
    delete(sstclm)
    if (isatt(years,"clm").and.years@clm)then
        sstspacing = 0.5
        sstmax     = 30.
        sstmin     = 26
        sstcolorstart= 11
    else
        sst = read_ersst(allyears)
        sstclm = clmMonTLL(sst)
        sstneed = sstneed - dim_avg_n_Wrap(sstclm({mons-1},:,:),0) ;; mon = 0-11
        sstspacing = 0.2
        sstmax     =  0.6
        sstmin     = -0.6
        sstcolorstart= 1
    end if

    ;; read wind to vertical wind shear(200hPa-850hPa)
    shearneed = read_wind_shear(years,mons)
    if (True.or.isatt(years,"clm").and.years@clm)then
        vsspacing = 2.
    else
        shearneed  = shearneed - read_wind_shear(allyears,mons)
        vsspacing = 1.
    end if

    ;; ploting
    wks = gsn_open_wks("ps",filename)

    gsn_define_colormap(wks,"temp_19lev") 

    res = True
    res@tiMainString = title
    res@tiMainJust = "CenterLeft"
    res@tiMainPosition      = "Left"
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@gsnDraw      = False
    res@gsnFrame     = False
    res@mpCenterLonF         = 180.               ; center plot at 180

    res@mpMinLonF            = 110.               ; select a subregion
    res@mpMaxLonF            = 180.
    res@mpMinLatF            =  00. 
    res@mpMaxLatF            =  40.
    res@cnLevelSelectionMode = "ManualLevels" 
    res@cnMinLevelValF = sstmin
    res@cnMaxLevelValF = sstmax
    res@cnLevelSpacingF  = sstspacing
    ;res@lbLabelStride  =  5
    res@gsnSpreadColorStart = sstcolorstart
    res@gsnSpreadColors = True


    res2 = True
    ;res2@cnLevelSelectionMode  = "ExplicitLevels"
    ;res2@cnLevels = (/27,30/)
    res2@cnLevelSpacingF  = vsspacing
    res2@cnLineLabelInterval = 1
    res2@cnInfoLabelOn = False
    res2@cnMonoLineThickness  = False
    res2@gsnContourNegLineDashPattern = 1
    res2@gsnContourZeroLineThicknessF = 3.
    ;res2@cnLineThicknesses =(/2.,1.5,1./)

    plot = gsn_csm_contour_map_overlay(wks,sstneed,shearneed,res,res2)
    draw(plot)
    frame(wks)
    print("plot: "+filename)

    return True
end

undef("plot_sst_vshear_SONdJJA")
function plot_sst_vshear_SONdJJA(years[*],filename,title)
begin
    load "res_years.ncl"
    ;; read sst
    sst = read_ersst(years)
    sstclm = clmMonTLL(sst)
    sstJJA = dim_avg_n_Wrap(sstclm({JJA-1},:,:),0) ;; mon = 0-11
    sstSON = dim_avg_n_Wrap(sstclm({SON-1},:,:),0) ;; mon = 0-11
    sstneed = sstSON
    sstneed = sstSON - sstJJA
    delete(sst)
    delete(sstclm)
    delete(sstJJA)
    delete(sstSON)

    sstspacing = 0.4
    sstmax     =  1.2
    sstmin     = -1.2
    sstcolorstart= 1

    ;; read wind to vertical wind shear(200hPa-850hPa)
    shearJJA  = read_wind_shear(years,JJA)
    shearSON  = read_wind_shear(years,SON)
    shearneed = shearSON
    shearneed = shearSON - shearJJA

    vsspacing = 1.

    ;; ploting
    wks = gsn_open_wks("ps",filename)

    gsn_define_colormap(wks,"temp_19lev") 

    res = True
    res@tiMainString = title
    res@tiMainJust = "CenterLeft"
    res@tiMainPosition      = "Left"
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@gsnDraw      = False
    res@gsnFrame     = False
    res@mpCenterLonF         = 180.               ; center plot at 180

    res@mpMinLonF            = 110.               ; select a subregion
    res@mpMaxLonF            = 180.
    res@mpMinLatF            =  00. 
    res@mpMaxLatF            =  40.
    res@cnLevelSelectionMode = "ManualLevels" 
    res@cnMinLevelValF = sstmin
    res@cnMaxLevelValF = sstmax
    res@cnLevelSpacingF  = sstspacing
    ;res@lbLabelStride  =  5
    res@gsnSpreadColorStart = sstcolorstart
    res@gsnSpreadColors = True


    res2 = True
    ;res2@cnLevelSelectionMode  = "ExplicitLevels"
    ;res2@cnLevels = (/27,30/)
    res2@cnLevelSpacingF  = vsspacing
    res2@cnLineLabelInterval = 1
    res2@cnInfoLabelOn = False
    res2@cnMonoLineThickness  = False
    res2@gsnContourNegLineDashPattern = 1
    res2@gsnContourZeroLineThicknessF = 3.
    ;res2@cnLineThicknesses =(/2.,1.5,1./)

    plot = gsn_csm_contour_map_overlay(wks,sstneed,shearneed,res,res2)
    draw(plot)
    frame(wks)
    print("plot: "+filename)

    return True
end
