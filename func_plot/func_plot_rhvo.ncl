load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;load "func_read_era_all.ncl"
load "func_read_reanalysis.ncl"

function plot_rhvo(years[*],mons[*],title,filename)
begin
    if(isatt(years,"RHlev"))then
        RHlev = years@RHlev
    else
        RHlev = 700
    end if
    if(isatt(years,"VORlev"))then
        VORlev = years@VORlev
    else
        VORlev = 850
    end if
    ;; read rh
    vn = "RH"
    vn@lev = RHlev
    rh = read_reanalysis_monclm(years,vn)
    prh = dim_avg_n_Wrap(rh({mons},:,:),0)
    ;; read vo and turn to absolute vorticity
    vn = "VOR"
    vn@lev = VORlev
    vor = read_reanalysis_monclm(years,vn)
    pvor = dim_avg_n_Wrap(vor({mons},:,:),0)
    if (isatt(years,"absvor"))then ;; to abs. vor.
        f = pvor
        dims = dimsizes(f)
        if(isatt(years,"clm").and.(.not.years@clm))then
            f = 0
        else
            do j = 0, dims(0)-1
                f(j,:) = (2.*2*3.14/86400.) * sin(f&lat(j)/180*3.14)
            end do
        end if
        pvor = pvor+f
    end if

    scale = 1.e05
    pvor = pvor*scale

    ;; set plot res
    res = True
    res@tiMainString = title
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
    res@cnMaxLevelValF    =  2. 
    res@cnMinLevelValF    = -2.
    res@cnLevelSpacingF   = .5
    res@gsnSpreadColors = True
    res@cnSmoothingOn = True
    res@cnSmoothingDistanceF = 0.005

    res2 = True
    res2@cnLevelSelectionMode = "ManualLevels"
    res2@cnMaxLevelValF    = 100. 
    res2@cnMinLevelValF    =   0.
    res2@cnLevelSpacingF   =   5.

    res2@cnInfoLabelOn = False

      ;; anomaly
        if(isatt(years,"clm").and.(.not.years@clm))then
            res2@cnLevelSpacingF   =   1.
            res@cnMaxLevelValF    =  0.5 
            res@cnMinLevelValF    = -0.5
            res@cnLevelSpacingF   =  0.1
        end if
    ;; plot
    wks = gsn_open_wks("ps",filename) 
    gsn_define_colormap(wks,"BlWhRe")

    pvor@long_name= ""
    pvor@units    = ""
    prh@long_name = ""
    prh@units    = ""
    plot = gsn_csm_contour_map_overlay(wks,pvor,prh,res,res2)
    res2@cnLineThicknessF   = 3.5
    res2@cnLevelSelectionMode = "ExplicitLevels"
    res2@cnLevels = 60.
    res2@gsnFrame = False
    res2@gsnDraw  = False
    res2@lbLabelBarOn = False
    res2@cnInfoLabelOn =False
    res2@cnLineLabelsOn = False
    plot2 = gsn_csm_contour(wks,prh,res2)
    overlay(plot,plot2)
    draw(plot)
    frame(wks)
    print("plot: "+filename)

    return True
end
