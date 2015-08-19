load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plot_addavgmark")
function plot_addavgmark(wks,plot,dd[*][*],res)
begin
    load "res_parameter.ncl"
    m = new((/5/),"float") ; WNP,IO,SH,EP,ATL
    m@title = (/ "WNP",  "IO",  "SH", "EP", "ATL"/)
    m@lon   = (/  170.,   80.,  120., 220.,  300./)
    m@lat   = (/   20.,   30.,  -25.,  30.,   30./)
    m(WNP)  = sum(dd({lat|0:60},{lon|100:210}))
    m(IO)   = sum(dd({lat|0:60},{lon| 30:100}))
    m(SH)   = sum(dd({lat|-60:0},{lon|:}))
    m(EP)   = sum(dd({lat|0:60},{lon|210:290}))
    m(ATL)  = sum(dd({lat|0:60},{lon|290:350}))
    tres       =  True
    tres@txFontHeightF = 0.01
    plot2 = gsn_add_text(wks,plot,m@title+" "+sprintf("%4.2f",m),m@lon,m@lat,tres)
    return plot2
end

undef("plot_global")
function plot_global(dd[*][*],title)
begin
    if (.not. isatt(dd,"avgmark"))then
        dd@avgmark = False
    end if
    res = True
    res@mpMinLatF = -40.
    res@mpMaxLatF =  40.
    res@mpMinLonF =  30.
    res@mpMaxLonF = 350.
    res@mpCenterLonF = 180.
    res@cnFillOn = True
    res@cnLinesOn= True 
    res@gsnSpreadColors = True 
    res@gsnSpreadColorStart =  12
    res@tiMainString = title

    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMaxLevelValF =     4
    res@cnMinLevelValF =  0000
    res@cnLevelSpacingF=     0.8

    res@gsnFrame = False
    res@gsnDraw  = False
    dd = where(dd.eq.0,dd-0.001,dd)

    wks  = gsn_open_wks("ps","fig")
    gsn_define_colormap(wks,"nrl_sirkes")
    plot = gsn_csm_contour_map_ce(wks,dd,res)
    if (dd@avgmark)then
        plot2 = plot_addavgmark(wks,plot,dd,res)
    end if
    draw(plot)
    frame(wks)
    return True
end
