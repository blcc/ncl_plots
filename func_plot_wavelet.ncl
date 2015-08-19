load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

function plot_wavelet(wave)
local jtot, n
begin
    periodList    = (/1,2,4,8,16,32,64,128,256,512,1024,2048,4096/)
    jtot                          = dimsizes(wave@gws)
    n                             = dimsizes(wave@coi)
    power                         = onedtond(wave@power,(/jtot,n/))
    power!0                       = "period"
    power&period                  = wave@period
    power!1                       = "time"
    power&time                    = fspan(0,n-1,n)

    power@long_name               = "Power Spectrum"
    power@units                   = "1/unit-freq"

    SIG                           = power
    SIG                           = power/conform(power,wave@signif,0)
    SIG@long_name                 = "Significance"
    SIG@unit                      = " "

    pres                          = True
    pres@gsnDraw                  = False
    pres@gsnFrame                 = False
    pres@cnFillOn                 = True
    pres@cnFillMode               = "RasterFill"
    pres@cnRasterSmoothingOn      = True
    pres@cnLinesOn                = False
    pres@cnLineLabelsOn           = False
    pres@cnInfoLabelOn            = False
    pres@gsnSpreadColors          = True
    pres@trYReverse               = True
    pres@tmYLMode                 = "Explicit"
    pres@tmYLValues               = periodList
    pres@tmYLLabels               = flt2string(periodList)
    pres@tmLabelAutoStride        = True
    pres@vpHeightF                = .4
    pres@vpWidthF                 = .7

    pres@gsnRightString           = "Wavelet Power"
    pres@tiYAxisString            = "Period"

    resS                          = True
    resS@trYReverse               = True
    resS@tmYLMode                 = "Explicit"
    resS@tmYLValues               = periodList
    resS@tmYLLabels               = flt2string(periodList)
    resS@gsnDraw                  = False
    resS@gsnFrame                 = False
    resS@cnInfoLabelOn            = False
    resS@cnLinesOn                = False
    resS@cnLineLabelsOn           = False
    resS@cnFillScaleF             = 0.5
    resS@gsnLeftString            = " "
    resS@gsnRightString           = " "

    wks =  gsn_open_wks("ps","wavelet")
    gsn_define_colormap(wks,"BlAqGrYeOrReVi200")    ; choose colormap

    plot    = new(2,graphic)
    plot(0) = gsn_csm_contour(wks,power,pres)
    iplot   = gsn_csm_contour(wks,SIG,resS)
    iplot   = ShadeGtContour(iplot,1.,17)
    overlay(plot(0),iplot)

    gws                   = wave@gws
    resL                  = True
    resL@gsnFrame         = False
    resL@gsnDraw          = False
    resL@trYAxisType      = "LogAxis"
    resL@tmYLValues       = periodList
    resL@tmYLLabels       = flt2string(periodList)

    plotLeft = gsn_csm_xy(wks,gws,power&period,resL)

    draw(plot(0))
    frame(wks)

    return(True)
end

