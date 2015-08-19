load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;; for drawNDCGrid()
    ;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    ;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
    ;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("plot_coloredScatter")
function plot_coloredScatter(x[*],y[*],c[*],title[1]:string,filename[1]:string)
begin
    res = True
    res@tiMainString      = title
    res@tiMainJust = "TopLeft"
    res@tiMainFont      = 21
    res@tiMainFontHeightF = 0.018
    res@tiMainPosition      = "Left"
    res@xyMarker  = 1
    res@xyMarkLineMode    = "Markers"  
    res@gsnFrame = False
    res@gsnDraw = False
    
    wks = gsn_open_wks("ps",filename)
    ;gsn_define_colormap(wks,"temp_19lev")
    gsn_define_colormap(wks,"so4_21")
    nc = 21 ; 
    print("plot "+filename)

    xmaxmin = nice_mnmxintvl(min(x),max(x),nc,True) ;; nc whatever
    ymaxmin = nice_mnmxintvl(min(y),max(y),nc,True)
    cmaxmin = nice_mnmxintvl(min(c),max(c),nc,True)
        ;; set x,y and color range if not setted
        xmin = xmaxmin(0)
        if(isatt(x,"min"))then
            xmin = x@min
        end if
        xmax = xmaxmin(1)
        if(isatt(x,"max"))then
            xmax = x@max
        end if
        ymin = ymaxmin(0)
        if(isatt(y,"min"))then
            ymin = y@min
        end if
        ymax = ymaxmin(1)
        if(isatt(y,"max"))then
            ymax = y@max
        end if
        cmin = cmaxmin(0)
        if(isatt(c,"min"))then
            cmin = c@min
        end if
        cmax = cmaxmin(1)
        if(isatt(c,"max"))then
            cmax = c@max
        end if
        cint = (cmax - cmin)/nc

    np = dimsizes(c)

    tres = res
    tres@xyMarkerColor= 0
    plot = gsn_csm_xy(wks,(/xmin,xmax/),(/ymin,ymax/),tres)   

    pres =True
    pres@gsMarkerIndex = 16 
    do i =0, np-1
        if(.not.(ismissing(x(i)).or.ismissing(y(i)).or.ismissing(c(i))))then
            pres@gsMarkerColor  = floattoint((c(i)-cmin)/cint)+1
            if(pres@gsMarkerColor .lt. 2)then
                pres@gsMarkerColor = 2
            end if
            if(pres@gsMarkerColor .gt. nc)then
                pres@gsMarkerColor = nc
            end if
            gsn_polymarker(wks,plot,x(i),y(i),pres)
        end if
    end do
    lbres                           = True
    lbres@lbAutoManage              = False         ; we control label bar
    lbres@lbFillColors              = ispan(2,nc,1); use nice strided colors  
    lbres@lbPerimOn                 = False       ; Turn off labelbar perimeter.
    lbres@lbMonoFillPattern         = True          ; one pattern, all solid
    lbres@vpWidthF             = 0.1                 ; size
    lbres@vpHeightF            = 0.6
    lbres@lbLabelFontHeightF   = 0.012               ; label font height

    lblabels = flt2string(fspan(cmin,cmax,nc))
    lbx = 0.8
    lby = 0.8
    gsn_labelbar_ndc(wks,nc,lblabels,lbx,lby,lbres) 
    ;;drawNDCGrid(wks)


    draw(plot)
    frame(wks)
    return True
end

function plot_colorScatter(title[1]:string,filename[1]:string,x[*],y[*],c[*],ires:logical)
begin
    if(all(ismissing(x)).or.all(ismissing(y)))then
        print(" colorScatter(): "+filename+" all missing value")
        return True
    end if
    res = True
    res@tiMainString      = title
    res@tiMainJust = "TopLeft"
    res@tiMainFont      = 21
    res@tiMainFontHeightF = 0.018
    res@tiMainPosition      = "Left"
    res@xyMarker  = 1
    res@xyMarkLineMode    = "Markers"  
    res@gsnFrame = False
    res@gsnDraw = False
    
    copy_VarAtts(ires,res)

    stdshading = True
    stdshading = False
    if(num(.not.ismissing(x)).lt.3)then
        stdshading = False
    end if
    if(stdshading)then
        reg = regline(x,y)
        reg@rstd = 1 ;; for test
        regx = (/min(x),max(x),max(x),min(x)/)
        regy = reg@yintercept+reg*regx
        regy = (/regy(0)-reg@rstd,regy(1)-reg@rstd,regy(2)+reg@rstd,regy(3)+reg@rstd/)

        ssres = True
            ssres@tfPolyDrawOrder = "PreDraw"
            ssres@gsFillColor = "SlateBlue"
        regres = True
            regres@txFontHeightF = 0.015    
            regres@txJust        = "CenterLeft"
    end if


    wks = gsn_open_wks("ps",filename)
    ;gsn_define_colormap(wks,"temp_19lev")
    gsn_define_colormap(wks,"so4_21")
    nc = 21 ; 
    print("plot "+filename)


    xmaxmin = nice_mnmxintvl(min(x),max(x),nc,True) ;; nc whatever
    ymaxmin = nice_mnmxintvl(min(y),max(y),nc,True)

    if(dimsizes(c).gt.2)then
        cmaxmin = nice_mnmxintvl(min(c),max(c),nc,True)
        monocolor = False
    else 
        cmaxmin = (/c,c,c/)
        monocolor = True
    end if
        ;; set x,y and color range if not setted
        xmin = xmaxmin(0)
        if(isatt(x,"min"))then
            xmin = x@min
        end if
        xmax = xmaxmin(1)
        if(isatt(x,"max"))then
            xmax = x@max
        end if
        ymin = ymaxmin(0)
        if(isatt(y,"min"))then
            ymin = y@min
        end if
        ymax = ymaxmin(1)
        if(isatt(y,"max"))then
            ymax = y@max
        end if
        cmin = cmaxmin(0)
        if(isatt(c,"min"))then
            cmin = c@min
        end if
        cmax = cmaxmin(1)
        if(isatt(c,"max"))then
            cmax = c@max
        end if
        cint = (cmax - cmin)/nc

    np = dimsizes(x)
    tres = res
    tres@xyMarkerColor= 0
    plot = gsn_csm_xy(wks,(/xmin,xmax/),(/ymin,ymax/),tres)   

    if(stdshading)then
        printVarSummary(regx)
        printVarSummary(regy)
        gsn_polygon(wks,plot,regx,regy,ssres)
        gsn_text_ndc(wks,"y = "+reg+"x + "+reg@yintercept,0.5,0.75,regres)
        gsn_text_ndc(wks,"rstd =  "+reg@rstd,0.5,0.73,regres)
    end if

    pres =True
    pres@gsMarkerIndex = 16 
    pres@tfPolyDrawOrder   = "Postdraw"
    if(monocolor)then
        pres@gsMarkerColor  = toint(c)
    end if
    do i =0, np-1
        if(.not.(ismissing(x(i)).or.ismissing(y(i))))then
            if((.not.monocolor).and.(ismissing(c(i))))then
                continue
            end if
            if(.not.monocolor)then
                pres@gsMarkerColor  = toint((c(i)-cmin)/cint+1)
                if(pres@gsMarkerColor .lt. 2)then
                    pres@gsMarkerColor = 2
                end if
                if(pres@gsMarkerColor .gt. nc)then
                    pres@gsMarkerColor = nc
                end if
            end if
            gsn_polymarker(wks,plot,x(i),y(i),pres)
        end if
    end do
    if(.not.monocolor)then
        lbres                           = True
            lbres@lbAutoManage              = False         ; we control label bar
            lbres@lbFillColors              = ispan(2,nc,1); use nice strided colors  
            lbres@lbPerimOn                 = False       ; Turn off labelbar perimeter.
            lbres@lbMonoFillPattern         = True          ; one pattern, all solid
            lbres@vpWidthF             = 0.1                 ; size
            lbres@vpHeightF            = 0.6
            lbres@lbLabelFontHeightF   = 0.012               ; label font height
        lblabels = flt2string(fspan(cmin,cmax,nc))
        lbx = 0.8
        lby = 0.8
        gsn_labelbar_ndc(wks,nc,lblabels,lbx,lby,lbres) 
    end if
    ;;drawNDCGrid(wks)


    draw(plot)
    frame(wks)
    if(isvar("reg"))then
        return reg
    end if
    return True
end
