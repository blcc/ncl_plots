function setting_plot(varname,var,level,plottags[*]:string)
begin
    ;; setting plot resource including to var, level and plot type
    ;; return res includ color map: res@COLORMAP
    ;; level = 1000, 500, 200, "sfc" etc.
    ;; plot tags =  "contour", "line", 
    ;;              "timelat", "lonz", "mapce", "mapploar",
    ;;              "climate", "diff", "anomaly"

    res = True

    ;;; add conditions below
    ;;; add conditions above

    if (varname .eq. "")then ;; default automatic
        if(min(var)*max(var) .gt.0)then
            mnmxint = nice_mnmxintvl(min(var),max(var),21,True)
        else
            symMinMaxPlt(var,21,True,res)
        end if
        if (any(plottype.eq."contour"))then ;; for contour default
            res@COLORTABLE = 
            
        end if
    end if
end
