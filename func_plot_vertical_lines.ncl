load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


function create_std_xz_polygon(avgline[*],stdline[*],z[*])
begin
    ;; create polygon coordinate for gsn_polygon
    ;; avgline +- stdline in x,z asix
    ;; return coords(2,nz*2) for x,y

    left = avgline-stdline
    right = avgline+stdline

    ;; 
    x = array_append_record(left,right(::-1),0)
    y = array_append_record(z,z(::-1),0)
    dims = dimsizes(x)

    coords = reshape(x,(/2,dims/))
    coords(0,:) = tofloat(x)
    coords(1,:) = tofloat(y)
    return coords
end

function plot_vertical_lines(title[1]:string,filename[1]:string,x[*][*]:numeric,z[*]:numeric,ires)
begin
    ;; assume x is (nline,lev)
    print(" plot vertical_lines: "+filename)
    xdims = dimsizes(x)
    nline = xdims(0)
    nlev =  xdims(1)
    res = ires
        res@gsnLeftString = title
        res@trYReverse        = True
    if(isatt(res,"nostd"))then
        nostd = res@nostd
        delete(res@nostd)
    else
        nostd = False
    end if

    wks = gsn_open_wks("ps",filename)

    if(nline.gt.100.and.(.not.nostd))then ;; too many lines, draw average line and std range
        avgline = dim_avg_n(x,0)
        stdline = dim_stddev_n(x,0)
        std1poly = create_std_xz_polygon(avgline,stdline,z)
        std2poly = create_std_xz_polygon(avgline,stdline*2,z)
            res@gsnFrame = False
            res@gsnDraw  = False
        res@xyLineColor  = "black"
        res@xyLineThicknessF = 3.0

        plot = gsn_csm_xy(wks,avgline,z,res)

        shadres = True
            shadres@gsFillColor = "cyan"
        gsn_polygon(wks,plot,std2poly(0,:),std2poly(1,:),shadres)
            shadres@gsFillColor = "cyan4"
        gsn_polygon(wks,plot,std1poly(0,:),std1poly(1,:),shadres)

        ;; add x = 0 ref line
        yrefres = True
            yrefres@gsLineColor = "red"
        x0 = z*0
        gsn_polyline(wks,plot,x0,z,yrefres)

        draw(plot)
        frame(wks)
        return True
    end if
    

    plot = gsn_csm_xy(wks,x,z,res)

    if(isatt(res,"gsnDraw").and.(.not.res@gsnDraw))then
        draw(plot)
    end if
    if(isatt(res,"gsnFrame").and.(.not.res@gsnDraw))then
        frame(wks)
    end if

    return True
end

