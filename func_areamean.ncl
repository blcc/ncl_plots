
function areamean(x)
begin
    ;; areamean for fix lat/lon grid
    dims = dimsizes(x)
    ndim = dimsizes(dims)
    lat = x&$x!(ndim-2)$
    if(dims(ndim-2).ne.dimsizes(lat))then
        print("areamean(): lat size mis match")
        exit
    end if

    ;; from http://www.ncl.ucar.edu/Document/Functions/Built-in/wgt_areaave.shtml
    rad    = 4.0*atan(1.0)/180.0
    clat = cos(lat*rad)
    xavg = wgt_areaave(x,clat,1.0,0)
    do i = 0,ndim-3
        xavg!i = x!i
        xavg&$xavg!i$ = x&$x!i$
    end do
    return xavg
end

