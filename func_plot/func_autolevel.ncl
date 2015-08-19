;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("roughly_value")
function roughly_value(ovalue)
begin
    valueChar = stringtochar(flt2string(ovalue))
    pointLoc  = ind(valueChar .eq. ".")
    if (ismissing(pointLoc)) then
        valOrder  = dimsizes(valueChar) - 1
    else
        valOrder  = dimsizes(valueChar(:pointLoc-1)) - 1
    end if
    prec = 2
    rvalue  = stringtofloat(chartostring(valueChar(0)))*(10^(valOrder-1))   \
             +stringtofloat(chartostring(valueChar(1)))*(10^(valOrder-2))

    return rvalue
end

undef ("autolevel")
function autolevel(array,ncolors,opt)
;; opt = 0  output all(max,min,spacing) value in a array
;;       1  output max value
;;       2  output min
;;       3  output spacing
;; ncolor   number of colors you wish in figure
begin
    ncolor = ncolors
    true_max = max(array)
    true_min = min(array)
    true_spc = (true_max - true_min) / ncolor

    rmax     = roughly_value(true_max)
;    rmin     = roughly_value(true_min)
    rmin     = roughly_value(0.8*rmax)
;    rspc     = roughly_value(true_spc)
    rspc     = (rmax - rmin) / ncolor
    print("rmax/rmin/rspc: "+rmax+" "+rmin+" "+rspc)

    if(opt .eq. 0)then
        return (/rmax,rmin,rspc/)
    end if
    if(opt .eq. 1)then
        return rmax
    end if
    if(opt .eq. 2)then
        return rmin
    end if
    if(opt .eq. 3)then
        return rspc
    end if

end

