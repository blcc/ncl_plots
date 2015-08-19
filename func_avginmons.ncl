
undef("avginmons")
function avginmons(var,mons)
begin
    dims = dimsizes(var)
    if(dims(0).ne.12)then
        print("avginmons() should use var[12][*][*] or var[12][*][*][*]")
        exit
    end if
    nm = dimsizes(mons)
    ndims = dimsizes(dimsizes(var))
    if(min(var&$var!0$).eq.0)then
        var&month = ispan(1,12,1)
    end if
    if(nm.gt.1)then
        if(ndims .eq. 3)then
            ovar = dim_avg_n_Wrap(var(mons-1,:,:),0)
        end if
        if(ndims .eq. 4)then
            ovar = dim_avg_n_Wrap(var(mons-1,:,:,:),0)
        end if
    else
        if(ndims .eq. 3)then
            ovar = var(mons-1,:,:)
        end if
        if(ndims .eq. 4)then
            ovar = var(mons-1,:,:,:)
        end if
    end if

    return ovar
end
undef("annual_sum")
function annual_sum(monts[*],mon[*])
begin
    if(.not.any(mon.eq.0))then
        yrts = monts(::12)
        yrts = 0
        do m = 0,11
            mm = m+1
            if(any(mm.eq.mon))then
                yrts = yrts + monts(m::12)
            end if
        end do
        nm = dimsizes(mon)
    else
        yrts = month_to_annual(monts,0)
        nm = 12
    end if
    yrts@nmonth = nm
    return yrts
end
undef("annual_sum3d")  ;; t,y,x
function annual_sum3d(monts[*][*][*],mon[*])
begin
    if(.not.any(mon.eq.0))then
        yrts = monts(::12,:,:)
        yrts = 0
        dims = dimsizes(yrts)
        do m = 0,11
            mm = m+1
            if(.not.any(mm.eq.mon))then
                monts(m::12,:,:) = monts@_FillValue
            end if
        end do
        do y = 0, dims(0)-1
            yrts(y,:,:) = dim_avg_n(monts(y*12:y*12+11,:,:),0)
        end do
        nm = dimsizes(mon)
        yrts = yrts*nm
    else
        yrts = month_to_annual(ts,0)
        nm = 12
    end if
    yrts@nmonth = nm
    return yrts
end
