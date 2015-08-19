load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("rvke")
function rvke(vknot,todo)  ; from Yu et al. 2009, Revised ACE index
begin
    racename = (/"RACE","race","rvke","RvKE"/)
    acename  = (/"ACE" ,"ace","TyKE","tyke","ke"/)
    tyhrname = (/"TyHr","hr"  ,"Hr","TLS","TCD","tcd"/)
    tygenname= (/"TyGen","Tygen","gen","TGN"/)
    pdiname  = (/"pdi","PDI"/)

    if (any(todo.eq.racename))then
        v   = vknot
        vc  = 35.
        a   = 0.51
        rc  = (v/vc)^(1/a)
        ke  = (v*v/(rc*rc)) * (0.5 + ((rc^(2-(2*a))-1)/(1-a)))
        return ke
    end if
    if (todo .eq. "rvketail")then
        v   = vknot
        vc  = 35.
        a   = 0.51
        rc  = (v/vc)^(1/a)
        ke  = (1/(rc*rc)) * (0.5 + ((rc^(2-(2*a))-1)/(1-a)))
        return ke
    end if
    if (any(todo.eq.acename))then
        ke  = vknot*vknot
        return ke
    end if
    if (any(todo.eq.pdiname))then
        ;ke  = vknot*vknot*vknot
        v = vknot*0.51
        ke  = v*v*v*6*3600  ;; units m**3 s**-2
        return ke
    end if
    if (todo.eq."rpdi" .or. todo.eq."RPDI")then
        ;;v   = vknot
        v   = vknot*0.51
        vc  = 35.*0.51   ;; m/s
        a   = 0.51
        rc  = (v/vc)^(1/a)
        ke  = (v*v*v/(rc*rc)) * (0.4 + ((2*rc^(2-(3*a))-2)/(2-3*a)))
        ke  = ke*6*3600  ;; followed PDI
        return ke
    end if
    if (todo.eq."rpditail")then
        v   = vknot
        vc  = 35.
        a   = 0.51
        rc  = (v/vc)^(1/a)
        ke  = (1/(rc*rc)) * (0.4 + ((2*rc^(2-(3*a))-2)/(2-3*a)))
        return ke
    end if
    if (any(todo.eq.tyhrname))then
        ke  = 6.
        return ke
    end if
    if (todo.eq."day".or. todo.eq."TLSday")then
        ke  = 0.25
        return ke
    end if
    print("error in todo: "+todo)
    exit
end

undef("global_ix")
function global_ix(lon,reso)
begin
    if(ismissing(lon))then
        return lon
    end if
    nx = floattoint(360./reso )
    ix = floattoint((lon + (0.5*reso))/reso)
    ix = floattoint(lon/reso)
    if ( ix.eq. nx)then
        ix = 0
    end if
    return ix
end
undef("global_iy")
function global_iy(lat,reso)
begin
    if(ismissing(lat))then
        return lat
    end if
    ny = floattoint(180./reso +1)
    iy = floattoint((lat+90)/reso)
    if ( iy.ge. ny)then
        print("LAT wrong: "+lat)
    end if
    return iy
end

undef("isnext6hr")
function isnext6hr(y1,m1,d1,h1,y2,m2,d2,h2)
begin
    time1 = ut_inv_calendar( y1, m1, d1, h1, 00, 00,  "hours since 1800-01-01 00:00", 0 )
    time2 = ut_inv_calendar( y2, m2, d2, h2, 00, 00,  "hours since 1800-01-01 00:00", 0 )
   ; print(y1+"/"+m1+"/"+d1+"/"+h1)
   ; print(y2+"/"+m2+"/"+d2+"/"+h2)
   ; print("t1t2: "+time1+" "+time2)
    if ((time2-time1).eq.6)then
        return True
    end if
    return False
end

undef("global_MonGen")
function global_MonGen(tyrec[*][*])  ; JTWC only (TYNUM)
begin
    load "res_parameter.ncl"
    reso = 5.

    if(.not. isatt(tyrec,"threshold"))then  ; for max Vmax
        threshold = 33.         
    else
        print("threshold "+tyrec@threshold)
        threshold = tyrec@threshold
    end if
    if(.not. isatt(tyrec,"limit"))then
        limit = 999.
    else
        print("limit "+tyrec@limit)
        limit = tyrec@limit
    end if

    yb = floattoint(min(tyrec(:,YYYY)))
    ye = floattoint(max(tyrec(:,YYYY)))
    yyyymm = yyyymm_time(yb, ye, "integer")
    dims = dimsizes(tyrec)
    nx = floattoint(360./reso )
    ny = floattoint(180./reso )
    nt = floattoint((ye-yb+1)*12)
    lats = latGlobeFo(ny,"lats","latitude","degrees_north")
    lons = lonGlobeFo(nx,"lons","longitude","degrees_east")
    mongen = new((/nt,ny,nx/),"float")
    mongen!0 = "time"
    mongen!1 = "lat"
    mongen!2 = "lon"
    mongen&time = yyyymm
    mongen&lat = lats
    mongen&lon = lons

    mongen = 0.
    t1 = new(4,integer) ; yyyy,mm,dd,hh
    t2 = new(4,integer)
    t1 = 0
    t2 = 0
    oldloc = (/0.,0./)
    unrecty = True 
    oldtynum = 0
    maxlastdata = max(ind(.not.ismissing(tyrec(:,TYNUM))))
    do i = 0, dims(0)-1
        if (.not. any(ismissing(tyrec(i,1:))))then
            ix    = global_ix(tyrec(i,LON),reso)
            iy    = global_iy(tyrec(i,LAT),reso)
            it    = floattoint((tyrec(i,YYYY)-yb)*12  + tyrec(i,MM)) -1
            tynum = floattoint((tyrec(i,TYNUM)))
            loc   = (/tyrec(i,LAT),tyrec(i,LON)/)
            vmax  = tyrec(i,VMAX)
            t2(0) = floattoint(tyrec(i,YYYY))
            t2(1) = floattoint(tyrec(i,MM  ))
            t2(2) = floattoint(tyrec(i,DD  ))
            t2(3) = floattoint(tyrec(i,HH  ))
              ;if(.not.isnext6hr(t1(0),t1(1),t1(2),t1(3),t2(0),t2(1),t2(2),t2(3)) \
              ;.or. .not.all(loc.eq.oldloc))then
            if((oldtynum.ne.tynum))then ;  .or. t1(0).ne.t2(0)))then
                unrecty = True
                oldtynum = tynum
            end if
            if(.not.ismissing(vmax) .and. unrecty .and. vmax.ge.34 )then
                unrecty = False
                i1 = i   ;; for typhoon start and end position (i)
                if(     tyrec(i,TYNUM).eq.tyrec(maxlastdata-1,TYNUM)   \
                  .and. tyrec(i,YYYY) .eq.tyrec(maxlastdata-1,YYYY)   )then
                    i2 = dims(0)-1
                else 
                    i2 = i1+ min(ind(.not.(tyrec(i:,TYNUM).eq.tyrec(i,TYNUM))))
                end if
                ;;print("i1 "+i1+":"+tyrec(i1,TYNUM)+":"+tyrec(i1,YYYY))
                ;;print("i2 "+i2+":"+tyrec(i2,TYNUM)+":"+tyrec(i2,YYYY))
                maxvmax = max(tyrec(i1:i2,VMAX))
                if(maxvmax.ge.threshold .and. maxvmax.le.limit)then
                    mongen(it,iy,ix) = mongen(it,iy,ix)+1
                    ;if(tyrec(i,LAT).ge.30.)then
                    ;    print(tyrec(i,:))
                    ;end if
                end if
            end if
            t1 = t2
            oldloc = loc
        end if
    end do
    return mongen
end

undef("global_MonRvKE")
function global_MonRvKE(tyrec[*][*])
begin
    load "res_parameter.ncl"
    reso = 5.
    yb = floattoint(min(tyrec(:,YYYY)))
    ye = floattoint(max(tyrec(:,YYYY)))
    yyyymm = yyyymm_time(yb, ye, "integer")
    dims = dimsizes(tyrec)
    if(.not. isatt(tyrec,"todo"))then
        tyrec@todo = "rvke"
    else
        print("cal "+tyrec@todo)
    end if
    if (tyrec@todo.eq."gen".or. tyrec@todo.eq."TGN")then
        mongen = global_MonGen(tyrec)
        return mongen
    end if
    if(.not. isatt(tyrec,"threshold"))then
        tyrec@threshold = 33.
    else
        print("threshold "+tyrec@threshold)
    end if
    if(.not. isatt(tyrec,"limit"))then
        tyrec@limit = 999.
    else
        print("limit "+tyrec@limit)
    end if
   
    nx = floattoint(360./reso )
    ny = floattoint(180./reso )
    nt = floattoint((ye-yb+1)*12)
    lats = latGlobeFo(ny,"lats","latitude","degrees_north")
    lons = lonGlobeFo(nx,"lons","longitude","degrees_east")
    monrvke = new((/nt,ny,nx/),"float")
    monrvke!0 = "time"
    monrvke!1 = "lat"
    monrvke!2 = "lon"
    monrvke&time = yyyymm
    monrvke&lat = lats
    monrvke&lon = lons

    monrvke = 0.
    oldy = 0.
    oldm = 0.
    oldd = 0.
    oldh = 0.
    oldlat = 0.
    oldlon = 0.
    got = 0
    do i = 0, dims(0)-1
        if(False .and. i.lt.dims(0)-1 .and. .not. any(ismissing(tyrec(i,1:)).or.ismissing(tyrec(i+1,1:))) .and. i.ne.dims(0)-1 .and. (tyrec(i,HH).eq.tyrec(i+1,HH)).and.(tyrec(i,DD).eq.tyrec(i+1,DD)))then
            print("something susp. in i "+i)
            print("i  : "+tyrec(i,0)+" "+tyrec(i,YYYY)+":"+tyrec(i,MM)+":"+tyrec(i,DD)+":"+tyrec(i,HH))
            print("i+1: "+tyrec(i+1,0)+" "+tyrec(i+1,YYYY)+":"+tyrec(i+1,MM)+":"+tyrec(i+1,DD)+":"+tyrec(i+1,HH))
        end if
        if(     .not.  any(ismissing(tyrec(i,1:))) \
           .and.  tyrec(i,VMAX).ge.tyrec@threshold       \  ; some reasons...
           .and.  tyrec(i,VMAX).le.tyrec@limit           \  ; some reasons...
           .and.   (     tyrec(i,YYYY).ne.oldy    \
                    .or. tyrec(i,MM  ).ne.oldm    \
                    .or. tyrec(i,DD  ).ne.oldd    \
                    .or. tyrec(i,HH  ).ne.oldh    \
                    .or. tyrec(i,LAT ).ne.oldlat  \
                    .or. tyrec(i,LON ).ne.oldlon  \
                   ) \
           )then
            ix = global_ix(tyrec(i,LON),reso)
            iy = global_iy(tyrec(i,LAT),reso)
            it = floattoint((tyrec(i,YYYY)-yb)*12  + tyrec(i,MM)) -1
            monrvke(it,iy,ix) = monrvke(it,iy,ix)+rvke(tyrec(i,VMAX),tyrec@todo)
                if (isatt(tyrec,"debug"))then
                    print(tyrec(i,VMAX)+" "+tyrec(i,YYYY)+" "+tyrec(i,TYNUM)+" "+tyrec(i,LAT)+"/"+tyrec(i,LON)+" "+tyrec(i,YYYY)+" "+tyrec(i,MM)+" "+tyrec(i,DD)+" "+tyrec(i,HH))
                end if
                if (     tyrec(i,LON).ge.160 .and. tyrec(i,LON).le.170 \
                    .and.tyrec(i,LAT).ge.20  .and. tyrec(i,LAT).le. 25 )then
                    ;print(tyrec(i,VMAX)+" "+tyrec(i,YYYY)+" "+tyrec(i,TYNUM)+" "+tyrec(i,LAT)+"/"+tyrec(i,LON)+" "+tyrec(i,YYYY)+" "+tyrec(i,MM)+" "+tyrec(i,DD)+" "+tyrec(i,HH)+" rvke:"+rvke(tyrec(i,VMAX),tyrec@todo))
                    ;print("it,iy,ix: "+it+","+iy+","+ix)
                end if
            got = 1
        end if
            ;oldy = tyrec(i,YYYY)
            ;oldm = tyrec(i,MM)
            ;oldd = tyrec(i,DD)
            ;oldh = tyrec(i,HH)
            ;oldlat = tyrec(i,LAT)
            ;oldlon = tyrec(i,LON)
                if (isatt(tyrec,"debug"))then
                    print(tyrec(i,VMAX)+" "+tyrec(i,YYYY)+" "+tyrec(i,TYNUM)+" "+tyrec(i,LAT)+"/"+tyrec(i,LON)+" "+tyrec(i,YYYY)+" "+tyrec(i,MM)+" "+tyrec(i,DD)+" "+tyrec(i,HH)+" got="+got)
                end if
            got = 0
    end do
    if(tyrec@todo.eq."ke" .or. tyrec@todo.eq."rvke")then
        ;monrvke = monrvke/1.
    end if
    if(tyrec@todo.eq."pdi")then
        ;monrvke = monrvke/1.
    end if
    return monrvke
end


undef("timeseries")
function timeseries(ty[*][*][*])
begin
    ;; ty(monthly time,lat,lon) -> ts(yearly time)
    sumty = dim_sum_Wrap(dim_sum_Wrap(ty))
    ravety= runave(sumty,12,0)
    yearlyty = ravety(5::12)
    yearlyty = yearlyty*12
    yearlyty!0 = "year"
    yearlyty&year = ty&time(5::12)
    return yearlyty
end

undef("global_DailyGen")
function global_DailyGen(tyrec[*][*],yyyymmdd[*],limit,threshold)
begin
    load "res_parameter.ncl"
    yb = floattoint(min(tyrec(:,YYYY)))
    ye = floattoint(max(tyrec(:,YYYY)))
    dims = dimsizes(tyrec)
    reso = 5.
    nx = floattoint(360./reso )
    ny = floattoint(180./reso )
    iy = 0
    lats = latGlobeFo(ny,"lats","latitude","degrees_north")
    lons = lonGlobeFo(nx,"lons","longitude","degrees_east")
    dailygen = new((/dimsizes(yyyymmdd),ny,nx/),"float")
    dailygen!0 = "time"
    dailygen!1 = "lat"
    dailygen!2 = "lon"
    dailygen&time = yyyymmdd
    dailygen&lat = lats
    dailygen&lon = lons
    dailygen = 0
    tyend     = 0
    do i = 0, dims(0)-1
        if(tyend .ge. i)then
            continue
        end if
        nrec = min(ind((tyrec(i,YYYY).ne.tyrec(i:,YYYY) .or. tyrec(i,TYNUM).ne.tyrec(i:,TYNUM))))
        nrec = nrec-1
        tyend = i+nrec
        if(ismissing(nrec))then
            tyend = dims(0)-1
        end if
        ;;print(i+":"+tyend+"  "+nrec)
        maxvmax = max(tyrec(i:tyend,VMAX))
        genrec  = min(ind(tyrec(i:tyend,VMAX).ge.34.)) +i
        if(ismissing(genrec))then
            ti = i+1
            ttyend = tyend+1
            ;print("err genrec between "+ti+":"+ttyend)
            continue
        end if
        geny    = floattoint(tyrec(genrec,YYYY))
        genm    = floattoint(tyrec(genrec,MM))
        gend    = floattoint(tyrec(genrec,DD))
        genlat  = tyrec(genrec,LAT)
        genlon  = tyrec(genrec,LON)
        dofy    = day_of_year(geny,genm,gend)
        gendate = geny*10000+genm*100+gend
        if(typeof(iy).ne.typeof(global_iy(genlat,reso)))then
            print(typeof(iy))
            print(typeof(global_iy(genlat,reso)))
            print(genlat+"  "+ti+" "+genrec+" "+ttyend)
            print(min(ind(tyrec(i:tyend,VMAX).ge.34.)))
            exit
        end if
        iy      = global_iy(genlat,reso)
        ix      = global_ix(genlon,reso)
        if(maxvmax.ge.threshold .and. maxvmax.le.limit)then
            dailygen({gendate},iy,ix) = dailygen({gendate},iy,ix)+1
        end if
    end do
    return dailygen
end

undef("global_DailyRvKE")
function global_DailyRvKE(tyrec[*][*])
begin
    load "res_parameter.ncl"
    reso = 5.
    yb = floattoint(min(tyrec(:,YYYY)))
    ye = floattoint(max(tyrec(:,YYYY)))
    yyyymmdd = yyyymmdd_time(yb, ye, "integer")
    dims = dimsizes(tyrec)
    if(.not. isatt(tyrec,"todo"))then
        tyrec@todo = "rvke"
    else
        print("cal "+tyrec@todo)
    end if
    if(.not. isatt(tyrec,"threshold"))then
        tyrec@threshold = 33.
    else
        print("threshold "+tyrec@threshold)
    end if
    if(.not. isatt(tyrec,"limit"))then
        tyrec@limit = 999.
    else
        print("limit "+tyrec@limit)
    end if
    if (tyrec@todo.eq."gen")then
        dailygen = global_DailyGen(tyrec,yyyymmdd,tyrec@limit,tyrec@threshold)
        return dailygen
    end if
   
    nx = floattoint(360./reso )
    ny = floattoint(180./reso )
    nt = dimsizes(yyyymmdd)   ;floattoint((ye-yb+1)*12)
    lats = latGlobeFo(ny,"lats","latitude","degrees_north")
    lons = lonGlobeFo(nx,"lons","longitude","degrees_east")
    dailyrvke = new((/nt,ny,nx/),"float")
    dailyrvke!0 = "time"
    dailyrvke!1 = "lat"
    dailyrvke!2 = "lon"
    dailyrvke&time = yyyymmdd
    dailyrvke&lat = lats
    dailyrvke&lon = lons

    dailyrvke = 0.
    oldy = 0.
    oldm = 0.
    oldd = 0.
    oldh = 0.
    oldlat = 0.
    oldlon = 0.
    got = 0
    do i = 0, dims(0)-1
        if(     .not.  any(ismissing(tyrec(i,1:))) \
           .and.  tyrec(i,VMAX).ge.tyrec@threshold       \  ; some reasons...
           .and.  tyrec(i,VMAX).le.tyrec@limit           \  ; some reasons...
           .and.   (     tyrec(i,YYYY).ne.oldy    \
                    .or. tyrec(i,MM  ).ne.oldm    \
                    .or. tyrec(i,DD  ).ne.oldd    \
                    .or. tyrec(i,HH  ).ne.oldh    \
                    .or. tyrec(i,LAT ).ne.oldlat  \
                    .or. tyrec(i,LON ).ne.oldlon  \
                   ) \
           )then
            ix = global_ix(tyrec(i,LON),reso)
            iy = global_iy(tyrec(i,LAT),reso)
            nowdate = tyrec(i,YYYY)+sprinti("%2.2d",floattoint(tyrec(i,MM)))+sprinti("%2.2d",floattoint(tyrec(i,DD)))
            it = ind(yyyymmdd.eq.stringtoint(nowdate))
            dailyrvke(it,iy,ix) = dailyrvke(it,iy,ix)+rvke(tyrec(i,VMAX),tyrec@todo)
                if (isatt(tyrec,"debug"))then
                    print(tyrec(i,VMAX)+" "+tyrec(i,YYYY)+" "+tyrec(i,TYNUM)+" "+tyrec(i,LAT)+"/"+tyrec(i,LON)+" "+tyrec(i,YYYY)+" "+tyrec(i,MM)+" "+tyrec(i,DD)+" "+tyrec(i,HH))
                end if
            got = 1
        end if
                if (isatt(tyrec,"debug"))then
                    print(tyrec(i,VMAX)+" "+tyrec(i,YYYY)+" "+tyrec(i,TYNUM)+" "+tyrec(i,LAT)+"/"+tyrec(i,LON)+" "+tyrec(i,YYYY)+" "+tyrec(i,MM)+" "+tyrec(i,DD)+" "+tyrec(i,HH)+" got="+got)
                end if
            got = 0
    end do
    if(tyrec@todo.eq."ke" .or. tyrec@todo.eq."rvke")then
        dailyrvke = dailyrvke/10000.
    end if
    return dailyrvke
end

undef("vmax_accu")
function vmax_accu(tyrec[*][*],years,mons)
begin
    load "res_parameter.ncl"  ;  TYNUM, YYYY, MM, DD, HH, LAT, LON, VMAX
    minVmax = 35 ; knots
    maxVmax = floattoint(max(tyrec(:,VMAX)))
    Vmaxs   = ispan(minVmax,maxVmax,5)  ; unit in 5 knots 
    counts  = Vmaxs
    counts  = 0
    counts!0 = "vmax"
    counts&vmax = Vmaxs
    ;; right year -> logical year
    lyear = tyrec(:,YYYY).eq.years(0)
    do i =1, dimsizes(years)-1
        lyear = lyear.or.tyrec(:,YYYY).eq.years(i)
    end do
    ;; right month -> logical mon
    lmon = tyrec(:,MM).eq.mons(0)
    do i =1, dimsizes(mons)-1
        lmon = lmon.or.tyrec(:,MM).eq.mons(i)
    end do
    do i =0, dimsizes(counts)-1
        v  = Vmaxs(i)
        counts({v}) = num((      tyrec(:,VMAX).eq.v)        \
                           .and.(lyear)                     \
                           .and.(lmon)                      )
    end do
    return counts
end

undef("emanuelCorrecting")
function emanuelCorrecting(tyrec[*][*]) ;; correcting method from Emanuel 2005 in Nature
begin
    ;; from Emanuel 2005 in Nature
    ;; Vp = 0.228*V**(1.288)   in m/s unit
    ;; Vpp = w*V + (1-w)*Vp 
    ;;     where w =    0.3        0.8       1.0
    ;;               1959-1966  1967-1972 1973-now

    load "res_parameter.ncl"  ;  TYNUM, YYYY, MM, DD, HH, LAT, LON, VMAX
    ;; knots to m/s
    tyrecm = tyrec
    tyrecm(:,VMAX) = tyrecm(:,VMAX)*0.51

    ;; all Vp
    allV  = tyrecm(:,VMAX)
    allVp = 0.228*(allV^1.288)
    yyyyall = tyrecm(:,YYYY)

    tyrecmvpp = tyrecm
    tyrecmvpp(:,VMAX) = tyrecmvpp@_FillValue
    tyrecmvpp(:,VMAX) = where( yyyyall.ge.1959.and.yyyyall.le.1966, \
                               allV*0.3 + 0.7*allVp, tyrecmvpp(:,VMAX))
    tyrecmvpp(:,VMAX) = where( yyyyall.ge.1967.and.yyyyall.le.1972, \
                               allV*0.8 + 0.2*allVp, tyrecmvpp(:,VMAX))
    tyrecmvpp(:,VMAX) = where( yyyyall.ge.1973, \
                               allV                , tyrecmvpp(:,VMAX))

    ;; back to knots
    tyrecmvpp(:,VMAX) = tyrecmvpp(:,VMAX)/0.51

    return tyrecmvpp
end
