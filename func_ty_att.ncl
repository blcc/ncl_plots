load "func_isnotmissing.ncl"
undef ("Get_ty_genesis")
function Get_ty_genesis(tydata[*][*][*][8])
;;tydata = new((/maxDataYear,maxTy,maxDataNum,8/),"float")    ;; num,YYYY,MM,DD,HH,lat,lon,Vmax
begin
    dims = dimsizes(tydata)
    maxyear = dims(0)
    maxty   = dims(1)
	maxdata = dims(2)
	tygen   = new((/maxyear,maxty,8/),"float")
	yeardim = tydata!0
	tygen!0 = yeardim
	tygen&$yeardim$ = tydata&$yeardim$
	do y = 0 , maxyear-1
    	do ty = 0, maxty-1
			ypt = tydata(y,ty,:,5)
			xpt = tydata(y,ty,:,6)
			vmax= tydata(y,ty,:,7)
			do i = 0, maxdata-1
				if((.not.ismissing(vmax(i))).and.(vmax(i).ge.34.))then
					tygen(y,ty,:) = tydata(y,ty,i,:)
					break
				end if
			end do
    	end do
    end do
	
	return tygen
end

undef("aveGenLocat")
function aveGenLocat(tygen[*][*][8])
begin
	dsize = dimsizes(tygen)
	maxyear = dsize(0)
	maxty   = dsize(1)
	aveloc = new((/5,3/),"float") ; (all,PhS,SCS,other,missing),(lat,lon,count)
	aveloc = 0.
	aveloc(0,0) = avg(tygen(:,:,5))
	aveloc(0,1) = avg(tygen(:,:,6))
	aveloc(0,2) = maxyear+maxty
	do y = 0, maxyear-1
	do ty = 0, maxty-1
		lon = tygen(y,ty,6)
		lat = tygen(y,ty,5)
		basin = GenBasin(lat,lon)
		if (basin .eq. "PhS")then
			aveloc(1,0) = aveloc(1,0) + lat
			aveloc(1,1) = aveloc(1,1) + lon
			aveloc(1,2) = aveloc(1,2) + 1
		end if
		if (basin .eq. "SCS")then
			aveloc(2,0) = aveloc(2,0) + lat
			aveloc(2,1) = aveloc(2,1) + lon
			aveloc(2,2) = aveloc(2,2) + 1
		end if
		if (basin .eq. "OutBound")then
			aveloc(3,0) = aveloc(3,0) + lat
			aveloc(3,1) = aveloc(3,1) + lon
			aveloc(3,2) = aveloc(3,2) + 1
		end if
		if (basin .eq. "Missing")then
			aveloc(4,2) = aveloc(4,2) + 1
		end if
	end do
	end do
	do i = 1,2
	do j = 0,1
		if (.not.(aveloc(i,2).eq.0))then
			aveloc(i,j) = aveloc(i,j)/aveloc(i,2)
		end if
	end do
	end do
	return aveloc
end

undef("GenBasin")
function GenBasin(lat,lon)
begin
	basin = "OutBound"
	if (ismissing(lat) .or. ismissing(lon))then
		basin = "Missing"
	else
		;;if ((lon.gt.125.).and.(lon.le.150.))then
		if ((lon.gt.125.).and.(lon.le.210.))then
			basin = "PhS"
		end if
		if ((lon.le.120.) .and. (lon.ge.110.))then
			basin = "SCS"
		end if
	end if
	return basin
end

undef("countGBasin")
function countGBasin(lats,lons)
begin
	maxty    = dimsizes(lats)
	maxlons  = dimsizes(lons)
	if(maxty .ne. maxlons)then
        ;print("function countGBasin(lats,lons), dimsizes mismatch")
	end if

	Basins = new(5,"integer") ;; all, PhS, SCS, Other, missing
	Basins = 0
	Basins!0 = "Counts"
	BNames = (/"All","PhS","SCS","OutBound","Missing"/)
	Basins&Counts = BNames

	do i = 0, maxty-1
		Basin = GenBasin(lats(i),lons(i))
		if(Basin.ne."Missing")then
			Basins(0) = Basins(0) +1
		end if
		if (Basin.eq."PhS")then
			Basins(1) = Basins(1) +1
		end if
		if (Basin.eq."SCS")then
			Basins(2) = Basins(2) +1
		end if
		if (Basin.eq."OutBound")then
			Basins(3) = Basins(3) +1
		end if
		if (Basin.eq."Missing")then
			Basins(4) = Basins(4) +1
		end if
	end do
	return Basins
end

undef("filterMonths")
function filterMonths(tydata[*][*][*][8],mm[*])
begin
    dims = dimsizes(tydata)
    maxyear = dims(0)
    maxty   = dims(1)
	maxdata = dims(2)
	do y = 0, maxyear-1
	do ty = 0, maxty-1
		if((.not.ismissing(tydata(y,ty,0,0)) .and. (ismissing(ind(mm.eq.tydata(y,ty,0,2) )))))then
			tydata(y,ty,:,:) = tydata@_FillValue
		end if
	end do
	end do
	return tydata
end

undef("filterBasins")
function filterBasins(tydata[*][*][*][8],basin)
begin
    dims = dimsizes(tydata)
    maxyear = dims(0)
    maxty   = dims(1)
	maxdata = dims(2)
	do y = 0, maxyear-1
	do ty = 0, maxty-1
    do i = 0, maxdata-1
		if(.not.ismissing(tydata(y,ty,i,6)))then
            lat = tydata(y,ty,i,5)
            lon = tydata(y,ty,i,6)
            ;;if(isgenbasin(lat,lon,basin))then
            ;;print(basin+" needed "+GenBasin(lat,lon)+" got")
            ;;print("lat/lon:"+lat+"/"+lon)
            if(GenBasin(lat,lon).eq.basin)then
                ;;print("Right")
                break
            end if
			tydata(y,ty,:,:) = tydata@_FillValue
            break
		end if
	end do
	end do
	end do
	return tydata
end

undef("filterYears")
function filterYears(tydata[*][*][*][8],yy[*])
begin
    dims = dimsizes(tydata)
	tyy  = tydata
    maxyear = dims(0)
    maxty   = dims(1)
	maxdata = dims(2)
	ystart  = tydata(0,0,0,1)
	do y = 0, maxyear-1
	do ty = 0, maxty-1
		if((.not.ismissing(tyy(y,ty,0,0)) .and. (ismissing(ind(yy.eq.tyy(y,ty,0,1) )))))then
			tyy(y,ty,:,:) = tyy@_FillValue
		end if
	end do
	end do
	return tyy
end

undef("filterNeartw")
function filterNeartw(tyrec[*][*][*][8])
begin
    tydata = tyrec
    dims = dimsizes(tydata)
	maxyear = dims(0)
	maxty   = dims(1)
	maxdata = dims(2)
	do y = 0, maxyear-1
	do ty =0, maxty-1
		ypt = tydata(y,ty,:,5)
		xpt = tydata(y,ty,:,6)
		if(.not.isneartaiwan(ypt,xpt))then
			tydata(y,ty,:,:) = tydata@_FillValue
		end if
    end do
    end do
 return tydata
end

undef("filterNearsouth")

function filterNearsouth(tyrec[*][*][*][8])
begin
    tydata = tyrec
    dims = dimsizes(tydata)
	maxyear = dims(0)
	maxty   = dims(1)
	maxdata = dims(2)
	do y = 0, maxyear-1
	do ty =0, maxty-1
		ypt = tydata(y,ty,:,5)
		xpt = tydata(y,ty,:,6)
		if(.not.isnearsouth(ypt,xpt))then
			tydata(y,ty,:,:) = tydata@_FillValue
		end if
    end do
    end do
 return tydata
end

undef("alignTY")
function alignTY(tcrec[*][*][*][8])
begin
	tydata = tcrec
    tydata@_FillValue = -999.
	dims = dimsizes(tydata)
	maxyear	= dims(0)
	maxty	= dims(1)
	maxdata	= dims(2)
	do y = 0, maxyear-1
	do ty= 0, maxty-1
		tymaxI= maxind(tydata(y,ty,:,7))       ; max wind
		if(ismissing(tymaxI))then
				continue
		end if
		do i = tymaxI, maxdata-1
			if((ismissing(tydata(y,ty,i,7))).or.(tydata(y,ty,i,7).lt.35.).or.tydata(y,ty,i,7).ge.200.)then
				tydata(y,ty,i:,:) = tydata@_FillValue
                break
            end if
		end do 
		do j = 0, tymaxI
			if(ismissing(tydata(y,ty,j,7)))then
				tydata(y,ty,:,:) = tydata@_FillValue
				break
			end if
			if (tydata(y,ty,j,7).ge.35.)then
                ptmp = tydata(y,ty,j:i-1,:)
				tydata(y,ty,0:i-j-1,:) = (/ptmp/)
				tydata(y,ty,i-j:,:)    = tydata@_FillValue

				delete(ptmp)
				break
			end if
		end do
	end do
	end do
	return tydata
end

undef("yearly_gen_counts")
function yearly_gen_counts(tygen[*][*][8])
begin
	dims    = dimsizes(tygen)
	maxyear = dims(0)
	maxty   = dims(1)
	counts  = new((/maxyear,5/),"float"); (all,PhS,SCS,other,missing),counts
	yeardim = tygen!0
	counts!0 = yeardim
	counts&$yeardim$ = tygen&$yeardim$
	do y = 0, maxyear-1
		lats   = tygen(y,:,5)
		lons   = tygen(y,:,6)
		yearcounts = countGBasin(lats,lons)
		counts(y,:) = yearcounts
	end do
	return counts
end

undef("tyturns")
function tyturns(tydata[*][*][*][8])
begin
    regN    = 4
    dims    = dimsizes(tydata)
    tyturn  = new(dims(0:2),"float")
    do y = 0, dims(0)-1
    do ty= 0, dims(1)-1
        do i = regN, dims(2)-1-regN
            if (any(ismissing(tydata(y,ty,i-regN:i+regN-1,5))))then
            else
                tyturn(y,ty,i) = traj_angle(tydata(y,ty,i-regN:i+regN-1,5),tydata(y,ty,i-regN:i+regN-1,6))
            end if
        end do
    end do
    end do
    return tyturn
end

undef("slop2angle")
function slop2angle(m1,m2,xdir1,xdir2)
begin
    a1 = atan(m1)
    a2 = atan(m2)
    if(xdir1 .lt.0)then
        a1 = a1+3.14159
    end if
    if(xdir2 .lt.0)then
        a2 = a2+3.14159
    end if
    ang = a2-a1
    ang = ang*180./3.14159
    return ang
end

undef("traj_angle")
function traj_angle(lats[*],lons[*])
begin
    dims = dimsizes(lats)
    if(dims.ne.dimsizes(lons)) then
            exit
    end if
    per  = dims/2
    if(any(ismissing(lats)))then
        angle = new(1,"float")
    else
        m1 = regline(lons(0:per-1),lats(0:per-1))
        xdir1 = lons(per-1) - lons(0)
        m2 = regline(lons(per:),lats(per:))
        xdir2 = lons(dims-1)  - lons(per)
        angle = slop2angle(m1,m2,xdir1,xdir2)
    end if
    return angle
end

undef("Tyrec2TyHr")
function Tyrec2TyHr(tcrec[*][*][*][8])
begin
    reso = 5.
    load "res_map.ncl"
    latB = res@mpMinLatF
    latE = res@mpMaxLatF
    lonB = res@mpMinLonF 
    lonE = res@mpMaxLonF
    nx   = floattoint((lonE-lonB)/reso)+2
    lons = fspan(lonB-(reso/2.),lonE+(reso/2.),nx)
    lons@units = "degree_east"
    ny   = floattoint((latE-latB)/reso)+2
    lats = fspan(latB-(reso/2.),latE+(reso/2.),ny)
    lats@units = "degrees_north"
    tyhr = new((/ny,nx/),"float")
    tyhr = 0.
    tyhr!0 = "lat"
    tyhr&lat = lats
    tyhr!1 = "lon"
    tyhr&lon = lons
    dims = dimsizes(tcrec)
    do y  = 0,dims(0) -1
    do ty = 0,dims(1) -1
    do d  = 0,dims(2) -1
        if(.not.ismissing(tcrec(y,ty,d,5)))then
            ypt = tcrec(y,ty,d,5)
            xpt = tcrec(y,ty,d,6)
            vmax= tcrec(y,ty,d,7)
            j   = floattoint((ypt-latB)/reso) +1
            i   = floattoint((xpt-lonB)/reso) +1
            if(j.lt.ny.and.j.ge.0.and.i.lt.nx.and.i.ge.0.and.vmax.ge.34)then
                tyhr(j,i) = tyhr(j,i) + 6.
            end if
            if(j.eq.4.and.i.eq.2)then
            end if
        end if
    end do
    end do
    end do
    tyhr = tyhr/(2007-1965+1)
    return tyhr
end

undef("Tyrec2TyHrmon")
function Tyrec2TyHrmon(tcrec[*][*][*][8])
begin
    reso = 5.
    load "res_map.ncl"
    latB = res@mpMinLatF
    latE = res@mpMaxLatF
    lonB = res@mpMinLonF 
    lonE = res@mpMaxLonF
    nx   = floattoint((lonE-lonB)/reso)+2
    lons = fspan(lonB-(reso/2.),lonE+(reso/2.),nx)
    lons@units = "degree_east"
    ny   = floattoint((latE-latB)/reso)+2
    lats = fspan(latB-(reso/2.),latE+(reso/2.),ny)
    lats@units = "degrees_north"
    dims = dimsizes(tcrec)
    nt   = dims(0)*12
    tyhr = new((/nt,ny,nx/),"float")
    tyhr = 0.
    tyhr!0 = "time"
    tyhr!1 = "lat"
    tyhr&lat = lats
    tyhr!2 = "lon"
    tyhr&lon = lons
    dims = dimsizes(tcrec)
    do y  = 0,dims(0) -1
    do ty = 0,dims(1) -1
    do d  = 0,dims(2) -1
        if(.not.ismissing(tcrec(y,ty,d,5)))then
            ypt = tcrec(y,ty,d,5)
            xpt = tcrec(y,ty,d,6)
            vmax= tcrec(y,ty,d,7)
            j   = floattoint((ypt-latB)/reso) +1
            i   = floattoint((xpt-lonB)/reso) +1
            if(j.lt.ny.and.j.ge.0.and.i.lt.nx.and.i.ge.0.and.vmax.ge.34)then
                t   = floattoint((y*12) + tcrec(y,ty,d,2))
                tyhr(t,j,i) = tyhr(t,j,i) + 6.
            end if
        end if
    end do
    end do
    end do
    
    dfile = "TyHrmon.nc"
    if (isfilepresent(dfile))then
        system("rm "+dfile)
    end if
    ff = addfile(dfile,"c")
    ff->tyhr = tyhr
    return tyhr
end

undef("Tyrec2TyOC")
function Tyrec2TyOC(tcrec[*][*][*][8])
begin
    dfile = "TyOC.nc"
    reso = 5.
    load "res_map.ncl"
    latB = res@mpMinLatF
    latE = res@mpMaxLatF
    lonB = res@mpMinLonF 
    lonE = res@mpMaxLonF
    nx   = floattoint((lonE-lonB)/reso)+2
    lons = fspan(lonB-(reso/2.),lonE+(reso/2.),nx)
    lons@units = "degree_east"
    ny   = floattoint((latE-latB)/reso)+2
    lats = fspan(latB-(reso/2.),latE+(reso/2.),ny)
    lats@units = "degrees_north"
    tyocp = new((/ny,nx/),"float")
    tyocp!0 = "lat"
    tyocp&lat = lats
    tyocp!1 = "lon"
    tyocp&lon = lons
    tyocp = 0.
    tyoc  = tyocp
    dims = dimsizes(tcrec)
    do y  = 0,dims(0) -1
    do ty = 0,dims(1) -1
    do d  = 0,dims(2) -1
        if(.not.ismissing(tcrec(y,ty,d,5)))then
            ypt = tcrec(y,ty,d,5)
            xpt = tcrec(y,ty,d,6)
            vmax= tcrec(y,ty,d,7)
            j   = floattoint((ypt-latB)/reso) +1
            i   = floattoint((xpt-lonB)/reso) +1
            if(j.lt.ny.and.j.ge.0.and.i.lt.nx.and.i.ge.0.and.vmax.ge.34)then
                tyocp(j,i) = 1.
            end if
        end if
    end do
    tyoc = tyoc + tyocp
    tyocp = 0.
    end do
    end do
    if (isfilepresent(dfile))then
        system("rm "+dfile)
    end if
    ff = addfile(dfile,"c")
    ff->tyoc = tyoc
    return tyoc
end

undef("Tyrec2TyOCmon")
function Tyrec2TyOCmon(tcrec[*][*][*][8])
begin
    dfile = "TyOCmon.nc"
    reso = 5.
    load "res_map.ncl"
    latB = res@mpMinLatF
    latE = res@mpMaxLatF
    lonB = res@mpMinLonF 
    lonE = res@mpMaxLonF
    nx   = floattoint((lonE-lonB)/reso)+2
    lons = fspan(lonB-(reso/2.),lonE+(reso/2.),nx)
    lons@units = "degree_east"
    ny   = floattoint((latE-latB)/reso)+2
    lats = fspan(latB-(reso/2.),latE+(reso/2.),ny)
    lats@units = "degrees_north"
    dims = dimsizes(tcrec)
    nt   = dims(0)*12
    tyocp = new((/nt,ny,nx/),"float")
    tyocp!0 = "time"
    tyocp!1 = "lat"
    tyocp&lat = lats
    tyocp!2 = "lon"
    tyocp&lon = lons
    tyocp = 0.
    tyoc  = tyocp
    dims = dimsizes(tcrec)
    do y  = 0,dims(0) -1
    do ty = 0,dims(1) -1
    do d  = 0,dims(2) -1
        if(.not.ismissing(tcrec(y,ty,d,5)))then
            ypt = tcrec(y,ty,d,5)
            xpt = tcrec(y,ty,d,6)
            vmax= tcrec(y,ty,d,7)
            mm  = tcrec(y,ty,d,2)
            t   = floattointeger(y*12+mm-1)
            j   = floattoint((ypt-latB)/reso) +1
            i   = floattoint((xpt-lonB)/reso) +1
            if(j.lt.ny.and.j.ge.0.and.i.lt.nx.and.i.ge.0.and.vmax.ge.34)then
                 ;print("t= "+t)
                tyocp(t,j,i) = 1.
            end if
        end if
    end do
    tyoc = tyoc + tyocp
    tyocp = 0.
    end do
    end do
    ;tyoc = tyoc/(2007-1965+1)
    if (isfilepresent(dfile))then
        system("rm "+dfile)
    end if
    ff = addfile(dfile,"c")
    ff->tyoc = tyoc
    return tyoc
end

undef("tyRecCont")
function tyRecCont(tcrec[*][*][*][8])
begin
    tycont = 0
    tycont@TyRecuCount = 0
    tcdims = dimsizes(tcrec)
    return tycont
end
undef("isgenbasin")
function isgenbasin (lat[*]:float, lon[*]:float, basin)
local i,j,k,genbasin
begin
    rightbasin = False
    do i = 0, dimsizes(lat)-1
        if (.not.ismissing(lat(i))  \
        )then
            if(GenBasin(lat(i),lon(i).eq.basin))then
                rightbasin = True
            end if
            break
        end if
    end do
    return(rightbasin)
end

undef("rvke")
function rvke(vknot)  ; from Yu et al. 2009, revised ACE index
begin
    v   = vknot
    vc  = 35.
    v   = where(v.lt.vc,v@_FillValue,v)
    a   = 0.51
    rc  = (v/vc)^(1/a)
    ke  = (v*v/(rc*rc)) * (0.5 + ((rc^(2-(2*a))-1)/(1-a)))
    return ke
end

undef("Tyrec2TyKEmon")
function Tyrec2TyKEmon(tcrec[*][*][*][8])
begin
    ;print("Tyrec2TyKEmon")
    dfile = "TyKEmon.nc"
    reso = 5.
    load "res_map.ncl"
    latB = res@mpMinLatF
    latE = res@mpMaxLatF
    lonB = res@mpMinLonF 
    lonE = res@mpMaxLonF
    nx   = floattoint((lonE-lonB)/reso)+2
    lons = fspan(lonB-(reso/2.),lonE+(reso/2.),nx)
    lons@units = "degree_east"
    ny   = floattoint((latE-latB)/reso)+2
    lats = fspan(latB-(reso/2.),latE+(reso/2.),ny)
    lats@units = "degrees_north"
    dims = dimsizes(tcrec)
    nt   = dims(0)*12
    tykep = new((/nt,ny,nx/),"float")
    tykep!0 = "time"
    tykep!1 = "lat"
    tykep&lat = lats
    tykep!2 = "lon"
    tykep&lon = lons
    tykep = 0.
    tyke  = tykep
    dims = dimsizes(tcrec)
    do y  = 0,dims(0) -1
    do ty = 0,dims(1) -1
    do d  = 0,dims(2) -1
        if(.not.ismissing(tcrec(y,ty,d,5)))then
            ypt = tcrec(y,ty,d,5)
            xpt = tcrec(y,ty,d,6)
            vmax= tcrec(y,ty,d,7)
            mm  = tcrec(y,ty,d,2)
            t   = floattointeger(y*12+mm-1)
            j   = floattoint((ypt-latB)/reso) +1
            i   = floattoint((xpt-lonB)/reso) +1
            if(j.lt.ny.and.j.ge.0.and.i.lt.nx.and.i.ge.0.and.vmax.ge.34)then
                ; print("t= "+t)
                tykep(t,j,i) = tykep(t,j,i) + (vmax*vmax)
            end if
        end if
    end do
    tyke = tyke + tykep
    tykep = 0.
    end do
    end do
    ;tyke = tyke/(2007-1965+1)
    if (isfilepresent(dfile))then
        system("rm "+dfile)
    end if
    ff = addfile(dfile,"c")
    ff->tyke = tyke
    return tyke
end
undef("Tyrec2TyGenmon")
function Tyrec2TyGenmon(tcrec[*][*][*][8])
begin
    ;print("Tyrec2TyGenmon")
    dfile = "TyGenmon.nc"
    reso = 5.
    yb   = floattoint(min(tcrec(:,:,:,1)))
    ye   = floattoint(max(tcrec(:,:,:,1)))
    load "res_map.ncl"
    latB = res@mpMinLatF
    latE = res@mpMaxLatF
    lonB = res@mpMinLonF 
    lonE = res@mpMaxLonF
    nx   = floattoint((lonE-lonB)/reso)+2
    lons = fspan(lonB-(reso/2.),lonE+(reso/2.),nx)
    lons@units = "degree_east"
    ny   = floattoint((latE-latB)/reso)+2
    lats = fspan(latB-(reso/2.),latE+(reso/2.),ny)
    lats@units = "degrees_north"
    dims = dimsizes(tcrec)
    times = yyyymm_time(yb,ye,"integer")
    nt   = dims(0)*12
    tygep = new((/nt,ny,nx/),"float")
    tygep!0 = "time"
    tygep&time = times
    tygep!1 = "lat"
    tygep&lat = lats
    tygep!2 = "lon"
    tygep&lon = lons
    tygep = 0.
    tyge  = tygep
    dims = dimsizes(tcrec)
    do y  = 0,dims(0) -1
    do ty = 0,dims(1) -1
    do d  = 0,dims(2) -1
        if(.not.ismissing(tcrec(y,ty,d,5)))then
            ypt = tcrec(y,ty,d,5)
            xpt = tcrec(y,ty,d,6)
            vmax= tcrec(y,ty,d,7)
            mm  = tcrec(y,ty,d,2)
            t   = floattointeger(y*12+mm-1)
            j   = floattoint((ypt-latB)/reso) +1
            i   = floattoint((xpt-lonB)/reso) +1
            if(j.lt.ny.and.j.ge.0.and.i.lt.nx.and.i.ge.0.and.vmax.ge.34)then
                ; print("t= "+t)
                tygep(t,j,i) = tygep(t,j,i) + 1
                break
            end if
        end if
    end do
    tyge = tyge + tygep
    tygep = 0.
    end do
    end do
    if (isfilepresent(dfile))then
        system("rm "+dfile)
    end if
    ff = addfile(dfile,"c")
    ff->tyge = tyge
    return tyge
end

undef("Tyrec2TyRvKEmon")
function Tyrec2TyRvKEmon(tcrec[*][*][*][8])
begin
    ;print("Tyrec2TyRvKEmon")
    dfile = "TyRvKEmon.nc"
    reso = 5.
    load "res_map.ncl"
    latB = res@mpMinLatF
    latE = res@mpMaxLatF
    lonB = res@mpMinLonF 
    lonE = res@mpMaxLonF
    nx   = floattoint((lonE-lonB)/reso)+2
    lons = fspan(lonB-(reso/2.),lonE+(reso/2.),nx)
    lons@units = "degree_east"
    ny   = floattoint((latE-latB)/reso)+2
    lats = fspan(latB-(reso/2.),latE+(reso/2.),ny)
    lats@units = "degrees_north"
    dims = dimsizes(tcrec)
    nt   = dims(0)*12
    tykep = new((/nt,ny,nx/),"float")
    tykep!0 = "time"
    tykep&time = yyyymm_time(1965, 2008, "integer")
    tykep!1 = "lat"
    tykep&lat = lats
    tykep!2 = "lon"
    tykep&lon = lons
    tykep = 0.
    tyke  = tykep
    do y  = 0,dims(0) -1
    do ty = 0,dims(1) -1
    do d  = 0,dims(2) -1
        if(.not.ismissing(tcrec(y,ty,d,5)))then
            ypt = tcrec(y,ty,d,5)
            xpt = tcrec(y,ty,d,6)
            vmax= tcrec(y,ty,d,7)
            mm  = tcrec(y,ty,d,2)
            t   = floattointeger(y*12+mm-1)
            j   = floattoint((ypt-latB)/reso) +1
            i   = floattoint((xpt-lonB)/reso) +1
            if(j.lt.ny.and.j.ge.0.and.i.lt.nx.and.i.ge.0.and.vmax.ge.34)then
                tykep(t,j,i) = tykep(t,j,i) + rvke(vmax)
            end if
        end if
    end do
    tyke = tyke + tykep
    tykep = 0.
    end do
    end do
    ;tyke = tyke/(2007-1965+1)
    if (isfilepresent(dfile))then
        system("rm "+dfile)
    end if
    ff = addfile(dfile,"c")
    ff->tyrvke = tyke
    return tyke
end
