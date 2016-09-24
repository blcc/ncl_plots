load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "func_read_era_all.ncl"
;load "func_read_ersst.ncl"
;load "func_read_cfsr.ncl"
external emf "/data1/home/pgchiu/blc_plot_function/PI_wrapit/pcmin_2013_wrapit.so"

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

undef("cal_vpot")
function cal_vpot(T[*],R[*],PSL[1],SST[1],Plev[*])  ;; temperature(C),mixing ratio(g/kg),Psfc(hPa),SST
begin
    ; SUBROUTINE PCMIN(SST,PSL, P, T, R,NA,N,PMIN,VMAX,IFL)
    ; see readme.txt    r   r  rd rd rd  i i   r    r   i
    ; emf::PCMIN(SST,PSL, P, T, R,NA,N,PMIN,VMAX,IFL)
    
    ;inputs
    if (any(T.ge.100))then
        Tc = T-273.15
        Tc@units = "degree C"
    else
        Tc = T
    end if
    Tk = Tc+273.15
    if (.not.(ismissing(SST).or.ismissing(R(0))))then
        if (any(SST.ge.100))then
            SSTc = SST-273.15
            T@units = "degree C"
        else
            SSTc = SST
        end if
    else 
        resu = new(1,"float")
        resu@capee = resu
        resu@capem = resu
        resu@cape_ast = resu
        return resu
    end if
    SSTk = SSTc+273.15

    if (max(Plev).ge.10000)then ; if Pa then to hPa
        P = Plev*.01
    else
        P = Plev*1.
    end if
    if (PSL.ge.5000)then ; if Pa then to hPa
        PSL = PSL*.01
    end if

    if (isMonotonic(P).eq.1)then ;; from bottom to top
        Tc = Tc(::-1)
        Tk = Tc+273.15
        R  = R(::-1)
        P  = P(::-1)
    end if
        
    NA = dimsizes(P)
    N = NA
    
    ;outputs
    PMIN = new(1,"float")
    VMAX = new(1,"float")
    IFL  = new(1,"integer")
    vpot = new(1,"float")
    P := tofloat(P)

    SIG =  0.0
    CAPEE = new(1,"float")
    CAPEM = new(1,"float")
    CAPES = new(1,"float")
    TOB   = new(1,"float")
    IFLAG = new(1,"integer")
    ;cal this points
    emf::PCMIN(SSTc,PSL, P, Tc, R,NA,N,PMIN,VMAX,IFL)

    ;; for CAPEs
    Rgg = R/1000

    ;; env CAPE
    emf::CAPE(Tk(0),Rgg(0),P(0),Tk,Rgg,P,NA,N,SIG,CAPEE,TOB,IFLAG)

    ;; CAPE at radius of maximum winds
    ppmin = min((/PMIN,1000/))
    RP=0.622*Rgg(0)*PSL/(ppmin*(0.622+Rgg(0))-Rgg(0)*PSL)
    emf::CAPE(Tk(0),RP,ppmin ,Tk,Rgg,P,NA,N,SIG,CAPEM,TOB,IFLAG)

    ;; saturation CAPE at radius of maximum winds
    ppmin = min((/PMIN,1000/))
    ES0 = 6.112*exp(17.67*SSTc/(243.5+SSTc))
    RP=0.622*ES0/(ppmin-ES0)
    emf::CAPE(SSTk,RP ,ppmin,Tk,Rgg,P,NA,N,SIG,CAPES,TOB,IFLAG)


    if(IFL.ne.1 .and. IFL.ne.0)then    ; 1 is ok, 0 is hypercane(not converge)
        ;print("IFL = "+IFL)
    end if
    if (IFL .eq.1 )then  ; 1 is ok, 0 is hypercane(not converge)
        vpot = VMAX
    else
        vpot = vpot@_FillValue
    end if
    
    vpot@long_name = "tropical cyclone intensity potential"
    if (IFL .eq.1 )then  ; 1 is ok, 0 is hypercane(not converge)
        vpot@capee = CAPEE
        vpot@capem = CAPEM
        vpot@cape_ast = CAPES
    else
        vpot@capee = vpot@_FillValue
        vpot@capem = vpot@_FillValue
        vpot@cape_ast = vpot@_FillValue
    end if

    return vpot 
end

function cal_vpot_LL(T[*][*][*],R[*][*][*],PSL[*][*],SST[*][*],Plev[*],Rlist[1]:logical)  ;; temperature(C),mixing ratio(g/kg),Psfc(hPa),SST
begin
    ;; assume T(z,lat,lon) and so on.
    ;; Rlist means return list
    dims = dimsizes(T) ;; assume all data get same dims
    vpot = new(dims((/1,2/)),typeof(T))
    copy_VarCoords(PSL,vpot)
    capee = vpot
    capem = vpot
    cape_ast = vpot

    do j = 0, dims(1)-1
    do i = 0, dims(2)-1
        if(ismissing(PSL(j,i)).or.ismissing(SST(j,i)))then
            continue
        end if
        ibot = ind(.not.ismissing(T(:,j,i)))
        if(any(ismissing(ibot)))then
            continue
        end if
        vpot(j,i) = cal_vpot(T(ibot:,j,i),R(ibot:,j,i),PSL(j,i),SST(j,i),Plev(ibot:))
        capee(j,i) = vpot@capee
        capem(j,i) = vpot@capem
        cape_ast(j,i) = vpot@cape_ast
    end do
    end do

    if(Rlist)then
        return [/vpot,capee,capem,cape_ast/]
    else
        return vpot
    end if
end

function cal_vpot_TLL(T[*][*][*][*],R[*][*][*][*],PSL[*][*][*],SST[*][*][*],Plev[*],Rlist[1]:logical)  ;; temperature(C),mixing ratio(g/kg),Psfc(hPa),SST
begin
    ;; Rlist not done yet
    vpot = SST
    vpot = vpot@_FillValue
    dims = dimsizes(T)
    do t = 0, dims(0)-1
        ;print("cal_vpot_TLL(): t = "+t)
        vpot(t,:,:) = cal_vpot_LL(T(t,:,:,:),R(t,:,:,:),PSL(t,:,:),SST(t,:,:),Plev,False)
    end do
    copy_VarCoords(T(:,0,:,:),vpot)

    return vpot
end

        ;; obsolute ;;undef("read_vpotDaily")
        ;; obsolute ;;function read_vpotDaily(year[1]:integer)
        ;; obsolute ;;begin
        ;; obsolute ;;    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 
        ;; obsolute ;;    cache = "gpi/Vpot_daily_"+year+".nc"
        ;; obsolute ;;    if(isfilepresent(cache))then
        ;; obsolute ;;        df = addfile(cache,"r")
        ;; obsolute ;;        vpot = df->vpot
        ;; obsolute ;;        vpot@capee = df->capee
        ;; obsolute ;;        vpot@capem = df->capem
        ;; obsolute ;;        vpot@cape_ast = df->cape_ast
        ;; obsolute ;;        return vpot
        ;; obsolute ;;    end if
        ;; obsolute ;;    ny = dimsizes(year)
        ;; obsolute ;;    dataset = ""
        ;; obsolute ;;    if(isatt(year,"dataset"))then
        ;; obsolute ;;        dataset = years@dataset
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    vn = "Q"
        ;; obsolute ;;        Q = read_ERAdataDaily3D(year,vn)
        ;; obsolute ;;
        ;; obsolute ;;    R = Q
        ;; obsolute ;;    R = Q/(1-Q)
        ;; obsolute ;;    R = R*1000
        ;; obsolute ;;    R@units = "g kg**-1"
        ;; obsolute ;;    R@long_name = "Mixing ratio"
        ;; obsolute ;;
        ;; obsolute ;;    vn = "T"
        ;; obsolute ;;        T = read_ERAdataDaily3D(year,vn)
        ;; obsolute ;;
        ;; obsolute ;;    vn = "P"
        ;; obsolute ;;        SLP = read_ERAdataDaily(year,vn)
        ;; obsolute ;;
        ;; obsolute ;;    Plev = T&$T!1$
        ;; obsolute ;;
        ;; obsolute ;;    lat = T&$T!2$
        ;; obsolute ;;    lon = T&$T!3$
        ;; obsolute ;;    ny  = dimsizes(lat)
        ;; obsolute ;;    nx  = dimsizes(lon)
        ;; obsolute ;;
        ;; obsolute ;;    oSST = read_oisstdaily(year)
        ;; obsolute ;;    SST = linint2_Wrap(oSST&lon,oSST&lat,oSST(:,0,:,:),True,lon,lat,0)
        ;; obsolute ;;    sstdim = dimsizes(SST)
        ;; obsolute ;;
        ;; obsolute ;;    vpot = SLP
        ;; obsolute ;;    vpot = vpot@_FillValue
        ;; obsolute ;;    capem = vpot
        ;; obsolute ;;    cape_ast = vpot
        ;; obsolute ;;    capee = vpot
        ;; obsolute ;;    vpot@long_name = "Potential intensity of Tropical Cyclone"
        ;; obsolute ;;    capem@long_name = "radius of max wind CAPE in vpot"
        ;; obsolute ;;    cape_ast@long_name = "radius of max wind saturate CAPE in vpot"
        ;; obsolute ;;    capee@long_name = "env. CAPE in vpot"
        ;; obsolute ;;        print("cal vpot") ;; daily cache needed.
        ;; obsolute ;;    do m = 0,sstdim(0)-1
        ;; obsolute ;;        mm = m+1
        ;; obsolute ;;        daycache = "gpi/vpot_"+year+"-"+sprinti("%3.3d",mm)+".nc"
        ;; obsolute ;;        if(isfilepresent(daycache))then
        ;; obsolute ;;            df = addfile(daycache,"r")
        ;; obsolute ;;            vpot(m,:,:)  = df->vpot
        ;; obsolute ;;            capee(m,:,:) = df->capee
        ;; obsolute ;;            capem(m,:,:) = df->capem
        ;; obsolute ;;            cape_ast(m,:,:) = df->cape_ast
        ;; obsolute ;;            continue
        ;; obsolute ;;        end if
        ;; obsolute ;;        print("cal vpot day: "+mm+"/"+sstdim(0))
        ;; obsolute ;;        system("date")
        ;; obsolute ;;        do j = 0,ny-1
        ;; obsolute ;;        do i = 0,nx-1
        ;; obsolute ;;            vpot(m,j,i) = cal_vpot(T(m,:,j,i),R(m,:,j,i),SLP(m,j,i),SST(m,j,i),Plev)
        ;; obsolute ;;            capee(m,j,i) = vpot@capee
        ;; obsolute ;;            capem(m,j,i) = vpot@capem
        ;; obsolute ;;            cape_ast(m,j,i) = vpot@cape_ast
        ;; obsolute ;;
        ;; obsolute ;;            delete(vpot@capem)
        ;; obsolute ;;            delete(vpot@capee)
        ;; obsolute ;;            delete(vpot@cape_ast)
        ;; obsolute ;;        end do
        ;; obsolute ;;        end do
        ;; obsolute ;;        df = addfile(daycache,"c")
        ;; obsolute ;;        df->vpot  = vpot(m,:,:)
        ;; obsolute ;;        df->capee = capee(m,:,:)
        ;; obsolute ;;        df->capem = capem(m,:,:)
        ;; obsolute ;;        df->cape_ast = cape_ast(m,:,:)
        ;; obsolute ;;
        ;; obsolute ;;    end do
        ;; obsolute ;;    
        ;; obsolute ;;    ;vpot!0 = "time"
        ;; obsolute ;;    ;vpot&time = ispan(1,sstdim(0),1)
        ;; obsolute ;;
        ;; obsolute ;;    df = addfile(cache,"c")
        ;; obsolute ;;    df->vpot  = vpot
        ;; obsolute ;;    df->capee = capee
        ;; obsolute ;;    df->capem = capem
        ;; obsolute ;;    df->cape_ast = cape_ast
        ;; obsolute ;;
        ;; obsolute ;;    vpot@capee = capee
        ;; obsolute ;;    vpot@capem = capem
        ;; obsolute ;;    vpot@cape_ast = cape_ast
        ;; obsolute ;;
        ;; obsolute ;;    return vpot
        ;; obsolute ;;end 
        ;; obsolute ;;
        ;; obsolute ;;undef("read_vpotclm12")
        ;; obsolute ;;function read_vpotclm12(years[1])
        ;; obsolute ;;begin
        ;; obsolute ;;    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 
        ;; obsolute ;;    ny = dimsizes(years)
        ;; obsolute ;;    dataset = ""
        ;; obsolute ;;    if(isatt(years,"dataset"))then
        ;; obsolute ;;        dataset = years@dataset
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    vn = "Q"
        ;; obsolute ;;    if(dataset .eq. "cfsr")then
        ;; obsolute ;;        Q = read_cfsr_monly(years,vn)
        ;; obsolute ;;    else
        ;; obsolute ;;        Q = read_ERAdataMon3Dclm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    R = Q
        ;; obsolute ;;    R = Q/(1-Q)
        ;; obsolute ;;    R = R*1000
        ;; obsolute ;;    R@units = "g kg**-1"
        ;; obsolute ;;    R@long_name = "Mixing ratio"
        ;; obsolute ;;
        ;; obsolute ;;    vn = "T"
        ;; obsolute ;;    if(dataset .eq. "cfsr")then
        ;; obsolute ;;        T = read_cfsr_monly(years,vn)
        ;; obsolute ;;    else
        ;; obsolute ;;        T = read_ERAdataMon3Dclm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    vn = "P"
        ;; obsolute ;;    if(dataset .eq. "cfsr")then
        ;; obsolute ;;        SLP = read_cfsr_monly(years,"SLP")
        ;; obsolute ;;    else
        ;; obsolute ;;        SLP = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    Plev = T&$T!1$
        ;; obsolute ;;
        ;; obsolute ;;    lat = T&$T!2$
        ;; obsolute ;;    lon = T&$T!3$
        ;; obsolute ;;    ny  = dimsizes(lat)
        ;; obsolute ;;    nx  = dimsizes(lon)
        ;; obsolute ;;
        ;; obsolute ;;    oSST = read_ersst(years)
        ;; obsolute ;;    ;oSST = read_hadsst2(years)
        ;; obsolute ;;    clmoSST = clmMonTLL(oSST)
        ;; obsolute ;;    clmoSST&month = ispan(1,12,1)
        ;; obsolute ;;    SST = linint2_Wrap(oSST&lon,oSST&lat,clmoSST,True,lon,lat,0)
        ;; obsolute ;;    if(dataset .eq. "eraimsstP2")then
        ;; obsolute ;;        SST = SST+2.
        ;; obsolute ;;    end if
        ;; obsolute ;;    if(dataset .eq. "eraimsstP1")then
        ;; obsolute ;;        SST = SST+1.
        ;; obsolute ;;    end if
        ;; obsolute ;;    if(dataset .eq. "eraimsstM1")then
        ;; obsolute ;;        SST = SST-1.
        ;; obsolute ;;    end if
        ;; obsolute ;;    if(dataset .eq. "eraimsstM2")then
        ;; obsolute ;;        SST = SST-2.
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    vpot = SLP
        ;; obsolute ;;    vpot = vpot@_FillValue
        ;; obsolute ;;    capem = vpot
        ;; obsolute ;;    cape_ast = vpot
        ;; obsolute ;;    capee = vpot
        ;; obsolute ;;    vpot@long_name = "Potential intensity of Tropical Cyclone"
        ;; obsolute ;;    capem@long_name = "radius of max wind CAPE in vpot"
        ;; obsolute ;;    cape_ast@long_name = "radius of max wind saturate CAPE in vpot"
        ;; obsolute ;;    capee@long_name = "env. CAPE in vpot"
        ;; obsolute ;;        print("cal vpot")
        ;; obsolute ;;    do j = 0,ny-1
        ;; obsolute ;;    do i = 0,nx-1
        ;; obsolute ;;    do m = 0,11
        ;; obsolute ;;        vpot(m,j,i) = cal_vpot(T(m,:,j,i),R(m,:,j,i),SLP(m,j,i),SST(m,j,i),Plev)
        ;; obsolute ;;        capee(m,j,i) = vpot@capee
        ;; obsolute ;;        capem(m,j,i) = vpot@capem
        ;; obsolute ;;        cape_ast(m,j,i) = vpot@cape_ast
        ;; obsolute ;;        delete(vpot@capem)
        ;; obsolute ;;        delete(vpot@capee)
        ;; obsolute ;;        delete(vpot@cape_ast)
        ;; obsolute ;;    end do
        ;; obsolute ;;    end do
        ;; obsolute ;;    end do
        ;; obsolute ;;    
        ;; obsolute ;;    vpot@capee = capee
        ;; obsolute ;;    vpot@capem = capem
        ;; obsolute ;;    vpot@cape_ast = cape_ast
        ;; obsolute ;;    vpot!0 = "time"
        ;; obsolute ;;    vpot&time = yyyymm_time(years,years,"integer")
        ;; obsolute ;;
        ;; obsolute ;;    return vpot
        ;; obsolute ;;end 
        ;; obsolute ;;
        ;; obsolute ;;undef("read_vpotclm12_clmsst")
        ;; obsolute ;;function read_vpotclm12_clmsst(years[*])
        ;; obsolute ;;begin
        ;; obsolute ;;    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 
        ;; obsolute ;;
        ;; obsolute ;;    vn = "Q"
        ;; obsolute ;;    Q = read_ERAdataMon3Dclm(years,vn)
        ;; obsolute ;;
        ;; obsolute ;;    R = Q
        ;; obsolute ;;    R = Q/(1-Q)
        ;; obsolute ;;    R = R*1000
        ;; obsolute ;;    R@units = "g kg**-1"
        ;; obsolute ;;    R@long_name = "Mixing ratio"
        ;; obsolute ;;
        ;; obsolute ;;    vn = "T"
        ;; obsolute ;;    T = read_ERAdataMon3Dclm(years,vn)
        ;; obsolute ;;
        ;; obsolute ;;    vn = "P"
        ;; obsolute ;;    SLP = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;
        ;; obsolute ;;    Plev = T&$T!1$
        ;; obsolute ;;
        ;; obsolute ;;    lat = T&$T!2$
        ;; obsolute ;;    lon = T&$T!3$
        ;; obsolute ;;    ny  = dimsizes(lat)
        ;; obsolute ;;    nx  = dimsizes(lon)
        ;; obsolute ;;
        ;; obsolute ;;    oSST = read_ersst(years)
        ;; obsolute ;;    ;oSST = read_hadsst2(ispan(1965,2008,1))
        ;; obsolute ;;    clmoSST = clmMonTLL(oSST)
        ;; obsolute ;;    clmoSST&month = ispan(1,12,1)
        ;; obsolute ;;    SST = linint2_Wrap(oSST&lon,oSST&lat,clmoSST,True,lon,lat,0)
        ;; obsolute ;;
        ;; obsolute ;;    vpot = SLP
        ;; obsolute ;;    vpot = vpot@_FillValue
        ;; obsolute ;;    vpot@long_name = "Potential intensity of Tropical Cyclone"
        ;; obsolute ;;    do j = 0,ny-1
        ;; obsolute ;;    do i = 0,nx-1
        ;; obsolute ;;    do m = 0,11
        ;; obsolute ;;        vpot(m,j,i) = cal_vpot(T(m,:,j,i),R(m,:,j,i),SLP(m,j,i),SST(m,j,i),Plev)
        ;; obsolute ;;    end do
        ;; obsolute ;;    end do
        ;; obsolute ;;    end do
        ;; obsolute ;;
        ;; obsolute ;;    return vpot
        ;; obsolute ;;end 
        ;; obsolute ;;
        ;; obsolute ;;undef("read_vpotclm_obslated")
        ;; obsolute ;;function read_vpotclm_obslated(years[*],omons)
        ;; obsolute ;;begin
        ;; obsolute ;;    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 
        ;; obsolute ;;    if(any(omons.eq.0))then
        ;; obsolute ;;        mons = ispan(1,12,1)
        ;; obsolute ;;    else
        ;; obsolute ;;        mons = omons
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    vn = "Q"
        ;; obsolute ;;    Q = read_ERAdataMon3Dclm(years,vn)
        ;; obsolute ;;    Qclm = avginmons(Q,mons)
        ;; obsolute ;;    R = Qclm
        ;; obsolute ;;    R = Qclm/(1-Qclm)
        ;; obsolute ;;    R = R*1000
        ;; obsolute ;;    R@units = "g kg**-1"
        ;; obsolute ;;    R@long_name = "Mixing ratio"
        ;; obsolute ;;
        ;; obsolute ;;    vn = "T"
        ;; obsolute ;;    T = read_ERAdataMon3Dclm(years,vn)
        ;; obsolute ;;    Tclm = avginmons(T,mons)
        ;; obsolute ;;
        ;; obsolute ;;    vn = "P"
        ;; obsolute ;;    SLP = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;    SLPclm = avginmons(SLP,mons)
        ;; obsolute ;;
        ;; obsolute ;;    Plev = T&$T!1$
        ;; obsolute ;;
        ;; obsolute ;;    lat = T&$T!2$
        ;; obsolute ;;    lon = T&$T!3$
        ;; obsolute ;;    ny  = dimsizes(lat)
        ;; obsolute ;;    nx  = dimsizes(lon)
        ;; obsolute ;;
        ;; obsolute ;;    oSST = read_ersst(years)
        ;; obsolute ;;    ;oSST = read_hadsst2(years)
        ;; obsolute ;;    clmoSST = clmMonTLL(oSST)
        ;; obsolute ;;    clmoSST&month = ispan(1,12,1)
        ;; obsolute ;;    SST = linint2_Wrap(oSST&lon,oSST&lat,clmoSST,True,lon,lat,0)
        ;; obsolute ;;    SSTclm = avginmons(SST,mons)
        ;; obsolute ;;
        ;; obsolute ;;    vpot = SLPclm
        ;; obsolute ;;    vpot = vpot@_FillValue
        ;; obsolute ;;    vpot@long_name = "Potential intensity of Tropical Cyclone"
        ;; obsolute ;;    do j = 0,ny-1
        ;; obsolute ;;    do i = 0,nx-1
        ;; obsolute ;;        vpot(j,i) = cal_vpot(Tclm(:,j,i),R(:,j,i),SLPclm(j,i),SSTclm(j,i),Plev)
        ;; obsolute ;;    end do
        ;; obsolute ;;    end do
        ;; obsolute ;;        ;;    system("rm -f vpot.nc")
        ;; obsolute ;;        ;;    of = addfile("vpot.nc","c")
        ;; obsolute ;;        ;;    of->vpot = vpot
        ;; obsolute ;;        ;;    of->Tclm = Tclm
        ;; obsolute ;;        ;;    of->R = R
        ;; obsolute ;;        ;;    of->SLPclm = SLPclm
        ;; obsolute ;;        ;;    of->SSTclm = SSTclm
        ;; obsolute ;;
        ;; obsolute ;;                ;print(max(vpot({0:40},{110:180})))
        ;; obsolute ;;                ;wks = gsn_open_wks("ps","v")
        ;; obsolute ;;                ;res = True
        ;; obsolute ;;                ;res@mpCenterLonF = 180.
        ;; obsolute ;;                ;p = gsn_csm_contour_map_ce(wks,vpot,res)
        ;; obsolute ;;                ;p = gsn_csm_contour_map_ce(wks,Tclm(0,:,:),res)
        ;; obsolute ;;                ;p = gsn_csm_contour_map_ce(wks,R(0,:,:),res)
        ;; obsolute ;;                ;p = gsn_csm_contour_map_ce(wks,SLPclm,res)
        ;; obsolute ;;                ;p = gsn_csm_contour_map_ce(wks,oSST(7,:,:),res)
        ;; obsolute ;;
        ;; obsolute ;;    return vpot
        ;; obsolute ;;end 
        ;; obsolute ;;
        ;; obsolute ;;undef("cal_gpiclm")
        ;; obsolute ;;function cal_gpiclm(years[*],omons) ;; Genesis potential index, Carmago et al. 2007
        ;; obsolute ;;begin
        ;; obsolute ;;    if(any(omons.eq.0))then
        ;; obsolute ;;        mons = ispan(1,12,1)
        ;; obsolute ;;    else
        ;; obsolute ;;        mons = omons
        ;; obsolute ;;    end if
        ;; obsolute ;;    if(isatt(years,"dataset"))then
        ;; obsolute ;;        dataset = years@dataset
        ;; obsolute ;;    else
        ;; obsolute ;;        dataset = "eraim"
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    ;; 850hPa absolute vorticity
        ;; obsolute ;;    vn = "VOR"
        ;; obsolute ;;    vn@lev = 850 ;; hPa
        ;; obsolute ;;    if(dataset.eq."cfsr")then
        ;; obsolute ;;        tmpall = read_cfsr_monly(years,"AVOR")
        ;; obsolute ;;        ETA    = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;        delete(tmpall)
        ;; obsolute ;;    else
        ;; obsolute ;;        VOR = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;        f = VOR
        ;; obsolute ;;        dims = dimsizes(f)
        ;; obsolute ;;        f = 0
        ;; obsolute ;;        pi = 3.14159
        ;; obsolute ;;        do j = 0, dims(1)-1
        ;; obsolute ;;            ;f(:,j,:) = (2.)*(2*pi/86400.) * sin(f&lat(j)*pi/180)
        ;; obsolute ;;            f(:,j,:) = (2.) *0.000072921   * sin(f&lat(j)*pi/180)
        ;; obsolute ;;                  ;;    2    omega          sin(lat)
        ;; obsolute ;;        end do
        ;; obsolute ;;        ETA = VOR
        ;; obsolute ;;        ETA = VOR + f
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    ;; 700hPa relative humidity
        ;; obsolute ;;    vn = "RH"
        ;; obsolute ;;    vn@lev = 700
        ;; obsolute ;;    if(dataset.eq."cfsr")then
        ;; obsolute ;;        tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;        RH    = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;        delete(tmpall)
        ;; obsolute ;;    else
        ;; obsolute ;;        RH = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    ;; Tropical cyclone potential intensity
        ;; obsolute ;;    Vpot12 = read_vpotclm12(years)
        ;; obsolute ;;    Vpot   = avginmons(Vpot12,mons)
        ;; obsolute ;;    delete(Vpot12)
        ;; obsolute ;;
        ;; obsolute ;;    ;; 200-850 hPa wind shear magnitude
        ;; obsolute ;;    vn = "U"
        ;; obsolute ;;    vn@lev = 200
        ;; obsolute ;;    if(dataset.eq."cfsr")then
        ;; obsolute ;;        tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;        Umon    = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;        delete(tmpall)
        ;; obsolute ;;    else
        ;; obsolute ;;        Umon = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;    U = avginmons(Umon,mons)
        ;; obsolute ;;    vn@lev = 850
        ;; obsolute ;;    if(dataset.eq."cfsr")then
        ;; obsolute ;;        tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;        Umonlo = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;        delete(tmpall)
        ;; obsolute ;;    else
        ;; obsolute ;;        Umonlo = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;    U = U - avginmons(Umonlo,mons)
        ;; obsolute ;;
        ;; obsolute ;;    vn = "V"
        ;; obsolute ;;    vn@lev = 200
        ;; obsolute ;;    if(dataset.eq."cfsr")then
        ;; obsolute ;;        tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;        Vmon   = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;        delete(tmpall)
        ;; obsolute ;;    else
        ;; obsolute ;;        Vmon = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;    V = avginmons(Vmon,mons)
        ;; obsolute ;;    vn@lev = 850
        ;; obsolute ;;    if(dataset.eq."cfsr")then
        ;; obsolute ;;        tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;        Vmonlo = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;        delete(tmpall)
        ;; obsolute ;;    else
        ;; obsolute ;;        Vmonlo = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;    end if
        ;; obsolute ;;    V = V - avginmons(Vmonlo,mons)
        ;; obsolute ;;    Vshear = V
        ;; obsolute ;;    Vshear = sqrt(U*U+V*V)
        ;; obsolute ;;
        ;; obsolute ;;    ;; GPI terms
        ;; obsolute ;;    cETA = avginmons(ETA,mons)
        ;; obsolute ;;    cRH  = avginmons(RH,mons)
        ;; obsolute ;;    cVpot= Vpot
        ;; obsolute ;;    cVshear = Vshear
        ;; obsolute ;;
        ;; obsolute ;;    cETA = abs(100000*cETA)^1.5
        ;; obsolute ;;    cRH  = (cRH/50)^3
        ;; obsolute ;;    cVpot= (cVpot/70)^3
        ;; obsolute ;;    cVshear = (1+0.1*cVshear)^(-2)
        ;; obsolute ;;
        ;; obsolute ;;    ;; GPI
        ;; obsolute ;;    GPI = Vpot  ;; for coordinates
        ;; obsolute ;;    GPI = GPI@_FillValue
        ;; obsolute ;;    GPI = cETA * cRH * cVpot * cVshear
        ;; obsolute ;;
        ;; obsolute ;;    dimsa = dimsizes(GPI)
        ;; obsolute ;;
        ;; obsolute ;;    GPI@long_name = "Genesis potential index"
        ;; obsolute ;;    if(dimsizes(years).eq.1.and.dimsizes(mons).eq.12)then
        ;; obsolute ;;        delete(U)
        ;; obsolute ;;        delete(V)
        ;; obsolute ;;        vn = "U"
        ;; obsolute ;;        vn@lev = 200
        ;; obsolute ;;        if(dataset.eq."cfsr")then
        ;; obsolute ;;            tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;            U      = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;            delete(tmpall)
        ;; obsolute ;;        else
        ;; obsolute ;;            U = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;        end if
        ;; obsolute ;;
        ;; obsolute ;;        vn@lev = 850
        ;; obsolute ;;        if(dataset.eq."cfsr")then
        ;; obsolute ;;            tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;            U      = U - tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;            delete(tmpall)
        ;; obsolute ;;        else
        ;; obsolute ;;            U = U - read_ERAdataMonClm(years,vn)
        ;; obsolute ;;        end if
        ;; obsolute ;;        vn = "V"
        ;; obsolute ;;        vn@lev = 200
        ;; obsolute ;;        if(dataset.eq."cfsr")then
        ;; obsolute ;;            tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;            V      = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;            delete(tmpall)
        ;; obsolute ;;        else
        ;; obsolute ;;            V = read_ERAdataMonClm(years,vn)
        ;; obsolute ;;        end if
        ;; obsolute ;;        vn@lev = 850
        ;; obsolute ;;        if(dataset.eq."cfsr")then
        ;; obsolute ;;            tmpall = read_cfsr_monly(years,vn)
        ;; obsolute ;;            V      = V - tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        ;; obsolute ;;            delete(tmpall)
        ;; obsolute ;;        else
        ;; obsolute ;;            V = V - read_ERAdataMonClm(years,vn)
        ;; obsolute ;;        end if
        ;; obsolute ;;
        ;; obsolute ;;        cVshear12 = V
        ;; obsolute ;;        cVshear12 = sqrt(U*U+V*V)
        ;; obsolute ;;
        ;; obsolute ;;        cETA12 = ETA
        ;; obsolute ;;        cRH12  = RH
        ;; obsolute ;;        cVpot12 = read_vpotclm12(years)
        ;; obsolute ;;
        ;; obsolute ;;        cETA12    = abs(100000*cETA12)^1.5
        ;; obsolute ;;        cRH12     = (cRH12 /50)^3
        ;; obsolute ;;        cVpot12   = (cVpot12/70)^3
        ;; obsolute ;;        cVshear12 = (1+0.1*cVshear12)^(-2)
        ;; obsolute ;;
        ;; obsolute ;;        gpi12     = cVpot12
        ;; obsolute ;;        gpi12     = cETA12 * cRH12 * cVpot12 * cVshear12
        ;; obsolute ;;        gpi12@long_name = "Genesis potential index"
        ;; obsolute ;;
        ;; obsolute ;;        cgpifile = "gpi/gpi"+years+"comp.nc"
        ;; obsolute ;;        if(isatt(years,"dataset"))then
        ;; obsolute ;;            cgpifile = "gpi/gpi"+years@dataset+years+"comp.nc"
        ;; obsolute ;;        end if
        ;; obsolute ;;
        ;; obsolute ;;        system("rm -f "+cgpifile)
        ;; obsolute ;;        print("cgpi: "+cgpifile)
        ;; obsolute ;;        of = addfile(cgpifile,"c")
        ;; obsolute ;;        of->cETA = cETA12
        ;; obsolute ;;        of->cRH  = cRH12
        ;; obsolute ;;        of->cVpot= cVpot12
        ;; obsolute ;;        of->cVshear = cVshear12
        ;; obsolute ;;        of->cgpi = gpi12
        ;; obsolute ;;        return gpi12
        ;; obsolute ;;    end if
        ;; obsolute ;;
        ;; obsolute ;;    return GPI
        ;; obsolute ;;end
        ;; obsolute ;;
        ;; obsolute ;;
