load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "func_read_era_all.ncl"
load "func_read_ersst.ncl"
load "func_read_cfsr.ncl"
external emf "./pcmin_revised_wrapit.so"

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
    if (.not.ismissing(SST))then
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
    RP=0.622*Rgg(0)*PSL/(PSL*(0.622+Rgg(0))-Rgg(0)*PSL)
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

undef("read_vpotDaily")
function read_vpotDaily(year[1]:integer)
begin
    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 
    cache = "gpi/Vpot_daily_"+year+".nc"
    if(isfilepresent(cache))then
        df = addfile(cache,"r")
        vpot = df->vpot
        vpot@capee = df->capee
        vpot@capem = df->capem
        vpot@cape_ast = df->cape_ast
        return vpot
    end if
    ny = dimsizes(year)
    dataset = ""
    if(isatt(year,"dataset"))then
        dataset = years@dataset
    end if

    vn = "Q"
        Q = read_ERAdataDaily3D(year,vn)

    R = Q
    R = Q/(1-Q)
    R = R*1000
    R@units = "g kg**-1"
    R@long_name = "Mixing ratio"

    vn = "T"
        T = read_ERAdataDaily3D(year,vn)

    vn = "P"
        SLP = read_ERAdataDaily(year,vn)

    Plev = T&$T!1$

    lat = T&$T!2$
    lon = T&$T!3$
    ny  = dimsizes(lat)
    nx  = dimsizes(lon)

    oSST = read_oisstdaily(year)
    SST = linint2_Wrap(oSST&lon,oSST&lat,oSST(:,0,:,:),True,lon,lat,0)
    sstdim = dimsizes(SST)

    vpot = SLP
    vpot = vpot@_FillValue
    capem = vpot
    cape_ast = vpot
    capee = vpot
    vpot@long_name = "Potential intensity of Tropical Cyclone"
    capem@long_name = "radius of max wind CAPE in vpot"
    cape_ast@long_name = "radius of max wind saturate CAPE in vpot"
    capee@long_name = "env. CAPE in vpot"
        print("cal vpot") ;; daily cache needed.
    do m = 0,sstdim(0)-1
        mm = m+1
        daycache = "gpi/vpot_"+year+"-"+sprinti("%3.3d",mm)+".nc"
        if(isfilepresent(daycache))then
            df = addfile(daycache,"r")
            vpot(m,:,:)  = df->vpot
            capee(m,:,:) = df->capee
            capem(m,:,:) = df->capem
            cape_ast(m,:,:) = df->cape_ast
            continue
        end if
        print("cal vpot day: "+mm+"/"+sstdim(0))
        system("date")
        do j = 0,ny-1
        do i = 0,nx-1
            vpot(m,j,i) = cal_vpot(T(m,:,j,i),R(m,:,j,i),SLP(m,j,i),SST(m,j,i),Plev)
            capee(m,j,i) = vpot@capee
            capem(m,j,i) = vpot@capem
            cape_ast(m,j,i) = vpot@cape_ast

            delete(vpot@capem)
            delete(vpot@capee)
            delete(vpot@cape_ast)
        end do
        end do
        df = addfile(daycache,"c")
        df->vpot  = vpot(m,:,:)
        df->capee = capee(m,:,:)
        df->capem = capem(m,:,:)
        df->cape_ast = cape_ast(m,:,:)

    end do
    
    ;vpot!0 = "time"
    ;vpot&time = ispan(1,sstdim(0),1)

    df = addfile(cache,"c")
    df->vpot  = vpot
    df->capee = capee
    df->capem = capem
    df->cape_ast = cape_ast

    vpot@capee = capee
    vpot@capem = capem
    vpot@cape_ast = cape_ast

    return vpot
end 

undef("read_vpotclm12")
function read_vpotclm12(years[1])
begin
    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 
    ny = dimsizes(years)
    dataset = ""
    if(isatt(years,"dataset"))then
        dataset = years@dataset
    end if

    vn = "Q"
    if(dataset .eq. "cfsr")then
        Q = read_cfsr_monly(years,vn)
    else
        Q = read_ERAdataMon3Dclm(years,vn)
    end if

    R = Q
    R = Q/(1-Q)
    R = R*1000
    R@units = "g kg**-1"
    R@long_name = "Mixing ratio"

    vn = "T"
    if(dataset .eq. "cfsr")then
        T = read_cfsr_monly(years,vn)
    else
        T = read_ERAdataMon3Dclm(years,vn)
    end if

    vn = "P"
    if(dataset .eq. "cfsr")then
        SLP = read_cfsr_monly(years,"SLP")
    else
        SLP = read_ERAdataMonClm(years,vn)
    end if

    Plev = T&$T!1$

    lat = T&$T!2$
    lon = T&$T!3$
    ny  = dimsizes(lat)
    nx  = dimsizes(lon)

    oSST = read_ersst(years)
    ;oSST = read_hadsst2(years)
    clmoSST = clmMonTLL(oSST)
    clmoSST&month = ispan(1,12,1)
    SST = linint2_Wrap(oSST&lon,oSST&lat,clmoSST,True,lon,lat,0)
    if(dataset .eq. "eraimsstP2")then
        SST = SST+2.
    end if
    if(dataset .eq. "eraimsstP1")then
        SST = SST+1.
    end if
    if(dataset .eq. "eraimsstM1")then
        SST = SST-1.
    end if
    if(dataset .eq. "eraimsstM2")then
        SST = SST-2.
    end if

    vpot = SLP
    vpot = vpot@_FillValue
    capem = vpot
    cape_ast = vpot
    capee = vpot
    vpot@long_name = "Potential intensity of Tropical Cyclone"
    capem@long_name = "radius of max wind CAPE in vpot"
    cape_ast@long_name = "radius of max wind saturate CAPE in vpot"
    capee@long_name = "env. CAPE in vpot"
        print("cal vpot")
    do j = 0,ny-1
    do i = 0,nx-1
    do m = 0,11
        vpot(m,j,i) = cal_vpot(T(m,:,j,i),R(m,:,j,i),SLP(m,j,i),SST(m,j,i),Plev)
        capee(m,j,i) = vpot@capee
        capem(m,j,i) = vpot@capem
        cape_ast(m,j,i) = vpot@cape_ast
        delete(vpot@capem)
        delete(vpot@capee)
        delete(vpot@cape_ast)
    end do
    end do
    end do
    
    vpot@capee = capee
    vpot@capem = capem
    vpot@cape_ast = cape_ast
    vpot!0 = "time"
    vpot&time = yyyymm_time(years,years,"integer")

    return vpot
end 

undef("read_vpotclm12_clmsst")
function read_vpotclm12_clmsst(years[*])
begin
    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 

    vn = "Q"
    Q = read_ERAdataMon3Dclm(years,vn)

    R = Q
    R = Q/(1-Q)
    R = R*1000
    R@units = "g kg**-1"
    R@long_name = "Mixing ratio"

    vn = "T"
    T = read_ERAdataMon3Dclm(years,vn)

    vn = "P"
    SLP = read_ERAdataMonClm(years,vn)

    Plev = T&$T!1$

    lat = T&$T!2$
    lon = T&$T!3$
    ny  = dimsizes(lat)
    nx  = dimsizes(lon)

    oSST = read_ersst(years)
    ;oSST = read_hadsst2(ispan(1965,2008,1))
    clmoSST = clmMonTLL(oSST)
    clmoSST&month = ispan(1,12,1)
    SST = linint2_Wrap(oSST&lon,oSST&lat,clmoSST,True,lon,lat,0)

    vpot = SLP
    vpot = vpot@_FillValue
    vpot@long_name = "Potential intensity of Tropical Cyclone"
    do j = 0,ny-1
    do i = 0,nx-1
    do m = 0,11
        vpot(m,j,i) = cal_vpot(T(m,:,j,i),R(m,:,j,i),SLP(m,j,i),SST(m,j,i),Plev)
    end do
    end do
    end do

    return vpot
end 

undef("read_vpotclm_obslated")
function read_vpotclm_obslated(years[*],omons)
begin
    ;; function cal_vpot(T[*],R[*],PSL,SST,P) 
    if(any(omons.eq.0))then
        mons = ispan(1,12,1)
    else
        mons = omons
    end if

    vn = "Q"
    Q = read_ERAdataMon3Dclm(years,vn)
    Qclm = avginmons(Q,mons)
    R = Qclm
    R = Qclm/(1-Qclm)
    R = R*1000
    R@units = "g kg**-1"
    R@long_name = "Mixing ratio"

    vn = "T"
    T = read_ERAdataMon3Dclm(years,vn)
    Tclm = avginmons(T,mons)

    vn = "P"
    SLP = read_ERAdataMonClm(years,vn)
    SLPclm = avginmons(SLP,mons)

    Plev = T&$T!1$

    lat = T&$T!2$
    lon = T&$T!3$
    ny  = dimsizes(lat)
    nx  = dimsizes(lon)

    oSST = read_ersst(years)
    ;oSST = read_hadsst2(years)
    clmoSST = clmMonTLL(oSST)
    clmoSST&month = ispan(1,12,1)
    SST = linint2_Wrap(oSST&lon,oSST&lat,clmoSST,True,lon,lat,0)
    SSTclm = avginmons(SST,mons)

    vpot = SLPclm
    vpot = vpot@_FillValue
    vpot@long_name = "Potential intensity of Tropical Cyclone"
    do j = 0,ny-1
    do i = 0,nx-1
        vpot(j,i) = cal_vpot(Tclm(:,j,i),R(:,j,i),SLPclm(j,i),SSTclm(j,i),Plev)
    end do
    end do
        ;;    system("rm -f vpot.nc")
        ;;    of = addfile("vpot.nc","c")
        ;;    of->vpot = vpot
        ;;    of->Tclm = Tclm
        ;;    of->R = R
        ;;    of->SLPclm = SLPclm
        ;;    of->SSTclm = SSTclm

                ;print(max(vpot({0:40},{110:180})))
                ;wks = gsn_open_wks("ps","v")
                ;res = True
                ;res@mpCenterLonF = 180.
                ;p = gsn_csm_contour_map_ce(wks,vpot,res)
                ;p = gsn_csm_contour_map_ce(wks,Tclm(0,:,:),res)
                ;p = gsn_csm_contour_map_ce(wks,R(0,:,:),res)
                ;p = gsn_csm_contour_map_ce(wks,SLPclm,res)
                ;p = gsn_csm_contour_map_ce(wks,oSST(7,:,:),res)

    return vpot
end 

undef("cal_gpiclm")
function cal_gpiclm(years[*],omons) ;; Genesis potential index, Carmago et al. 2007
begin
    if(any(omons.eq.0))then
        mons = ispan(1,12,1)
    else
        mons = omons
    end if
    if(isatt(years,"dataset"))then
        dataset = years@dataset
    else
        dataset = "eraim"
    end if

    ;; 850hPa absolute vorticity
    vn = "VOR"
    vn@lev = 850 ;; hPa
    if(dataset.eq."cfsr")then
        tmpall = read_cfsr_monly(years,"AVOR")
        ETA    = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        delete(tmpall)
    else
        VOR = read_ERAdataMonClm(years,vn)
        f = VOR
        dims = dimsizes(f)
        f = 0
        pi = 3.14159
        do j = 0, dims(1)-1
            ;f(:,j,:) = (2.)*(2*pi/86400.) * sin(f&lat(j)*pi/180)
            f(:,j,:) = (2.) *0.000072921   * sin(f&lat(j)*pi/180)
                  ;;    2    omega          sin(lat)
        end do
        ETA = VOR
        ETA = VOR + f
    end if

    ;; 700hPa relative humidity
    vn = "RH"
    vn@lev = 700
    if(dataset.eq."cfsr")then
        tmpall = read_cfsr_monly(years,vn)
        RH    = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        delete(tmpall)
    else
        RH = read_ERAdataMonClm(years,vn)
    end if

    ;; Tropical cyclone potential intensity
    Vpot12 = read_vpotclm12(years)
    Vpot   = avginmons(Vpot12,mons)
    delete(Vpot12)

    ;; 200-850 hPa wind shear magnitude
    vn = "U"
    vn@lev = 200
    if(dataset.eq."cfsr")then
        tmpall = read_cfsr_monly(years,vn)
        Umon    = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        delete(tmpall)
    else
        Umon = read_ERAdataMonClm(years,vn)
    end if
    U = avginmons(Umon,mons)
    vn@lev = 850
    if(dataset.eq."cfsr")then
        tmpall = read_cfsr_monly(years,vn)
        Umonlo = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        delete(tmpall)
    else
        Umonlo = read_ERAdataMonClm(years,vn)
    end if
    U = U - avginmons(Umonlo,mons)

    vn = "V"
    vn@lev = 200
    if(dataset.eq."cfsr")then
        tmpall = read_cfsr_monly(years,vn)
        Vmon   = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        delete(tmpall)
    else
        Vmon = read_ERAdataMonClm(years,vn)
    end if
    V = avginmons(Vmon,mons)
    vn@lev = 850
    if(dataset.eq."cfsr")then
        tmpall = read_cfsr_monly(years,vn)
        Vmonlo = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
        delete(tmpall)
    else
        Vmonlo = read_ERAdataMonClm(years,vn)
    end if
    V = V - avginmons(Vmonlo,mons)
    Vshear = V
    Vshear = sqrt(U*U+V*V)

    ;; GPI terms
    cETA = avginmons(ETA,mons)
    cRH  = avginmons(RH,mons)
    cVpot= Vpot
    cVshear = Vshear

    cETA = abs(100000*cETA)^1.5
    cRH  = (cRH/50)^3
    cVpot= (cVpot/70)^3
    cVshear = (1+0.1*cVshear)^(-2)

    ;; GPI
    GPI = Vpot  ;; for coordinates
    GPI = GPI@_FillValue
    GPI = cETA * cRH * cVpot * cVshear

    dimsa = dimsizes(GPI)

    GPI@long_name = "Genesis potential index"
    if(dimsizes(years).eq.1.and.dimsizes(mons).eq.12)then
        delete(U)
        delete(V)
        vn = "U"
        vn@lev = 200
        if(dataset.eq."cfsr")then
            tmpall = read_cfsr_monly(years,vn)
            U      = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
            delete(tmpall)
        else
            U = read_ERAdataMonClm(years,vn)
        end if

        vn@lev = 850
        if(dataset.eq."cfsr")then
            tmpall = read_cfsr_monly(years,vn)
            U      = U - tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
            delete(tmpall)
        else
            U = U - read_ERAdataMonClm(years,vn)
        end if
        vn = "V"
        vn@lev = 200
        if(dataset.eq."cfsr")then
            tmpall = read_cfsr_monly(years,vn)
            V      = tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
            delete(tmpall)
        else
            V = read_ERAdataMonClm(years,vn)
        end if
        vn@lev = 850
        if(dataset.eq."cfsr")then
            tmpall = read_cfsr_monly(years,vn)
            V      = V - tmpall(:,{vn@lev*100},:,:) ;; vertical unit is Pa 
            delete(tmpall)
        else
            V = V - read_ERAdataMonClm(years,vn)
        end if

        cVshear12 = V
        cVshear12 = sqrt(U*U+V*V)

        cETA12 = ETA
        cRH12  = RH
        cVpot12 = read_vpotclm12(years)

        cETA12    = abs(100000*cETA12)^1.5
        cRH12     = (cRH12 /50)^3
        cVpot12   = (cVpot12/70)^3
        cVshear12 = (1+0.1*cVshear12)^(-2)

        gpi12     = cVpot12
        gpi12     = cETA12 * cRH12 * cVpot12 * cVshear12
        gpi12@long_name = "Genesis potential index"

        cgpifile = "gpi/gpi"+years+"comp.nc"
        if(isatt(years,"dataset"))then
            cgpifile = "gpi/gpi"+years@dataset+years+"comp.nc"
        end if

        system("rm -f "+cgpifile)
        print("cgpi: "+cgpifile)
        of = addfile(cgpifile,"c")
        of->cETA = cETA12
        of->cRH  = cRH12
        of->cVpot= cVpot12
        of->cVshear = cVshear12
        of->cgpi = gpi12
        return gpi12
    end if

    return GPI
end


