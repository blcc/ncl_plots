external emf "$HOME/blc_plot_function/PI_wrapit/pcmin_2013_wrapit.so"

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
        P = tofloat(Plev*.01)
    else
        P = tofloat(Plev*1.)
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

undef("cal_vpotLL")
function cal_vpotLL(T[*][*][*],R[*][*][*],PSL[*][*],SST[*][*],Plev[*])
begin
    ;; frontend of cal_vpot() for data in (z,y,x)
    dims = dimsizes(T)
    if( any(dims.ne.dimsizes(R)).or. \
        any(dims(1:).ne.dimsizes(PSL)).or.\
        any(dims(1:).ne.dimsizes(SST)).or.\
        (dims(0).ne.dimsizes(Plev)))then
            print("cal_vpotLL(): dimensions not match")
        exit
    end if

    vpot = new(dims(1:),typeof(T))
    capee = new(dims(1:),typeof(T))
    capem = new(dims(1:),typeof(T))
    cape_ast = new(dims(1:),typeof(T))
    do j = 0,dims(1)-1
    do i = 0,dims(2)-1
        vpot1 = cal_vpot(T(:,j,i),R(:,j,i),PSL(j,i),SST(j,i),Plev)
        vpot(j,i) = (/vpot1/)
        capee(j,i) = vpot1@capee
        capem(j,i) = vpot1@capem
        cape_ast(j,i) = vpot1@cape_ast
    end do
    end do
    copy_VarCoords(PSL,vpot)
    vpot@capee = capee
    vpot@capem = capem
    vpot@cape_ast = cape_ast
    return vpot
end
