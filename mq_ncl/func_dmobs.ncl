;; by pgchiu Jan2010
;; mod for nuwa Sep2013

;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("myint2p")
function myint2p(ip,ix,op,opt)
begin
    if(all(ismissing(ix)))then
        ox = ix
        ox = ox@_FillValue
        return ox
    end if
    if(isMonotonic(ip).eq.1)then
        iip = ip(::-1)
        iix = ix(::-1)
    else
        iip = ip
        iix = ix
    end if
    ox = int2p_Wrap(iip,iix,op,opt)
    idims = dimsizes(ip)
    odims = dimsizes(op)
    return ox
end

undef("outputmmq")
function outputmmq(mmq[*][*][*],fn[1]:string)
begin
    pblktop = mmq@pblktop
    copy_VarCoords(mmq,pblktop)
    dm  = mmq(0,:,:)
    mq  = mmq(1,:,:)
    pbl = pblktop(0,:,:)
    ktop= pblktop(1,:,:)
    system ("rm -f "+fn)
    df = addfile(fn,"c")
    df->m = dm
    df->mq= mq
    df->pbl = pbl
    df->ktop= ktop
    return fn
end


undef("plot_p_h")
function plot_p_h(ip,ih,iq,tk,ireflev,iktoplev,title,filename)
begin
    h = ih
    h = ih/10000  ;; as a scale
    if(ireflev .gt. 10000)then
        reflev = ireflev/100
        ktoplev = iktoplev/100
    else
        reflev = ireflev
        ktoplev = iktoplev
    end if
    refh = h({reflev})
    p    = h&$h!0$
    mixr = iq
    mixr = iq/(1-iq)

    p@_FillValue = tk@_FillValue
    rh = relhum(tk,mixr,p*100) ;; units in  K, kg/kg, Pa
    rh!0 = "plev"
    rh&plev = ih&$ih!0$
    rlev = ip(10:)
    rlev = rlev(::-1)
    ;print(rh({rlev}))
    rrh = myint2p(rlev,rh({rlev}),p,1)   ;; re-interpolate rh avoid inconsistence
    ;print("rh,tk,mixr,p = "+rh+" , "+tk+" , "+mixr+" , "+p)

    res                   = True
    res@tiMainString      = title
    res@trYReverse        = True
    res@gsnFrame          = False     ; don't advance frame yet
    
    res2 = res
    res2@xyLineColors  = "blue"
    wks = gsn_open_wks("ps",filename)
    ;plot  = gsn_csm_xy (wks,h,p,res)
    ;plot  = gsn_csm_x2y (wks,h,rh,p,res,res2)
    plot  = gsn_csm_x2y (wks,h,rrh,p,res,res2)
    gsn_polyline(wks,plot,(/refh,refh/),(/max(p),min(p)/),False)
    gsn_polyline(wks,plot,(/0,refh/),(/reflev,reflev/),False)
    ;gsn_polyline(wks,plot,(/0,refh/),(/ktoplev,ktoplev/),False)
    frame(wks)
    print("plot p-h "+filename)
    return True
end

undef("find_cb")
function find_cb(tk[*],q[*],lev[1]:string)
begin
    ;; lev = "CCL" or "LCL"
    ;; but all CCL now
    i0 = min(ind(.not.ismissing(tk))) ;; value base
    it = max(ind(.not.ismissing(tk))) ;; value top
    minpbl = 1000.  ;; 900.   min PBL lev in hPa

    mixr = q
    mixr = q/(1-q)

    p = tk&$tk!0$
    p@_FillValue = tk@_FillValue
    rh0 = relhum(tk(i0),mixr(i0),p(i0)*100) ;; units in  K, kg/kg, Pa

    sfcdp = dewtemp_trh(tk(i0),rh0)
    do i = i0,it
        ;print("sfcdp, tk @ : "+sfcdp+" "+tk(i)+" "+p(i))
        if(sfcdp.ge.tk(i))then
            ipbl = i
            break
        end if
    end do
    if(isvar("ipbl"))then
        pbl = p(ipbl)
        if(pbl.ge.minpbl)
            pbl = minpbl
        end if
        ;print("pbl: "+pbl)
    else
        ;print("no pbl")
        pbl = minpbl
        pbl@nopbl = True
    end if
        ;print("pbl: "+pbl)
    return pbl
end

undef("find_ct")
function find_ct(h[*],pbl[1]) ;; cloud top
begin
    ;; h(cloud top) .eq. h(cloud base)
    maxktop = 100. ; 100 hPa means no high limit
    minktop = 400. ; lowest ktop 
    if(pbl .le. maxktop)then
        print("no cloud")
        minktop@noktop = True
        return minktop
    end if
    p= h&$h!0$
    maxiktop = ind(p.eq.maxktop)
    miniktop = ind(p.eq.minktop)
    ipbl = ind(p.eq.pbl)
    hpbl = h(ipbl)

    bottotop = False    ;; find same h from top or from bottom
    iktop = maxiktop
    if(bottotop)then
        do i = ipbl+1,maxiktop
            if(hpbl.le.h(i))then
                iktop = i
                break
            end if
        end do
    else
        iktop = 0
        nz = maxiktop - miniktop +1
        do i = 0,nz-1
            if(hpbl.ge.h(maxiktop-i))then
                iktop = maxiktop - i
                break
            end if
        end do
    end if
    ktop = p(iktop)
    ;print("ktop: "+ktop)
    return ktop
end

undef("cal_M_Mq")
function cal_M_Mq(ip[*],iz[*],it[*],iq[*],opt)
begin                ;;  Estimating the Gross Moist Stability of the Tropical Atomsphere, Jia-Yuh Yu and Chia Chou and J.David Neelin, J. of the Atmospheric Science, 1998, 55, 1354--1372
;..   p      : pressure [hPa] -> must be
;..   z      : geopotental hight [m]
;..   t      : temperature [K]
;..   q      : specific humidity [kg/kg]
    ;; interpolate
    linopt = -2    ;; linopt.eq.1 means linear, .ne.1 is log interpolate, .lt.0 means extrapolation
    P = fspan(1000,100,901) ; from sfc to top of atmosphere
    Z = myint2p(ip,iz,P,linopt)
    T = myint2p(ip,it,P,linopt)
    Q = myint2p(ip,iq,P,linopt)
    if (any(ismissing((/Z,T,Q/))))then
        ;print(T)
        ;print(ip+" "+it)
        ;;print("missing values after interpolate")
        resu = new(2,"float")
        resu = resu@_FillValue
        resu@ktop = resu@_FillValue
        resu@pbl  = resu@_FillValue
        return resu
    end if
    ;; units
    if (max(P) .lt. 5000)then
        P = P*100               ; hPa to Pa
    end if
    if (max(T).lt. 100)then
        T = T+273.15            ; degree C to K
    end if
    if (max(T).gt.400)then     ;;  well...
        ;;print(ip+" "+it)
        resu = new(2,typeof(T))
        resu = resu@_FillValue
        resu@ktop = resu@_FillValue
        resu@pbl  = resu@_FillValue
        return resu
    end if
    if (max(Z&$Z!0$).lt.5000.and.Z({500}).ge.10000)then  ; geopotential height @ hPa coord
        Z = Z/9.81
    end if
    if (max(Z&$Z!0$).gt.5000.and.Z({50000}).ge.10000)then  ; geopotential height @ Pa coord
        Z = Z/9.81
    end if
    ;; constants
    l   = 2.44*10^6     ; J Kg-1 K-1
    cp  = 1005          ; J Kg-1 K-1
    rv  = 461.51        ; J Kg-1 K-1
    g   = 9.81          ; m s-2
    kappa = 0.286       ; R/Cp

    ;; cal ga (as gamma profile)  \gamma = (\frac{dq_{sat}}{dT}_{\overline{T}})
    TCsfc = T(0)-273.15  ; T(0) is sfc temperature
        ;;e0 = 100 * 6.112 * exp((17.67*TCsfc)/(TCsfc+243.5))
    e0 = 611.2 * exp((17.67*TCsfc)/(TCsfc+243.5))
    ga = 0.622*l*l*e0/(rv*cp*T*T*P)*exp(l*((1/T(0))-(1/T))/rv)

    ;; cal moist-static energy profile (h(p))
    h = T
    h@long_name = "moist static energy"
    h = cp*T + l*Q + g*Z
    lq = l*Q
    cpt= cp*T
    gz = g*Z
    ;print("h = cp*T + l*Q + g*Z")
    ;print("h = "+cpt+" + "+lq+" + "+gz+" = "+h+" @P="+P/100+"hPa")

    ;; def cloud top(ktop) & base(pbl)
    pbl  = find_cb(T,Q,"CCL")
    ktop = find_ct(h,pbl)
    pbl  = pbl  * 100            ;hPa to Pa
    ktop = ktop * 100
            ;;if(isatt(pbl,"nopbl"))then    ;; ktop less than 300
            ;;if((pbl-ktop).le. 20)then
                ;print("h @ lev: "+h+" "+p+" hPa")
                ;print("hpbl: "+hpbl)
                ;print("ktop/pbl = "+ktop+" / "+pbl)
                ;plot = plot_p_h(ip,h,Q,T,pbl,ktop,"ktop too high","ph")
                ;sleep(2)
                ;printVarSummary(iz)
                ;;print("ip iz it iq")
                ;;print(ip+" "+iz+" "+it+" "+iq)
                ;exit
            ;;end if
    ipbl = ind(P.eq.pbl)
    iktop= ind(P.eq.ktop)
    if((iktop-ipbl).lt.10.or.isatt(pbl,"nopbl").or.isatt(ktop,"noktop"))then
        resu = new(2,"float")
        resu = 0.
        resu@ktop = 0.
        resu@pbl  = 0.
        return resu
    end if

    ;; vertical average of h between pbl and ktop -> hhat
    hhat = avg(h(ipbl:iktop))
    qhat = avg(Q(ipbl:iktop))*l   ;; as moist energy

    ;; cal A profile
    gasum = sum(1/(1+ga(ipbl:iktop))/P(ipbl:iktop))*100  ; dp = 100 Pa
        ;; \int_p^{p_b}{(1+\gamma)^{-1}d lnp}  interpolate resolution = 100Pa
    a = h   ;; quik create array
    a = a@_FillValue
    ;; above pbl
    a(ipbl:iktop) = (1/(1+ga(ipbl:iktop)))*exp(-1*kappa*gasum)
    ;; below pbl
    a(0:ipbl)     = (1/(1+ga(0:ipbl)))*(P(0:ipbl)/P(ipbl))*kappa
  
    ;; ahat
    ahat =  avg(a(0:iktop))

    ;; aplus = int_p^{p_0}{A(p`)}dlnp`
    aplus = a
    aplus = aplus@_FillValue
    do k = ipbl, iktop
        aplus(k) = sum(a(0:k)/P(0:k)*100)
    end do

    ;; aplushat : vertical average of aplus
    aplushat = avg(aplus(ipbl:iktop))
    
    ;; hahat and qahat: avg(aplus*h) & avg(aplus*q)
    hahat = avg(h(ipbl:iktop)*aplus(ipbl:iktop))
    qahat = avg(Q(ipbl:iktop)*aplus(ipbl:iktop))*l  ;; as moist energy

    ;print("h*aplus: "+h+" * "+aplus+" = "+h*aplus)
    ;print("h = "+cpt+" + "+lq+" + "+gz+" = "+h+" @P="+P/100+"hPa")
    ;print("hahat = "+hahat)
    ;print("aplushat*hhat = "+aplushat*hhat)

    deltam  = hahat - aplushat*hhat
    deltamq = -qahat + (aplushat*qhat)
    if (deltamq .le.0)then
        deltamq = 0
    end if
    if (deltam .le.0)then
        deltam = 0
    end if
    if(deltam .ge. 1000 .or. deltam.lt.0)then ;; should between 0 and 700
        print("m = hahat-aplushat*hhat")
        print("m = "+hahat+" - "+aplushat+" * "+hhat+" = "+deltam)
        ;print("aplus = "+aplus(ipbl:iktop)+" @P="+P(ipbl:iktop))
        print("m out of range")
        deltam = 0
        ;exit
    end if
    if(deltamq .ge. 2800 .or. deltamq.lt.0)then ;; should between 0 and 2400
        print("mq = -qahat + aplushat * qhat")
        print("mq = -"+qahat+" + "+aplushat+" * "+qhat+" = "+deltamq)
        print("mq out of range")
        deltamq = 0
        ;exit
    end if
    resu = new(2,"float")
    resu(0) = deltam
    resu(1) = deltamq
    resu@ktop = ktop/100 ;P(iktop)  ;; hPa
    resu@pbl  = pbl/100  ;P(ipbl)   ;; hPa
    return resu
end

undef("cal_m_mq2d")
function cal_m_mq2d(p[*],z[*][*][*],t[*][*][*],q[*][*][*],opt)
begin
    ;; assume lev*lat*lon
    dims = dimsizes(z)
    ny = dims(1)
    nx = dims(2)
    mmq = new((/4,ny,nx/),"float")
    mmq!1 = "lat"
    mmq!2 = "lon"
    mmq&lat = z&$z!1$
    mmq&lon = z&$z!2$
    iy = ny-1
    do j = 0, ny-1
    cputb = get_cpu_time()
    ;;system("date")
    ;print("j = "+j+" of "+iy)
    do i = 0, nx-1
        ;print("j,i = "+j+","+i+" of "+iy)
        dmmq = cal_M_Mq(p,z(:,j,i),t(:,j,i),q(:,j,i),opt)
        mmq(0,j,i) = dmmq(0)
        mmq(1,j,i) = dmmq(1)
        mmq(2,j,i) = dmmq@pbl
        mmq(3,j,i) = dmmq@ktop
        delete(dmmq)
    end do
    ;print("              in "+(get_cpu_time()-cputb)+" secs")
    end do
    ;system("date")
    ;print("Cal M,Mq done.")

    mmq!0 = "mmq"
    mmq&mmq = (/"M","Mq","pbl","ktop"/)

    return mmq
end

undef("cal_m_mqTLL")
function cal_m_mqTLL(p[*],z[*][*][*][*],t[*][*][*][*],q[*][*][*][*],opt)
begin
    ;; assume z,t,q(time,lev,lat,lon) 
    if(isatt(opt,"monprefix"))
        monprefix = opt@monprefix
    end if
    dims = dimsizes(z)
    nt   = dims(0)
    maxnt = 12 ;; cal will slow down in ncl, so restart needed
    maxnt = min((/maxnt,(96*72*12)/(dims(2)*dims(3))/))  ;; cal will slow down in ncl, so restart needed
    ;maxnt = nt
    datatime = z&$z!0$
    outdims=dims
    outdims(1) = 2
    mmq     = new(outdims,typeof(z))
    mmq!0   = "time"
    mmq&time= datatime
    pblktop = mmq    
    mmq!1 = "mmq"
    mmq&mmq = (/"M","Mq"/)

    calt  = 0
    do time = 0,nt-1
        if(isvar("monprefix"))then
            tmpfn = "ncfiles_mon/"+monprefix+"_"+datatime(time)+".nc"
            if(isfilepresent(tmpfn))then
                df = addfile(tmpfn,"r")
                mmq(time,0,:,:)     = df->m
                mmq(time,1,:,:)     = df->mq
                pblktop(time,0,:,:) = df->pbl
                pblktop(time,1,:,:) = df->ktop
                continue
            end if
        end if
        print("cal "+datatime(time))
        monmmq = cal_m_mq2d(p,z(time,:,:,:),t(time,:,:,:),q(time,:,:,:),opt)
        calt = calt+1
        mmq(time,:,:,:)     = monmmq
        pblktop(time,:,:,:) = monmmq@pblktop
        tmpfn = outputmmq(monmmq,tmpfn)
        delete(monmmq)
        if(calt.ge.maxnt)then
            exit
        end if
    end do
    delete(mmq@pblktop)
    mmq@pblktop = pblktop

    return mmq
end
