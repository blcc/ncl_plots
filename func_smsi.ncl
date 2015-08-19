load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
undef("wnponset")
function wnponset(opentaRain[*],pentau[*][*],pentav[*][*])
begin
        ;;1)低層風向出現西南風或南風,高層吹偏東風, and > 5 knots;
        ;;2)降水量高於6mm/day,
        ;;3)兩項條件均持續二個候（pentad）以上,則視為亞洲夏季季風已建立,並以第一候當作季風肇始日期;
        ;;4)反之,若兩個條件其中之一持續消失達6個候以上,則視為季風結束,並以季風建立期間的最後一侯當作季風的結束日期.
    nt = dimsizes(opentaRain)   ; penta
    pentaRain = fftfilter(opentaRain,12)
    lowlev  = 850.
    highlev = 100.
    onoffset = new(2,"integer")
    do i = 1, nt-1-6
        velo = sqrt((pentau(i:i+5,{lowlev})*pentau(i:i+5,{lowlev}))+(pentav(i:i+5,{lowlev})*pentav(i:i+5,{lowlev})))
        if(      pentaRain(i)       .ge.6.      \
            .and.pentaRain(i+1)     .ge.6.      \
            .and.pentau(i,{lowlev}) .ge.0.      \
            .and.pentav(i,{lowlev}) .ge.0.      \
            .and.pentau(i,{highlev}).le.0.      \
            .and.velo(0)            .ge.5.      \
            .and.velo(1)            .ge.5.      \
            .and.(ismissing(onoffset(0)))       \
            )then
            onoffset(0) = i
        end if

        if(    ( (.not.any(pentaRain(i:i+5).ge.6.))            \
            .or.(.not.any(pentau(i:i+5,{lowlev}).ge.0.))       \
            .or.(.not.any(pentav(i:i+5,{lowlev}).ge.0.))       \
            .or.(.not.any(pentau(i:i+5,{highlev}).le.0.))      \
            .or.(.not.any(velo.ge.5.))                    )    \
            .and.ismissing(onoffset(1))                        \
            .and.(.not.ismissing(onoffset(0)))                 \
            )then
            onoffset(1) = i
        end if
    end do
    return onoffset
end

undef("RRi")
function RRi(pentaRain)
begin
    ;; INPUT: pentaRain should be daily rain in 1 year
    ;; OUTPUT: smoothed pentad rain data
    ;; relative pentad mean rainfall rate, Wang and LinHo 2002
    ;; reconstructed pentad series which consists of the long-term mean and the first 12 harmonics
    ;; RRi = Ri - R_JAN

    R_JAN = avg(pentaRain(0:5))  ; JAN avg rain  mm/day

    fftrain = ezfftf(pentaRain)
    fftrain(:,12:) = 0     ; only first 12 harmonics
    Ri = ezfftb(fftrain,fftrain@xbar)

    rRRi = Ri - R_JAN
    return rRRi
end

undef("smsi")
function smsi(pentaRain[*],pentau[*][*],pentav[*][*])
begin
    ;; INPUT: pentaRain should be daily rain in 1 year(73)
    ;;        or monsoon season(not 73)
    ;; OUTPUT: Summer Monsoon Strength Index(SMSI)
    ;; SMSI = sum(pentaRain(onset:offset))
    
    nt = dimsizes(pentaRain)
    if (nt.gt.73)then
        print("pentaRain for smsi should be one year or monsoon season")
        exit
    end if

    if (nt.eq.73)then
        onoffset = wnponset(pentaRain,pentau,pentav)
        onoffp = onoffset +1
        print("onset/offset: "+onoffp(0)+"/"+onoffp(1))
        if(any(ismissing(onoffset)))then
            return 0
        end if
        monsoonRain = pentaRain(onoffset(0):onoffset(1))
    else
        monsoonRain = pentaRain
    end if

    res = sum(monsoonRain)
    return res
end

undef("fftfilter")
function fftfilter(var[*],nfft)
begin
    fftvar = ezfftf(var)
    fftvar(:,nfft:) = 0     ; only first nfft harmonics
    fed = ezfftb(fftvar,fftvar@xbar)
    return fed
end 

