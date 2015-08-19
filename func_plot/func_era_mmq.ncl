;; by pgchiu Jan2010

load "func_read_era_all.ncl"
load "func_dmobs.ncl"

undef("cal_era_mmq_clm")
function cal_era_mmq_clm(years[*]:integer,imons[*]:integer)
begin
    ;; only cal 40S to 40N M and Mq value
    if(any(imons .eq.0))then
        mons = ispan(1,12,1)
    else
        mons = imons
    end if
    nmon = dimsizes(mons)

    vn = "Z"
    z = read_ERAdataMon3Dclm(years,vn)  ;; read as 12month,lev,lat,lon
                                       ;; month should be 1-12, not 0-11

    p = z&$z!1$

    vn = "T"
    t = read_ERAdataMon3Dclm(years,vn)

    vn = "Q"
    q = read_ERAdataMon3Dclm(years,vn)
    if(nmon.eq.1)then
        print("reading and cal m mq of "+years+" "+mons)
        z1 = z({mons},:,:,:)
        t1 = t({mons},:,:,:)
        q1 = q({mons},:,:,:)
    else
        z1 = dim_avg_n_Wrap(z({mons},:,:,:),0)
        t1 = dim_avg_n_Wrap(t({mons},:,:,:),0)
        q1 = dim_avg_n_Wrap(q({mons},:,:,:),0)
    end if

    opt = True
    mmq = cal_m_mq2d(p                 \
                  ,z1(:,{-40:40},{:})  \
                  ,t1(:,{-40:40},{:})  \
                  ,q1(:,{-40:40},{:})  \
                  ,opt)
    pblktop = mmq@pblktop
    copy_VarCoords(mmq,pblktop)

        ;;        system("rm -f mmq.nc")
        ;;        of = addfile("mmq.nc","c")
        ;;        of->m = mmq(0,:,:)
        ;;        of->mq = mmq(1,:,:)
        ;;        of->pbl = pblktop(0,:,:)
        ;;        of->ktop= pblktop(1,:,:)
    return mmq
end
