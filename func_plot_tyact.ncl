load "func_dots.ncl"
load "func_read_rvke_sst.ncl"
load "func_read_era_all.ncl"

undef("plot_tyact")
function plot_tyact(years[*],season:string,title:string,filename:string)
begin
    ny = dimsizes(years)
    res  = True
        res@title = title
        res@filename = filename
        res@dotmax  = 0.3
    if(isatt(years,"clm"))then
        clm = True
        res@dotmax  = 0.6
    else
        clm = False
    end if

    ;; reading ssta from CPC data
    sstamon = readnino34mon(True)        ; time coord. is YYYYMM in float
    if(season.eq."JJASON")then
        sstasea = month_to_season(sstamon,"JJA")
        sstasea = (sstasea + month_to_season(sstamon,"SON"))/2
    else
        sstasea = month_to_season(sstamon,season)
    end if
    sstasea&time = floattoint(sstasea&time/100)      ; trans to YYYY in int
    ssta    = sstasea({years})             ; data really needed.

    ;; reading 850hPa height for WPSH
    vn      = "Z"
    vn@lev  = 850
    zmon    = read_ERAdataMonly(years,vn)
    zmonclm = clmMonTLL(zmon)
    if(season.eq."JJASON")then
        z       = month_to_season(zmonclm,"JJA")
        z       = (z + month_to_season(zmonclm,"SON"))/2
    else
        z       = month_to_season(zmonclm,season)
    end if
    z = z/9.8  ; geopotional height to geopotional height meter

    ;; reading RACE for dot plot. Should be ano.
    race    = readmonclmRACE(years)
    if(.not.clm)then
        raceclm = readmonclmRACE(ispan(1965,2008,1))
        race    = race - raceclm
    end if
    if(season.eq."JJASON")then
        racesea = month_to_season(race,"JJA")
        racesea = (racesea + month_to_season(race,"SON"))/2
    else    
        racesea = month_to_season(race,season)
    end if
    racesea = racesea/10000.

    ;; reading u,v for wind barb
    vn      = "U"
    vn@lev  = 850
    umon    = read_ERAdataMonly(years,vn)
    umonclm = clmMonTLL(umon)
    if(season.eq."JJASON")then
        u       = month_to_season(umonclm,"JJA")
        u       = (u + month_to_season(umonclm,"SON"))/2
    else
        u       = month_to_season(umonclm,season)
    end if
    vn      = "V"
    vmon    = read_ERAdataMonly(years,vn)
    vmonclm = clmMonTLL(vmon)
    if(season.eq."JJASON")then
        v       = month_to_season(vmonclm,"JJA")
        v       = (v + month_to_season(vmonclm,"SON"))/2
    else
        v       = month_to_season(vmonclm,season)
    end if
    ;; ploting
    a = plot_ssta_wpsh_race_wind(ssta,z(0,:,:),racesea(0,:,:),u(0,:,:),v(0,:,:),res)

    return True
end
