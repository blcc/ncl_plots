undef("isrecurve")
function isrecurve(lat[*]:float,lon[*]:float)
begin
    recu = False
    wlon = minind(lon)
    if(ismissing(wlon))then
        return recu
    end if
    if(.not. ismissing(lon(wlon+3)))then
        recu = True
    end if

    return recu
end
