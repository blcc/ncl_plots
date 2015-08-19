undef("yyyymm_years")
function yyyymm_years(years[*])
begin
    ny = dimsizes(years)
    nt = ny*12
    ym = new(nt,"integer")
    allmons = ispan(1,12,1)
    do y = 0,ny-1
        ym(y*12:y*12+11) = years(y)*100+allmons
    end do
    return ym
end


