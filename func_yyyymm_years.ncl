
undef("yyyymm_years")
function yyyymm_years(years[*]:integer,months[*]:integer)
begin
    ny = dimsizes(years)
    nm = dimsizes(months)
    y4m2 = new(ny*nm,"integer")
    do y = 0,ny-1
    do m = 0,nm-1
        y4m2(y*nm+m) = years(y)*100+months(m)
    end do
    end do
    return y4m2
end
undef("y4m2toY")
function y4m2toY(y4m2:integer)
begin
    return toint(y4m2/100)
end
undef("y4m2toM")
function y4m2toM(y4m2:integer)
begin
    return toint(mod(y4m2,100))
end
