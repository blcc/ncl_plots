
function ifisatt(var,attname,default)
begin
    if(isatt(var,attname))then
        rr = var@$attname$
    else
        rr = default
    end if
    return rr
end
