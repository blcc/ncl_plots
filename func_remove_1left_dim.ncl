
function remove_1left_dim(var)
begin ;; turn (1,z,y,x) -> (z,y,x)
    dims = dimsizes(var)
    ndim = dimsizes(dims)
    if(dims(0).ne.1)then
        return var
    end if
    if(ndim.eq.3)then
        return var(0,:,:)
    end if
    if(ndim.eq.4)then
        return var(0,:,:,:)
    end if
    print("remove_1left_dim(): check var dims")
    printVarSummary(var)
    exit
end
