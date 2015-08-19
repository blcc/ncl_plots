undef("put_array_into_bigger_one")
function put_array_into_bigger_one(big,ix,n[1]:integer)
begin
    xdims = dimsizes(ix)
    nxdims = dimsizes(xdims)
    bdims = dimsizes(big)
    nbdims = dimsizes(bdims)
    if((nbdims-1).eq.nxdims .and. all(xdims.eq.bdims(1:)))then
        x = ix
    else
        ;; interpolate ix to part of big
        ;; not done yet
        if(nxdims.eq.1)then
            x = linint1_Wrap(ix&$ix!0$,ix,False,big&$big!1$,0)
        end if
        if(nxdims.eq.2)then
            x = linint2_Wrap(ix&$ix!1$,ix&$ix!0$,ix,False,big&$big!2$,big&$big!1$,0)
        end if
        if(.not.isvar("x"))then
            print("put_array_into_bigger_one(): ix ndim is not 1 or 2 and not equal to big one")
            printVarSummary(ix)
            print(kkk)
            exit
        end if
    end if

    ndim = dimsizes(xdims)

    if(ndim.eq.1)then
        big(n,:) = x
    end if
    if(ndim.eq.2)then
        big(n,:,:) = x
    end if
    if(ndim.eq.3)then
        big(n,:,:,:) = x
    end if
    if(ndim.eq.4)then
        big(n,:,:,:,:) = x
    end if
    return big
end
undef("add_a_dim_left")
function add_a_dim_left(x,n)
begin
    dims = dimsizes(x)
    ndim = dimsizes(dims)
    newdims = array_append_record(n,dims,0)
    newx = new(newdims,typeof(x))
    newx = put_array_into_bigger_one(newx,x,0)
    return newx
end
