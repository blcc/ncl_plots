undef("expdist")
function expdist(x[*],ib)
begin
    ;; generate exponential distribution by given x and beta
    ;; beta(b) should be scalar and float/double
    ;; f(x) = exp(-x/b)
    if(typeof(ib).eq."integer")then
        b = int2flt(ib)
    else
        b = ib
    end if

    y = exp(-x/b)

    
    if(isdimnamed(x,0))then
        y!0 = x!0
        y&$y!0$ = x&$x!0$
    end if
    return y
end


undef("fitTyCountsExp")
function fitTyCountsExp(counts[*])
begin
    ;; fit exponential distribution by given x and beta
    ;; counts(0) must .eq. 1
    ;; return best b
    vmaxs = counts&$counts!0$
    x = int2flt((vmaxs-35))/35
    x!0 = "x"
    x&x = x
    olderr = 100.
    oldb   = 0.
    do i = -99, 99
        b = 1. + int2flt(i)/100
        ga= expdist(x,b)
        ga = ga*b
        err = counts-ga
        gaerr = sum(abs(err))/dimsizes(counts)
        if(olderr.gt.gaerr)then
            olderr = gaerr
            oldb   = b
        else
            ;print("rmse/beta: "+olderr+"/"+oldb)
            break
        end if
        ;print("rmse/beta: "+gaerr+"/"+b)
    end do
    return oldb
end
