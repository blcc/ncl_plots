
function ttest_data(alldata[*][*][*],seldata[*][*][*])
begin
    ;; assume data are in (t,lat,lon)
    ;; return in significance, less than 0.1 means good difference.

    avg1 = dim_avg_n(alldata,0)
    avg2 = dim_avg_n(seldata,0)
    var1 = dim_variance_n(alldata,0)
    var2 = dim_variance_n(seldata,0)
    s    = dimsizes(alldata)
    s1   = s(0)
    s    = dimsizes(seldata)
    s2   = s(0)

    iflag = False
    prob = ttest(avg1,var1,s1,avg2,var2,s2,iflag,False)

    return prob  
end
