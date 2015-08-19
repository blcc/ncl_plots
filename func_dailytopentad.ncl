undef("dailytopentad")
function dailytopentad(var[*])
begin
    nday = dimsizes(var)
    if (nday.eq.1464 .or. nday.eq.1460)then
        ravg = runave_Wrap(var,4,0)
        dailyvar = ravg(1::4)
        delete(ravg)
    end if
    if (nday.eq.365  .or. nday.eq.366)then
        dailyvar = var
    end if
    if (.not.isvar("dailyvar"))then
        printVarSummary(var)
        exit
    end if
    nday = dimsizes(dailyvar)
    ravg = runave_Wrap(dailyvar,5,0)
    pentad = ravg(2::5)
    delete(ravg)
    if (nday .eq. 366)then
        ravg = runave_Wrap(dailyvar(61:),5,0)
        pentad(11) = avg(dailyvar(55:60))
        pentad(12:) = ravg(2::5)
        delete(ravg)
    end if
    pentad@_FillValue = -99.99
    return pentad
end

undef("yrdailytopentadTL")
function yrdailytopentadTL(var[*][*])
begin
    dims = dimsizes(var)
    penvar = new((/73,dims(1)/),"float")
    penvar!0 = "time"
    penvar!1 = "lev"
    ;a = dailytopentad(var(:,0))
    do i =0,  dims(1)-1
        penvar(:,i) = dailytopentad(var(:,i))
    end do
    penvar&time = fspan(1,73,73)
    return penvar
end

undef("yrdailytopentadTLn")
function yrdailytopentadTLn(var)
begin
    ;; assume var t*z*y*x or t*y*x
    dims = dimsizes(var)
    dims(0) = 73
    penvar = new(dims,"float")
    ;a = dailytopentad(var(:,0))
    if(dimsizes(dims).eq.4)
        do k =0,  dims(3)-1
        do j =0,  dims(2)-1
        do i =0,  dims(1)-1
            penvar(:,i,j,k) = dailytopentad(var(:,i,j,k))
        end do
        end do
        end do
        penvar!1 = "lev"
        penvar&lev = var&$var!1$
        penvar!2 = "lat"
        penvar&lat = var&$var!2$
        penvar!3 = "lon"
        penvar&lon = var&$var!3$
    end if
    if(dimsizes(dims).eq.3)
        do j =0,  dims(2)-1
        do i =0,  dims(1)-1
            penvar(:,i,j) = dailytopentad(var(:,i,j))
        end do
        end do
        penvar!1 = "lat"
        penvar&lat = var&$var!1$
        penvar!2 = "lon"
        penvar&lon = var&$var!2$
    end if
    penvar!0 = "time"
    penvar&time = fspan(1,73,73)
    return penvar
end
