function vort2etaLL(vort)
begin
    f = vort
    dims = dimsizes(f)
    pi = 3.14159
    lat = tofloat(f&$f!0$)
    do j = 0, dims(0)-1
        f(j,:) = (2.) *0.000072921   * sin(lat(j)*pi/180)
            ;;    2    omega           sin(lat)
    end do
    ETA = vort
    ETA = vort  + f
    ETA@long_name = "absolute vorticity"
    return ETA
end
function vort2etaTLL(vort)
begin
    f = vort
    dims = dimsizes(f)
    pi = 3.14159
    lat = tofloat(f&$f!1$)
    do j = 0, dims(1)-1
        f(:,j,:) = (2.) *0.000072921   * sin(lat(j)*pi/180)
              ;;    2    omega           sin(lat)
    end do
    ETA = vort
    ETA = vort  + f
    ETA@long_name = "absolute vorticity"
    return ETA
end
function vort2etaTLLL(vort)
begin
    f = vort
    dims = dimsizes(f)
    pi = 3.14159
    lats = tofloat(f&$f!2$)
    do j = 0, dims(2)-1
        f(:,:,j,:) = (2.) *0.000072921   * sin(lats(j)*pi/180)
              ;;    2    omega           sin(lat)
    end do
    ETA = vort
    ETA = vort  + f
    ETA@long_name = "absolute vorticity"
    return ETA
end
