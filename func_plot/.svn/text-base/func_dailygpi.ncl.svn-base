
function cal_dailygpi(year[1]:integer)
begin
    ;; 850hPa absolute vorticity
    vn = "VOR"
    vn@lev = 850 ;; hPa
        VOR = read_ERAdataDaily(year,vn)
        f = VOR
        dims = dimsizes(f)
        f = 0
        pi = 3.14159
        do j = 0, dims(1)-1
            ;f(:,j,:) = (2.)*(2*pi/86400.) * sin(f&lat(j)*pi/180)
            f(:,j,:) = (2.) *0.000072921   * sin(f&lat(j)*pi/180)
                  ;;    2    omega          sin(lat)
        end do
        ETA = VOR
        ETA = VOR + f

    ;; 700hPa relative humidity
    vn = "RH"
    vn@lev = 700
        RH = read_ERAdataDaily(year,vn)

    ;; Tropical cyclone potential intensity
    Vpot = read_vpotDaily(year)

    ;; 200-850 hPa wind shear magnitude
    vn = "U"
    vn@lev = 200
            U = read_ERAdataDaily(year,vn)
    vn@lev = 850
            Ulo = read_ERAdataDaily(year,vn)
        U = U - Ulo
    
    vn = "V"
    vn@lev = 200
            V = read_ERAdataDaily(year,vn)
    vn@lev = 850
            Vlo = read_ERAdataDaily(year,vn)
        V = V - Vlo

    Vshear = V
    Vshear = sqrt(U*U+V*V)

    ;; GPI terms
    cETA = ETA
    cRH  = RH
    cVpot= Vpot
    cVshear = Vshear

    cETA = abs(100000*cETA)^1.5
    cRH  = (cRH/50)^3
    cVpot= (cVpot/70)^3
    cVshear = (1+0.1*cVshear)^(-2)

    ;; GPI(12,lat,lon)
    GPI = Vpot  ;; for coordinates
    GPI = GPI@_FillValue
    GPI = cETA * cRH * cVpot * cVshear

    GPI@long_name = "Genesis potential index"

    return GPI
end
