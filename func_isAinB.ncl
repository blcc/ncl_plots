function isAinB(a[*],b[*])
begin
    ;;; check every item of a is present in b
    ;;; return logical array in a size
    na = dimsizes(a)
    nb = dimsizes(b)

    isa = a.eq.0
    isa = False
    if(na.ge.nb)then
        do i = 0 , nb -1
            isa = isa.or.(a.eq.b(i))
        end do
    else
        do i = 0 , na -1
            isa(i) = any(a(i).eq.b)
        end do
    end if
    return isa
end
