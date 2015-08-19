undef("isnotmissing")
function isnotmissing(dd[*])
begin
	istrue = False
	dim = dimsizes(dd)
	do i = 0, dim-1
		if(.not. ismissing(dd(i)))then
			istrue = True
			break
		end if
	end do
	return istrue
end
