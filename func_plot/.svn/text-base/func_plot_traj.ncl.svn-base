load "func_ty_att.ncl"
undef ("checkRes")
function checkRes(resAtt,resname)
begin
	if(isdefined(resAtt))then
		print("res missing: "+resname)
	end if
	return True
end

undef ("plot_add_traj")
function plot_add_traj(wks:graphic, map:graphic, lat[*]:float,lon[*]:float,mode,color)
begin
   if(mode .eq. "traj")then
      pres                  = True               ; polyline resources
      pres@gsLineThicknessF = 2.0                ; line thickness
      pres@gsLineColor = color
            ;;print("add traj")
      gsn_polyline(wks,map,lon,lat,pres)      ; draw the traj
   end if

; add markers to the trajectories
      mres = True
      mres@gsMarkerIndex  = 16           ; marker style (circle)
      mres@gsMarkerSizeF  = 4.0          ; marker size
      mres@gsMarkerColor  = color      ; maker color
	  if (mode .eq. "traj")then
		  mres@gsMarkerColor = "black"
	  end if
            ;;print("add genesis locate")
	  if (mode .eq. "aveloc")then
        mres@gsMarkerSizeF  = 8.0          ; marker size
		do i = 0, dimsizes(lat)-1
			gsn_polymarker(wks,map,lon(i),lat(i),mres)
		end do
	  else
        gsn_polymarker(wks,map,lon(0),lat(0),mres) ; draw location
	  end if


; create a unique marker to indicate the start of the trajectory
;      first = True
;      first@gsMarkerSizeF  = 8.0         ; marker size
;      first@gsMarkerColor  = "black"     ; marker color
;      gsn_polymarker(wks,map,lon(0),lat(0),first) ; draw start of traj

  return True
end


undef ("plot_traj")
function plot_traj(tydata[*][*][*][8]:float,resTraj)
;;tydata = new((/maxDataYear,maxTy,maxDataNum,8/),"float")    ;; num,YYYY,MM,DD,HH,lat,lon,Vmax
begin
   debug = False
   dims = dimsizes(tydata)
   maxyear = dims(0)
   maxty   = dims(1)
   maxdata = dims(2)

   wks = gsn_open_wks("ps",resTraj@filename)               ; open workstation
   res            = True                          ; map resources
   res@gsnDraw    = False                         ; don't draw
   res@gsnFrame   = False                         ; don't advance frame
   res@vpWidthF   = 0.80                          ; make map bigger
   res@vpHeightF  = 0.80
   res@mpMaxLatF  =  50.                           ; select subregion
   res@mpMinLatF  =  10. 
   res@mpMinLonF  =  90.
   res@mpMaxLonF  = 160.
   res@mpCenterLonF = 180.
   res@mpFillDrawOrder = "PreDraw"
   res@mpFillOn = False
   res@tiMainString = resTraj@title

   colors= (/"red","blue","dark green","grey","magenta"/) ; line color

              ypt = tydata(0,0,:,5)
              xpt = tydata(0,0,:,6)
   map = gsn_csm_map_ce(wks,res)                         ; create map
   tynum   = 0
   twtynum = 0

;;  mode:All TC(include TD) Traj
   if (resTraj@mode .eq. "All TC Traj")then
        do y = 0, maxyear -1
          do ty = 0, maxty -1
            if(.not. ismissing(tydata(y,ty,0,0)))then
              ypt = tydata(y,ty,:,5)
              xpt = tydata(y,ty,:,6)
              addtraj = plot_add_traj(wks,map,ypt,xpt,"traj","dark green")
			  tynum = tynum+1
            end if
          end do
        end do
        draw(wks)
        frame(wks)
   end if
;;  mode:All Ty(exclude TD) Traj
   if (resTraj@mode .eq. "All TY Traj")then
        do y = 0, maxyear -1
          do ty = 0, maxty -1
            if(.not. ismissing(tydata(y,ty,0,0)))then
              ypt = tydata(y,ty,:,5)
              xpt = tydata(y,ty,:,6)
              vmax= tydata(y,ty,:,7)
			  if(.not. isneartaiwan(ypt,xpt))then
	              addtraj = plot_add_traj(wks,map,ypt,xpt,"traj","dark green")
			  end if
		  	tynum = tynum+1
			end if
		  end do
		end do
        do y = 0, maxyear -1
          do ty = 0, maxty -1
            if(.not. ismissing(tydata(y,ty,0,0)))then
              ypt = tydata(y,ty,:,5)
              xpt = tydata(y,ty,:,6)
              vmax= tydata(y,ty,:,7)
              tymaxI = maxind(vmax)
			  if(isneartaiwan(ypt,xpt))then
			    twtynum = twtynum+1
                addtraj = plot_add_traj(wks,map,ypt,xpt,"traj","red")
			    if(debug)then
				    print("TY near TW, source date = "+tydata(y,ty,0,1)+tydata(y,ty,0,2)+tydata(y,ty,0,3))
				    print("TY near TW, source vmax = "+vmax(0))
			    end if
			  end if
            end if
          end do
        end do
		;;tynum = tynum + twtynum
		print(""+twtynum+" / "+tynum)
		txres               = True                 ; text mods desired
		txres@txJust        = "CenterCenter"       ; Default is "CenterCenter"
		txres@txFontHeightF = 0.030 

		if (tynum.ne.twtynum)then
			eqn1 = ":F21::V20:twty:H-60::V-3:____:H-40::V-25::F21:ty"
			stynum   = ":F21:"+tynum
			stwtynum = ":F21:"+twtynum
			gsn_text_ndc(wks,eqn1,.65,.65,txres)
			gsn_text_ndc(wks,":F21:_____",     .80,.66,txres)
			gsn_text_ndc(wks,":F21:=",         .72,.65,txres)
			gsn_text_ndc(wks,stynum,           .80,.63,txres)
			gsn_text_ndc(wks,stwtynum,         .80,.67,txres)
		end if
       	draw(wks)
        frame(wks)
   end if

;;  mode:All Ty genesis location
   if (resTraj@mode .eq. "All TY Gene")then
		delete(ypt)
		delete(xpt)
		tygen = Get_ty_genesis(tydata)
        do y = 0, maxyear -1
          do ty = 0, maxty -1
			ypt = tygen(y,ty,5)
			xpt = tygen(y,ty,6)
            if(.not. ismissing(ypt))then
              addtraj = plot_add_traj(wks,map,ypt,xpt,"genesis","black")
            end if
          end do
        end do
		aveloc = aveGenLocat(tygen)
		addtraj = plot_add_traj(wks,map,aveloc(1:2,0),aveloc(1:2,1),"aveloc","red")
		sloc1 = ":F21:"+"avg. PhS lat/lon = "+aveloc(1,0)+"N / "+aveloc(1,1)+"E"
		sloc2 = ":F21:"+"avg. SCS lat/lon = "+aveloc(2,0)+"N / "+aveloc(2,1)+"E"
		txres               = True                 ; text mods desired
		txres@txJust        = "CenterLeft"       ; Default is "CenterCenter"
		txres@txFontHeightF = 0.020 
		gsn_text_ndc(wks,sloc1,.35,.67,txres)
		gsn_text_ndc(wks,sloc2,.35,.64,txres)

        draw(wks)
        frame(wks)
   end if
    
  return "Done"
end

undef("plot_angle_lat")
function plot_angle_lat(tydata[*][*][*][8],tyturn[*][*][*],res)
begin
    dims = dimsizes(tyturn)

    wks  = gsn_open_wks("ps",res@filename)
    resA            = True
    resA@gsnDraw    = False                         ; don't draw
    resA@gsnFrame   = False       
    resA@trYMinF    =  10.
    resA@trYMaxF    =  40.
    resA@trXMaxF    =  90.
    resA@trXMinF    = -90.
    resL               = True
    resL@gsLineColor   = "red"
    plot = gsn_csm_xy(wks,(/ 0., 0./),(/0.,60./),resA)
    n = 0 
    do y = 0,dims(0)-1
    do ty= 0,dims(1)-1
        DRAWIT = True
        do i = 0, dims(2)-3
            if ((.not.any(ismissing(tyturn(y,ty,i:i+2)))).and.DRAWIT)then
                str = unique_string("polyline")  
                yy =  tydata(y,ty,i:,5)
                xx =  tyturn(y,ty,i:)
                DRAWIT = False
                n = n+1
                ;;print(""+n)
                    plot@$str$ = gsn_add_polyline(wks,plot,xx,yy,resL)
                delete(xx)
                delete(yy)
            end if
        end do
    end do
    end do
    draw(plot)
    frame(wks)
    return True
end

undef("singleTraj")
function singleTraj(traj[*][*],turnangle[*])
local a,wks2,map2
begin
    dims = dimsizes(traj)
    tyn  = floattoint(traj(0,0))
    yyyy = floattoint(traj(0,1))
    mm   = floattoint(traj(0,2))
    dd   = floattoint(traj(0,3))
    filename = "Traj_"+sprinti("%4.4i",yyyy)+"_"+sprinti("%2.2i",mm)+"_"+sprinti("%2.2i",dd)+"_ty"+sprinti("%2.2i",tyn)
    print(filename)
    wks2 = gsn_open_wks("ps",filename)
   res            = True
   res@gsnDraw    = False                         ; don't draw
   res@gsnFrame   = False                         ; don't advance frame
   res@vpWidthF   = 0.80                          ; make map bigger
   res@vpHeightF  = 0.80
   res@mpMaxLatF  =  50.                           ; select subregion
   res@mpMinLatF  =  10.
   res@mpMinLonF  =  90.
   res@mpMaxLonF  = 160.
   res@mpCenterLonF = 180.
   res@mpFillDrawOrder = "PreDraw"
   res@mpFillOn = False
   res@tiMainString = filename

   ypt = traj(:,5)
   xpt = traj(:,6)
   map2 = gsn_csm_map_ce(wks2,res)                         ; create map

   str = unique_string("oneTraj")  
   a = True
   a@$str$ = plot_add_traj(wks2,map2,ypt,xpt,"traj","dark green")
   do i = 0, dims(0)-1
    if (.not.ismissing(ypt(i)))then
            color = "Green"
        if((.not.ismissing(turnangle(i))).and.turnangle(i).ge.30.)then
            color = "Orange"
        end if
        if((.not.ismissing(turnangle(i))).and.turnangle(i).le.-30.)then
            color = "Blue1"
        end if
        if((.not.ismissing(turnangle(i))).and.turnangle(i).ge.60.)then
            color = "OrangeRed"
        end if
        if((.not.ismissing(turnangle(i))).and.turnangle(i).le.-60.)then
            color = "Blue4"
        end if
        if((.not.ismissing(turnangle(i))).and.turnangle(i).ge.90.)then
            color = "Red4"
        end if
        if((.not.ismissing(turnangle(i))).and.turnangle(i).le.-90.)then
            color = "BlueViolet"
        end if
        a@$str$ = plot_add_traj(wks2,map2,ypt(i),xpt(i),"point",color)
    end if
   end do
   dotexample1 = plot_add_traj(wks2,map2,45.,140.,"point","Green")
   dotexample2 = plot_add_traj(wks2,map2,43.,140.,"point","Orange")
   dotexample3 = plot_add_traj(wks2,map2,41.,140.,"point","OrageRed")
   dotexample4 = plot_add_traj(wks2,map2,39.,140.,"point","Red4")
   dotexample5 = plot_add_traj(wks2,map2,43.,145.,"point","Blue1")
   dotexample6 = plot_add_traj(wks2,map2,41.,145.,"point","Blue4")
   dotexample7 = plot_add_traj(wks2,map2,39.,145.,"point","BlueViolet")
   draw(wks2)
   frame(wks2)
   return True

end

undef("plot_yearly_bar")
function plot_yearly_bar(ty,twty,res)
begin
	a = checkRes(res@title,"bar plot title")
	a = checkRes(res@filename,"bar plot filename")
	a = checkRes(res@basin,"basin to get")
	dims    = dimsizes(ty)
	maxyear = dims(0)
	coord   = ty!0
	years   = ty&$coord$
	data1 = new((/2,maxyear/),"float")
    if(res@basin .eq. "All")then
		i = 0
		title = "All "+res@title
	end if
    if(res@basin .eq. "PhS")then
		i = 1
		title = "PhS "+res@title
	end if
    if(res@basin .eq. "SCS")then
		i = 2
		title = "SCS "+res@title
	end if
	data1(0,:) = ty(:,i) ; All,PhS,SCS
	data1(1,:) = twty(:,i) 

	print("All ty/twty correlation: "+escorc(ty(:,0),twty(:,0)))
	print("PhS ty/twty correlation: "+escorc(ty(:,1),twty(:,1)))
	print("SCS ty/twty correlation: "+escorc(ty(:,2),twty(:,2)))


	filename = res@filename
	print("bar fig: "+filename)

	wks   = gsn_open_wks ("ps",filename)           ; open workstation
	Bres                   = True                      ; plot mods desired
	Bres@tiMainString      = ":F21:"+title                     ; add title
	Bres@xyLineThicknesses = (/1.0,2.0/)               ; make 2nd lines thicker
	Bres@xyLineColors      = (/"blue","red"/)          ; change line color
    ;Bres@xyLabelMode            = "Custom"             ; label a line
    ;Bres@xyExplicitLabels       = (/"twty","ty"/)      ; text to use
	;Bres@xyExplicitLegendLabels = (/"twty","ty"/)

    Bres@xyLineLabelFontHeightF = 0.020                ; font height
    Bres@xyLineLabelFontColor   = "black"              ; label color
	Bres@gsnYRefLine			= 0.
	Bres@gsnYRefLineColor		= Bres@xyLineColors(0)
	Bres@gsnYRefLineDashPattern = 1
	Bres2 = Bres
	Bres2@xyLineColors = (/"red","blue"/)
	Bres2@gsnYRefLineColor		= Bres2@xyLineColors(0)

	plot  = gsn_csm_xy (wks,years,data1,Bres) ; create plot

	data2 = dim_standardize_Wrap(data1,0)
	Bres@trYMaxF            =    2.                    ; axis max
	Bres@trYMinF            =   -5.                    ; axis min
	Bres2@trYMaxF           =    5.                    ; axis max
	Bres2@trYMinF           =   -2.                    ; axis min

		corr1 = escorc(data1(0,:),data1(1,:))
        eqn2 = "Cor="+corr1
		Bres@tiMainString = Bres@tiMainString+" "+eqn2
		corr2 = escorc(data2(0,:),data2(1,:))
		print("corr 1 2 : "+corr1+"  "+corr2)
	
	plot  = gsn_csm_xy2 (wks,years,data2(0,:),data2(1,:),Bres,Bres2) ; create plot
	print("year     tynum      twtynum")
	do i = 0, dimsizes(data2(0,:))-1
		print(years(i)+"    "+data2(0,i)+"    "+data2(1,i))
	end do

	return True
end

undef("plot_lon_angle")
function plot_lon_angle(tydata[*][*][*][8],tyturn[*][*][*],res)
begin
    dims = dimsizes(tyturn)

    wks  = gsn_open_wks("ps",res@filename)
    resA            = True
    resA@gsnDraw    = False                         ; don't draw
    resA@gsnFrame   = False       
    resA@trYMaxF    =  160.
    resA@trYMinF    = -160.
    resA@trXMinF    =  120.
    resA@trXMaxF    =  160.
    resL               = True
    resL@gsLineColor   = "red"
    plot = gsn_csm_xy(wks,(/90.,210./),(/0.,0./),resA)
    n = 0 
    do y = 0,dims(0)-1
    do ty= 0,dims(1)-1
        DRAWIT = True
        do i = 0, dims(2)-3
            if ((.not.any(ismissing(tyturn(y,ty,i:i+2)))).and.DRAWIT)then
                str01 = unique_string("polyline")  
                ;;print(str01)
                xx =  tydata(y,ty,i:,6)
                yy =  tyturn(y,ty,i:)
                DRAWIT = False
                n = n+1
                ;;print(""+n)
                    plot@$str01$ = gsn_add_polyline(wks,plot,xx,yy,resL)
                delete(xx)
                delete(yy)
            end if
        end do
    end do
    end do
    draw(plot)
    frame(wks)
    ;;;     ==========================   plot large angle ty traj
    do y = 0,dims(0)-1
    do ty= 0,dims(1)-1
        DRAWIT = True
        do i = 0, dims(2)-3
                xx =  tydata(y,ty,i:,6)
                yy =  tyturn(y,ty,i:)
                if ((.not.any(ismissing(tyturn(y,ty,i:i+2)))) \
                    .and. max(abs(yy)).ge. 30. \
                    .and. DRAWIT)then
                    print("angle > 30 :"+tydata(y,ty,i,0)+" "+tydata(y,ty,i,1)+" "+tydata(y,ty,i,2)+" "+tydata(y,ty,i,3)+" "+tydata(y,ty,i,4)+" "+tydata(y,ty,i,5)+" "+tydata(y,ty,i,6))
                    print("y "+y+"  ty "+ty)
                    a = singleTraj(tydata(y,ty,i:,:),yy)
                    DRAWIT = False
                    ;a = single_lon_ang(tydata(y,ty,:,:),tyturn(y,ty,:))
                    ;a = single_ang_lat(tydata(y,ty,:,:),tyturn(y,ty,:))
                end if
                delete(xx)
                delete(yy)
        end do
    end do
    end do

    return True
end
