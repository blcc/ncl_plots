;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
function isneartaiwan (lat[*]:float, lon[*]:float)
local i,j,k,neartaiwan,ntarea
begin
  debug = False
  neartaiwan = False
  ntarea = (/             37681,37682,37683,37684,37685,       \
                    36960,36961,36962,36963,36964,36965,       \
              36239,36240,36241,36242,36243,36244,36245,36246, \
              35519,35520,35521,35522,35523,35524,35525,35526, \
        34798,34799,34800,34801,34802,34803,34804,34805,34806, \
  34077,34078,34079,34080,34081,34082,34083,34084,34085,       \
  33357,33358,33359,33360,33361,33362,33363,33364,33365,       \
  32637,32638,32639,32640,32641,32642,32643,32644,             \
        31918,31919,31920,31921,31922,31923,31924,             \
        31198,31199,31200,31201,31202,31203,                   \
              30479,30480,30481,30482,30483,                   \
                    29760,29761,29762/)
  i = floattoint((lon/0.5))
  j = floattoint((lat/0.5))
  k = j*720+i       ; 720 = 360/0.5 
  
;;  print("k = "+k+" get1Dindex = "+get1Dindex(ntarea,k))

  do l = 0, dimsizes(k)-1
    if(.not. ismissing(get1Dindex(ntarea,k(l))))then
      neartaiwan = True
	  if(debug)then
	    print("TY near TW, source lat/lon:"+lon(0)+" / "+lat(0))
	  end if
      break
    end if
  end do

  return(neartaiwan)

end
