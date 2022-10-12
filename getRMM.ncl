; This file is part of the MJO_teleconnection module of the MDTF code package (see LICENSE.txt)

;-------------------------------------------------------------------------------------
;Calculate RMM indices for a given  model data
;Projects model data onto reanalysis MJO EOFs
;RMM indices based on Wheeler and Hendon (2004)
;EOF patterns downloaded from: http://cawcr.gov.au/staff/mwheeler/maproom/RMM/eof1and2.htm
;
;Model time domain: 1861-2005.  Only 1950 - 2005 used
;Regional domain: 15N - 15S 
;RMM indices use 250mb zonal wind and NOT 200mb like in reanalysis (200mb not available)
;
;June 22, 2015, updated from RMMcalc.ncl for Libby's data (extra models) March 16, 2016
;Her data has a separate file per level per variable of each model
;Written by Stephanie Henderson
;Jan 10 2018
;Modified by Bohar Singh for MDTF Diagnostic Package 
;Data units: u250,u850: m/s   olr: W/m2
;-------------------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

setfileoption("nc", "Format", "LargeFile")

;--------------Choose Model--------------
;previously used models (already have RMM indices for all of these)
;As Multimodel daignostic the RMM following model are already calculated and will be used in diagnostic  
;Model = "CNRM-CM5"      
;Model = "MRI-CGCM3"     ;don't have omega or olr_250 for this one
;Model = "MIROC5"        ;don't have 250mb OLR for this one
;Model = "IPSL-CM5A-MR"  ;don't have olr_250 for this one
;Model = "GFDL-ESM2M"  
;Model = "GFDL-ESM2G"	 ;noleap
;Model = "GFDL-CM3"       ;noleap.  olr only goes through 2004
;Model = "NorESM1-M"	 ;noleap
;Model = "bcc-csm1-1"
;Model = "CanESM2"	 ;365_day
 Model = "ERA5"

print("---------------------------------------------------------------------------------")

print("Model: " + Model)

;------------ Read in zonal wind and olr --------------------

diri1   = "/global/homes/l/linyaoly/ERA5/reanalysis/"         ;Observed Input data directory
diro    = "/global/homes/l/linyaoly/ERA5/reanalysis/"
; filolr  = "/pscratch/sd/l/linyaoly/ERA5/Unet/reanalysis/ERA5.olr.6hr1979to2022.nc"          ; input file
; filu850 = "/pscratch/sd/l/linyaoly/ERA5/Unet/reanalysis/ERA5.u850.6hr1979to2022.nc"          ; input file
; filu250 = "/pscratch/sd/l/linyaoly/ERA5/Unet/reanalysis/ERA5.u200.6hr1979to2022.nc"         ; input file
filolr  = "/global/homes/l/linyaoly/ERA5/reanalysis/ERA5.olr.day.1978to2022.nc"          ; input file
filu850 = "/global/homes/l/linyaoly/ERA5/reanalysis/ERA5.u850.day.1978to2022.nc"          ; input file
filu250 = "/global/homes/l/linyaoly/ERA5/reanalysis/ERA5.u200.day.1978to2022.nc"         ; input file
; filolr  = "/pscratch/sd/l/linyaoly/ERA5/Unet/reanalysis/ERA5.olr.day.nc"          ; input file
; filu850 = "/pscratch/sd/l/linyaoly/ERA5/Unet/reanalysis/ERA5.u850.day.nc"          ; input file
; filu250 = "/pscratch/sd/l/linyaoly/ERA5/Unet/reanalysis/ERA5.u200.day.nc"         ; input file
olrvar = "olr"
u250var = "u200"
u850var = "u850"

print("Input data path(obs) :" + diri1)
print("Output data path :" + diro)
print("Input data file for OLR :" + filolr)
print("Input data file for U850:" + filu850)
print("Input data file for U250:" + filu250)
print("---------------------------------------------------------------------------------")

; time range to remove seasonal and interanual variability 
ymdStrt = tointeger("1978"+"0903")                         ; start yyyymmdd
ymdLast = tointeger("2022"+"0531") 
yrStrt  = ymdStrt/10000
yrLast  = ymdLast/10000

; EOF analysis strat from year 
ymdstrteof = tointeger("1979"+"0101")
ymdlasteof = tointeger("1979"+"1231")


;---------------------------------------------------------------------------------------
print("Input data starts:" + yrStrt)

;-- Reading U250 data --
    a  = addfile (filu250 , "r")      
  ;  print(a)
     latt = a->lat
    ; print(latt)
     lat_ind = ind(latt .le. 15.5 .and. latt .ge. -15.5)  ;22 latitudes

    if(typeof(latt) .ne. "float")
      lat = tofloat(latt(lat_ind))
    else
      lat = latt(lat_ind)
    end if

    ; print(lat)
    
    lont = a->lon
 
    if(typeof(lont) .ne. "float")
      lon = tofloat(lont)
    else
      lon = lont
    end if
 
   ;--- read in time ---
   temp_time = tointeger(a->time)
  ;  temp_time@units = "hours since 1978-01-01 00:00:00"
   temp_time@units = "days since 1978-01-01 00:00:0.0"
   temp_time@calendar = "proleptic_gregorian"
   print(temp_time(1:10))
   calend = temp_time@calendar
   print(calend)
  ;  print("here" + temp_time@calendar)

  if(calend .ne. "360_day")  ;HadGEM2-CC uses a 360 day calendar
    ; print("here1")
    utc_date = cd_calendar(temp_time, 0)
    ymd     = cd_calendar(temp_time, -2)            ; yyyymmdd
    ; print("ymd:" + ymd)
    yeart   = tointeger(utc_date(:,0))
    montht  = tointeger(utc_date(:,1))
    dayst   = tointeger(utc_date(:,2))
    iStrt   = ind(ymd.eq.ymdStrt)            ; index start
    iLast   = ind(ymd.eq.ymdLast)              ; index last 
    print("CALENDAR DEBUG")   
    printMinMax(ymd,0)
    print("ymdStrt: "+ymdStrt+" iStrt: "+iStrt+" ymdLast: "+ymdLast+" iLast: "+iLast)
    print("Date range from year "+yeart(iStrt)+" to "+yeart(iLast))

    if(typeof(a->$u250var$) .ne. "float")
      u250t = short2flt(a->$u250var$(iStrt:iLast,lat_ind,:))     ;time,lev,lat,lon , 15N - 15S
    else 
      u250t = a->$u250var$(iStrt:iLast,lat_ind,:)
    end if 

    print("type of u250t: " + typeof(u250t))
    u250t@_FillValue = -999

    ; delete(temp_time)
    delete(yeart)
    delete(montht)
    delete(dayst)
  ;***********************************************************
  ; Read user specified time and create required yyyyddd                    
  ;***********************************************************
    time1    = temp_time(iStrt:iLast)             ; time:units = "hours since"
    delete(temp_time)
    TIME    = cd_calendar(time1, 0)             ; type float 

    year    = tointeger( TIME(:,0) )
    month   = tointeger( TIME(:,1) )
    days     = tointeger( TIME(:,2) ) 
    
  ;***********************************************************
    print("Starting date:"+" " +  "Year:"+ year(0) + "  " + "Month:"+ month(0) + "  " + "Day:"+ days(0))
    print("End date:"+"      " +  "Year:"+ year(dimsizes(year)-1) + "  " + "Month:"+ month(dimsizes(year)-1) + " " + "Day:"+ days(dimsizes(year)-1))
    
    ;clean up
    delete(utc_date)
  end if
   
  year@calendar = "gregorian"
  doy = day_of_year(year,month,days)
  time = tointeger(year * 1000) + tointeger(doy)

  ; print("doy: " + min(doy))

  print("First doy: " + min(time) + "  last doy: " + max(time))
  print("Total days: " + dimsizes(year))  ;20440 (4950 DJF) days
  print("-------------------------------------------")
  delete(year)
  delete(month)
  delete(days)
  delete(doy)
print("Reading data for 250mb wind")

  u250 = u250t(:,:,:); Reading data over the all avaliable time period
  u250@_FillValue = default_fillvalue("float")

delete(u250t)
print("U250 data has been read")
;----------- 850mb wind --------------
print("Reading data for 850mb wind")
;-------------------------------------

  b = addfile(filu850,"r")
 
 ;-check lats, lons, and time -
 
  latt2 = b->lat
  lont2 = b->lon
  temp_time = tointeger(b->time)
  ; temp_time@units = "hours since 1978-01-01 00:00:00"
  temp_time@units = "days since 1978-01-01 00:00:0.0"
  temp_time@calendar = "proleptic_gregorian"
  time2 = temp_time(iStrt:iLast)
;print(dimsizes(time2))
;print(dimsizes(time))
;print(time2(0))
;print(time(0))

  if(latt2(0) .eq. latt(0) .and. (dimsizes(latt2) .eq. dimsizes(latt)))
    print("850mb wind: latitudes match")
  else
    print("850mb wind: latitudes don't match")
    exit
  end if
  if(lont2(0) .eq. lont(0).and. (dimsizes(lont2) .eq. dimsizes(lont)))
    print("850mb wind: longitudes match")
  else
    print("850mb wind: longitudes don't match")
    exit
  end if
  if(time2(0) .eq. time1(0).and.(dimsizes(time2) .eq. dimsizes(time1)))
    print("850mb wind: times match")
  else
    print("850mb wind: times don't match")
    exit
  end if
  delete(latt2)
  delete(lont2)
  delete(time2)
 
 ;-read in data-
  if(typeof(b->$u850var$) .ne. "float")
    u850 = short2flt(b->$u850var$(iStrt:iLast,lat_ind,:))     ;time,lat,lon , 15N - 15S
  else 
    u850 = b->$u850var$(iStrt:iLast,lat_ind,:)    ;time,lat,lon , 15N - 15S
  end if 

  u850@_FillValue = default_fillvalue("float")
print("U850 data has been read")
print("-------------------------------------------")

   delete(b)

;-- olr --
print("Now reading OLR data")

  
  b = addfile(filolr,"r")

 ;-check lats, lons, and time -
 
  latt2 = b->lat
  lont2 = b->lon
  temp_time = tointeger(b->time)
  ; temp_time@units = "hours since 1978-01-01 00:00:00"
  temp_time@units = "days since 1978-01-01 00:00:0.0"
  temp_time@calendar = "proleptic_gregorian"
  time2 = temp_time(iStrt:iLast)


  if(latt2(0) .eq. latt(0) .and. (dimsizes(latt2) .eq. dimsizes(latt)))
    print("olr: latitudes match")
  else
    print("olr: latitudes don't match")
    exit
  end if
  if(lont2(0) .eq. lont(0).and.(dimsizes(lont2) .eq. dimsizes(lont)))
    print("olr: longitudes match")
  else
    print("olr: longitudes don't match")
    exit
  end if
  
 if(time2(0) .eq. time1(0).and. (dimsizes(time2) .eq. dimsizes(time1)))
    print("olr: times match")
  else
    print("olr: times don't match")
  end if  

  delete(latt2)
  delete(lont2)
  delete(latt)
  ; delete(lont)
  delete(time2)

  
 ;-read in data-
  if(typeof(b->$olrvar$) .ne. "float")
    olr = short2flt(b->$olrvar$(iStrt:iLast,lat_ind,:))     ;time,lat,lon , 15N - 15S
  else 
    olr = b->$olrvar$(iStrt:iLast,lat_ind,:)
  end if 
    olr = -olr/3600.
    olr@_FillValue = default_fillvalue("float")

  ;clean up
  delete(b)
  delete(lat_ind)
print("OLR data has been read")
print("olr "+olr(0,0,0))
print("-------------------------------------------")

;-------------------------------------------------------------------------
;------ remove mean of previous 120 days ------
print("remove interannual variability")
;-approximately removes interannual variability
;-not calculated for first 120 days (which are not used so will be tossed)
;-------------------------------------------------------------------------
;first find the interannual mean.  start on day 120
 
 u200_int = new(dimsizes(u250),typeof(u250),u250@_FillValue)
 u850_int = new(dimsizes(u850),typeof(u850),u850@_FillValue)
 olr_int =  new(dimsizes(olr),typeof(olr),olr@_FillValue)

;printVarSummary(u200_int)

;prescribe an interannual mean of zero for first 120 days, which will get tossed
 
 u200_int(:120-1,:,:) = 0.0 
 u850_int(:120-1,:,:) = 0.0 
 olr_int(:120-1,:,:) = 0.0 

 u200_int@_FillValue = default_fillvalue("float")
 u850_int@_FillValue = default_fillvalue("float")
 olr_int@_FillValue = default_fillvalue("float")

 do t = 120,dimsizes(time)-1
  l = t - 120   ;beginning index of the previous 120 days
    u200_int(t,:,:) = dim_avg_n(u250(l:t,:,:),0)
    u850_int(t,:,:) = dim_avg_n(u850(l:t,:,:),0)
    olr_int(t,:,:) = dim_avg_n(olr(l:t,:,:),0)  
    delete(l)
 end do

;remove interannual
 u200_p = u250 - u200_int
 u850_p = u850 - u850_int
 olr_p = olr - olr_int

 u200_p@_FillValue = default_fillvalue("float")
 u850_p@_FillValue = default_fillvalue("float")
 olr_p@_FillValue = default_fillvalue("float")
 print("Interannual variability has been removed")


;clean up
 delete(u200_int)
 delete(u850_int)
 delete(olr_int)
 delete(u250)
 delete(u850)
 delete(olr)

 u200_p!0 = "time"
 u200_p!1 = "lat" 
 u200_p!2 = "lon" 

 u850_p!0 = "time"
 u850_p!1 = "lat" 
 u850_p!2 = "lon"

 olr_p!0 = "time"
 olr_p!1 = "lat" 
 olr_p!2 = "lon"


print("-------------------------------------------")
;------ remove seasonal cycle ------
print("remove seasonal cycle")
print("calculate raw climatology then smooth it")

if(any(ismissing(u200_p)))
  print("u200 missing values")
else
  print("no missing values")
end if

if(any(ismissing(u850_p)))
  print("u850 missing values")
else
  print("no missing values")
end if

if(any(ismissing(olr_p)))
  print("olr missing values")
else
  print("no missing values")
end if

nsize = dimsizes(u200_p)   ; time, lat, lon
print("nsize: " + nsize)
tmp = time1(120:nsize(0)-1)
delete(time1)
time1 = tmp 
delete(tmp)
print("time0: " + time1(0))
print("timeend: " + max(time1))
tmp = time(120:nsize(0)-1)
delete(time)
time = tmp 
delete(tmp)
print("time0: " + time(0))
print("timeend: " + max(time))
; delete(time)
u200_p1 = u200_p(120:nsize(0)-1,:,:)
delete(u200_p)
u850_p1 = u850_p(120:nsize(0)-1,:,:)
delete(u850_p)
olr_p1  = olr_p(120:nsize(0)-1,:,:)
delete(olr_p)

;200mb wind
 rawclim = clmDayTLL(u200_p1,time)   ; (120:dimsizes(time)-1)
 print("check rawclim wind: "+rawclim(0:4,0,0))
 print("rawclim size must be 366: " + dimsizes(rawclim))
 smoothclim = smthClmDayTLL(rawclim,3)  ;remove first 3 harmonics
 print("smoothclim for u200: "+dimsizes(smoothclim)) 

;print(rawclim(:,10,10) + "     " + smoothclim(:,10,10))

;remove seasonal cycle
 u200_anom = calcDayAnomTLL(u200_p1,time,smoothclim) ;remove seasonal 
 delete(rawclim) 
 delete(smoothclim)
 delete(u200_p1)

;850mb wind
 rawclim = clmDayTLL(u850_p1,time)
 smoothclim = smthClmDayTLL(rawclim,3)

 u850_anom = calcDayAnomTLL(u850_p1,time,smoothclim)
 delete(rawclim) 
 delete(smoothclim)
 delete(u850_p1)

;OLR
 rawclim = clmDayTLL(olr_p1,time)
 smoothclim = smthClmDayTLL(rawclim,3)
 print("check rawclim olr: "+rawclim(0:4,0,0))
 print("rawclim size: " + dimsizes(rawclim))
 
 olr_anom = calcDayAnomTLL(olr_p1,time,smoothclim)
 print("olr_anom size: " + dimsizes(olr_anom))
 
 delete(rawclim) 
 delete(smoothclim)
 delete(olr_p1)
print("seasonal cycle has been removed")
print("-------------------------------------------")

;print(olr_p(150:200,10,10) + "    " + olr_anom(150:200,10,10))

;- average the data latitudes
print("averaging the data over 15S-15N")
u200_anom@_FillValue = default_fillvalue("float")
u850_anom@_FillValue = default_fillvalue("float")
olr_anom@_FillValue = default_fillvalue("float")

u200 = dim_avg_n(u200_anom,1)  ;time,lon
u850 = dim_avg_n(u850_anom,1)
olr = dim_avg_n(olr_anom,1)
print("Done with latitudinal averaging")
print("-------------------------------------------")
;clean up
delete(u200_anom)
delete(u850_anom)
delete(olr_anom)

;------ compute the temporal variance at each lon ------
var_olr  = dim_variance_n_Wrap( olr, 0)     ; (lon)
var_u850 = dim_variance_n_Wrap(u850, 0)
var_u200 = dim_variance_n_Wrap(u200, 0)
if(any(ismissing(var_olr)))
  print("olr missing values Again")
else
  print("no missing values for var_olr")
end if

if(any(ismissing(var_u850)))
  print("olr missing values Again")
else
  print("no missing values for var_olr")
end if

if(any(ismissing(var_u200)))
  print("olr missing values Again")
else
  print("no missing values for var_olr")
end if

;************************************************
; Compute the zonal mean of the temporal variance
;************************************************
zavg_var_olr  = dim_avg_n_Wrap( var_olr , 0)     
zavg_var_u850 = dim_avg_n_Wrap( var_u850, 0)
zavg_var_u200 = dim_avg_n_Wrap( var_u200, 0)

;************************************************
; Normalize by sqrt(avg_var*)
;************************************************
olr   =  olr/sqrt(zavg_var_olr )          ; (time,lon)
u850  = u850/sqrt(zavg_var_u850)
u200  = u200/sqrt(zavg_var_u200)

print("norm values:")
print(sqrt(zavg_var_olr )  )
print(sqrt(zavg_var_u850))
print(sqrt(zavg_var_u200))


dim_lon = 180      ; 2deg
dl = dim_lon * 2   ;360
dl2 = dim_lon * 3  ;540

;- combine all three fields into one array (3*lon, time)
RMM_field = new((/dl2,dimsizes(time)/),typeof(olr),olr@_FillValue)  ; (120:dimsizes(time)-1)
lon_field = new((/dl2/),typeof(lont),lont@_FillValue)
do ml=0,dim_lon-1
  RMM_field(ml,:) = (/  olr(:,ml) /)
  RMM_field(ml+dim_lon,:) = (/ u850(:,ml) /)
  RMM_field(ml+dl,:) = (/ u200(:,ml) /)
  lon_field(ml) = (/  lont(ml) /)
  lon_field(ml+dim_lon) = (/  lont(ml) /)
  lon_field(ml+dl) = (/  lont(ml) /)
end do

RMM_field@_FillValue = -999

;printVarSummary(RMM_field)
;clean up
delete(olr)
delete(u850)
delete(u200)
delete(dl)

printVarSummary(lont)
print(lont(0))
;************************************************
; Compute **combined** EOF; Sign of EOF is arbitrary
;************************************************
neof = 2 
delete(ymd)
; calculate EOF for the selected time
printVarSummary(time1)

; time1@units = "days since 1978-01-01 00:00:0.0"
; time1@calendar = "proleptic_gregorian"
ymd = cd_calendar(time1, -2)
istrteof = ind(ymd .eq. ymdstrteof)
ilasteof = ind(ymd .eq. ymdlasteof)

printVarSummary(ymd)

print("time1 " + time1(0))
print("istrteof: " + istrteof)
print("ilasteof: " + ilasteof)
print("eof starts: " + ymd(istrteof))
print("eof ends: " + ymd(ilasteof))

tmp = time1(istrteof:ilasteof)
delete(time1)
time1 = tmp
delete(tmp)
; time1@units = "days since 1978-01-01 00:00:0.0"
; time1@calendar = "proleptic_gregorian"
tmp = RMM_field(:,istrteof:ilasteof)
delete(RMM_field)
RMM_field = tmp 
delete(tmp)

printVarSummary(RMM_field)

eof_RMM_field    = eofunc_Wrap(RMM_field, neof, False)      ; (neof,3*mlon)

;************************************************
; Change the Sign of EOF to be consistent with WH04
;************************************************

; find the position of the maximum OLR
imax_olr_eof1 = maxind(eof_RMM_field(0,0:dim_lon-1))
lonmax_eof1 = lont(imax_olr_eof1)

if ( (lonmax_eof1 .ge. 100) .and. (lonmax_eof1 .lt. 180)) then
  eof_RMM_field(0,:) = -eof_RMM_field(0,:) 
end if 

imax_olr_eof2 = maxind(eof_RMM_field(1,0:dim_lon-1))
lonmax_eof2 = lont(imax_olr_eof2)

if ( (lonmax_eof2 .ge. 120) .and. (lonmax_eof2 .lt. 220)) then
  eof_RMM_field(1,:) = -eof_RMM_field(1,:) 
end if 


print("==============")
printVarSummary(eof_RMM_field)
printMinMax(eof_RMM_field, True)

eof_ts_RMM_field = eofunc_ts_Wrap(RMM_field,eof_RMM_field,False)   ; (neof,ntim)
delete(RMM_field)
print("==============")                                  
printVarSummary(eof_ts_RMM_field)
printMinMax(eof_ts_RMM_field, True)

pc = transpose(eof_ts_RMM_field)     ;(ntim,neof)
delete(eof_ts_RMM_field)
; ei_vectors = transpose(eof_RMM_field)
; ei_vectors@_FillValue=-999
; print(ei_vectors(0,:))
ei_values = eof_RMM_field@eval 

eofout =  transpose(eof_RMM_field)
header = (/"longitude,    EOF1,    EOF2"/)
hlist = [/header/]

print("eofout size: " + dimsizes(eofout))

alist = [/lon_field,eofout(:,0),eofout(:,1)/]
eilist = [/ei_values/]

write_table(diro+Model+"_RMM_EOF_ERA5only"+ymdstrteof+"to"+ymdlasteof+".txt", "w", hlist, "%s ")
write_table(diro+Model+"_RMM_EOF_ERA5only"+ymdstrteof+"to"+ymdlasteof+".txt", "a", alist, "%07i, %7.7f, %7.7f")
; write_table(diro+Model+"_RMM_EOF_ERA5only"+ymdStrt+".txt", "a", eilist, "%7.7f, %7.7f")

delete(eof_RMM_field)
print("eigenvalues:")
print(ei_values)

;- normalize by EOF-calculated std dev
 do m = 0,1
   pc(:,m) = pc(:,m)/sqrt(ei_values(m))
 end do

;--change due to Eric's suggestion --
;instead, normalize by the standard deviation of each model so that std dev of each RMM index is unity
do m = 0,1
   pc(:,m) = dim_standardize_n(pc(:,m),0,0)
end do
print("==============")                                  
printVarSummary(pc)
printMinMax(pc, True)
; ;print("PC1: " + pc(:,0) + "            PC2: " + pc(:,1))  ;20454 values


;------ calculate MJO phase ------
print("calculate MJO phase")
phase = (/ 202.5, 247.5, 292.5, 337.5,  22.5,  67.5, 112.5, 157.5 /) ; phi: Matt Wheeler's f77 code

  nPhase      = dimsizes(phase)
  phaBnd      = 180.0/nPhase
  angBnd      = new( (/2,nPhase/), "float")
  angBnd(0,:) = phase-phaBnd 
  angBnd(1,:) = phase+phaBnd
 ;print("phase: "+phase+"   angBnd: "+angBnd(0,:)+" - "+angBnd(1,:))

  r2d         = 180./(4.*atan(1.0))
  ang         = atan2(pc(:,1),pc(:,0))*r2d     ; phase space 
  nn          = ind(ang.lt.0)
  ang(nn)     = ang(nn) + 360          ; make 0 to 360
  ;print(ang)

delete(phase)
delete(nPhase)
delete(r2d)
delete(nn)
delete(phaBnd)

phase = new((/dimsizes(ang)/),integer)

do t = 0,dimsizes(ang)-1
 do i=0,7
  if(ang(t) .ge. angBnd(0,i) .and. ang(t) .le. angBnd(1,i))
    phase(t) = i+1
  end if
 end do
end do

;print("Angle: " + ang + "     Phase:" + phase)

delete(ang)
delete(angBnd)

;------calculate amplitude------
print("calculate amplitude")
amp = sqrt(pc(:,0)^2. + pc(:,1)^2.) 

;print("RMM1: " + pc(:,0) + "   RMM2: " + pc(:,1) + "   phase: " + phase + "   amplitude: " + amp)

;------write output file------

header = (/"Time,    RMM1,    RMM2,    phase,    amplitude"/)
hlist = [/header/]

alist = [/time1,pc(:,0),pc(:,1),phase,amp/]


 write_table(diro+Model+"_RMM_ERA5only"+ymdstrteof+"to"+ymdlasteof+".txt", "w", hlist, "%s ")
 write_table(diro+Model+"_RMM_ERA5only"+ymdstrteof+"to"+ymdlasteof+".txt", "a", alist, "%07i, %7.7f, %7.7f, %01i, %7.7f")

print("Data is saved into txt format at:"+diro)
print("In following file:")
print(Model+"_RMM_ERA5only"+ymdstrteof+"to"+ymdlasteof+".txt")
print("Process Completed")
print("==================================================")
