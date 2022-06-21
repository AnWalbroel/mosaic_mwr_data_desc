;+
;**********************
PRO PL_LEVEL1B_HATPRO,$
;**********************
;INPUT:
date=date,$                   ;STRING (yymmdd)
station=station,$             ;3 character STRING
ceilo_type=ceilo_type         ;either 'igm' or 'dwd' currently
; $Id: $
; Abstract:
; * reads quality controlled daily time series of 
;   iwv, lwp from HATPRO netcdf level1b data and plots them
; Authors:
; U. Loehnert (loehnert@meteo.uni-koeln.de)
; Date:
; 2008-10-27
; Dependencies:
; - get_next_date.pro
; - write_hatpro_level0b_nc.pro
; - ceilometer reading routines - if ceilometer available
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

dummy  = -9.e+33
iwv_min = 0.
iwv_max = 30.
IF station EQ 'amf' THEN iwv_max = 50.
IF station EQ 'ufs' THEN iwv_max = 20.
lwp_min = -50.
lwp_max = 600.

date_m = STRMID(date, 0, 4)
date_y = STRMID(date, 0, 2)
date_m_x = STRMID(date, 2, 2)
date_d = STRMID(date, 4, 2)
file_path  = '/data/data_hatpro/'+station + '/data/level1b/'+ date_m+'/'
out_path   = '/data/data_hatpro/'+station + '/plots/level1b/'+ date_m+'/'

GET_NEXT_DATE, '20'+date, yyyy_next, mm_next, dd_next
date_d_next = dd_next
date_y_next = STRMID(yyyy_next, 2, 2)
date_m_x_next = mm_next
date_m_next = STRMID(yyyy_next, 2, 2) + mm_next
date_next = STRMID(yyyy_next, 2, 2) + mm_next + dd_next

plot_file  = out_path + date + '_h_'+ station + '_l1b.ps'
 
n_in_ceilo = 0
in_name_ceilo = ''

IF station EQ 'amf' THEN BEGIN

 ceilo_path = '/data/data_hatpro/amf_other/ceilometer/'+ date_m+'/'
 spawn,'ls '+ceilo_path + '*' + date + '*', in_name_ceilo

 n_in_ceilo = N_ELEMENTS(in_name_ceilo)
 IF in_name_ceilo(0) EQ '' THEN n_in_ceilo = 0

ENDIF ELSE IF station EQ 'ufs' THEN BEGIN

 IF ceilo_type EQ 'dwd' THEN BEGIN

  ceilo_path = '/data/data_hatpro/ufs_ceilo/'+ date_m+'/'
  ceilo_path_next = '/data/data_hatpro/ufs_ceilo/'+ date_m_next+'/'
 
  dn = FLOAT(date_y+date_m_x+date_d)
  dn_next = FLOAT(date_y_next+date_m_x_next+date_d_next)

  IF dn LT 70417. THEN BEGIN
   spawn,'ls '+ceilo_path + '*'+date_y+date_m_x+date_d+'*.T02', in_name_ceilo1
   spawn,'ls '+ceilo_path_next + '*'+date_y_next+date_m_x_next+date_d_next+'*.T02', in_name_ceilo2
   in_name_ceilo = [in_name_ceilo1, in_name_ceilo2]
  ENDIF ELSE IF dn GE 70417. THEN BEGIN
   spawn,'ls '+ceilo_path + '*'+date_y+'-'+date_m_x+'-'+date_d+'*.t02', in_name_ceilo1
   spawn,'ls '+ceilo_path_next + '*'+date_y_next+'-'+date_m_x_next+'-'+date_d_next+'*.t02', in_name_ceilo2
   in_name_ceilo = [in_name_ceilo1, in_name_ceilo2]
  ENDIF

  ij = WHERE(in_name_ceilo NE '')
  IF ij(0) NE -1 THEN in_name_ceilo = in_name_ceilo(ij)

  help = SIZE(in_name_ceilo)
  n_in_ceilo = help(1)
  IF ij(0) EQ -1 THEN n_in_ceilo = 0

 ENDIF ELSE IF ceilo_type EQ 'igm' THEN BEGIN

  ceilo_path = '/projekt/TOSCA/ceilo_igm/data/'+ date_m+'/'
  spawn,'ls '+ceilo_path + date + '*.nc', in_name_ceilo
  n_in_ceilo = 1
 ENDIF
ENDIF ELSE IF station EQ 'sel' THEN BEGIN 
 ceilo_path = '/data/TR32/D2/data/ceilo/level0b/'+date_m+'/'
 in_name_ceilo = FILE_SEARCH(ceilo_path + date + '*.nc')
 n_in_ceilo = 1
ENDIF ELSE IF station EQ 'jue' THEN BEGIN
 ceilo_path = '/data/TR32/D2/data/ceilo/jenoptik/prelim/2009/' + date_m + '/'
 in_name_ceilo = FILE_SEARCH(ceilo_path + '*' + date + '*.nc')
 n_in_ceilo = 1
ENDIF ELSE n_in_ceilo = 0

PRINT,'Total number of CEILO files: ',n_in_ceilo, in_name_ceilo

n = 100000
lcb_all = FLTARR(n)
time_c = FLTARR(n)
start_ceilo = 0
IF in_name_ceilo(0) EQ '' THEN n_in_ceilo = 0

FOR i_f = 0,n_in_ceilo-1 DO BEGIN

  IF station EQ 'ufs' AND ceilo_type EQ 'dwd' THEN BEGIN
   READ_CEILO_UFS, 0, in_name_ceilo(i_f), time_ceilo, lcb
   n_ceilo = LONG(N_ELEMENTS(time_ceilo))

   IF n_ceilo GT 0 THEN BEGIN
    IF i_f EQ 0 THEN BEGIN
     time_c(start_ceilo:start_ceilo+n_ceilo-1l) = time_ceilo-1. ; stimmt die -1 immer ???
     lcb_all(start_ceilo:start_ceilo+n_ceilo-1l) = lcb
     start_ceilo = start_ceilo + n_ceilo
    ENDIF
    IF i_f EQ 1 THEN BEGIN
     time_c(start_ceilo:start_ceilo+n_ceilo-1l) = time_ceilo-1.+24.
     lcb_all(start_ceilo:start_ceilo+n_ceilo-1l) = lcb
     start_ceilo = start_ceilo + n_ceilo
    ENDIF
   ENDIF 

  ENDIF ELSE IF station EQ 'jue' THEN BEGIN
   s = '/data/TR32/D2/data/ceilo/jenoptik/idl'
   if strpos(!PATH, s ) lt 0 then !PATH = !PATH + ':' + s
   IF n_elements(in_name_ceilo) LE 1 THEN BEGIN
    read_jo_ceilo_nc, filename = in_name_ceilo(0), time = time_julian, first_cbh = lcb_c
    time_help = time_julian-floor(time_julian)
    time_c = time_help*24.
   ENDIF ELSE BEGIN
    FOR i=0, n_elements(in_name_ceilo)-1 DO BEGIN
     read_jo_ceilo_nc, filename = in_name_ceilo(i), time = time_julian, first_cbh = lcb_c
     IF i GT 0 THEN BEGIN
      time_julian_ = time_julian
      lcb_c_ = lcb_c 
      k=0
      j=0
      i=0
      time_julian_help = dblarr(n_elements(time_julian_) + n_elements(time_julian))
      lcb_c_help = dblarr(n_elements(time_julian_) + n_elements(time_julian))
      WHILE j LE n_elements(time_julian_)-1 OR k LE n_elements(time_julian)-1 DO BEGIN
       time_julian_help(i) = min([time_julian_(j), time_julian(k)])
       IF time_julian_(j) GT time_julian(k) THEN BEGIN
        lcb_c_help(i) = lcb_c(k)
        k=k+1
       ENDIF ELSE IF time_julian_(j) LT time_julian(k) THEN BEGIN
        lcb_c_help(i) = lcb_c_(j)
        j=j+1
       ENDIF ELSE IF time_julian_(j) EQ time_julian(k) THEN BEGIN 
        lcb_c_help(i) = lcb_c_(j)
        j=j+1
        k=k+1
       ENDIF
       i=i+1
      ENDWHILE
      IF j LT n_elements(time_julian_) THEN BEGIN 
       time_julian_help(i:i+n_elements(time_julian_)-j-1)=time_julian_(j+1:n_elements(time_julian_)-1)
       lcb_c_help(i:i+n_elements(time_julian_)-j-1)=lcb_c_(j+1:n_elements(time_julian_)-1)
      ENDIF ELSE IF k LT n_elements(time_julian) THEN BEGIN 
       time_julian_help(i:i+n_elements(time_julian_)-k-1)=time_julian(k+1:n_elements(time_julian)-1)
       lcb_c_help(i:i+n_elements(time_julian)-k-1)=lcb_c(k+1:n_elements(time_julian)-1)
      ENDIF
      index_ceilo = where(time_julian_help NE 0.0)
      time_julian = time_julian_help(index_ceilo)
      lcb_c = lcb_c_help(index_ceilo)
     ENDIF
    ENDFOR
   ENDELSE
   lcb_all = lcb_c
   n_ceilo = N_ELEMENTS(time_c)
   start_ceilo = n_ceilo
  ENDIF ELSE IF (station EQ 'ufs' AND ceilo_type EQ 'igm') OR station EQ 'sel' THEN BEGIN 

   READ_CT25K_LEVEL0B_NC, in_name_ceilo(0),$
                          time_c, height_ceilo, detection_status, status_ceilo, lcb_all, scb_all,$
                          tcb_all, vertical_vis, max_bscat, laser_pulse_energy, laser_temperature, receiver_sensitivity,$
                          window_contamination, tilt_angle, background_light, sum_backscatter,$
                          backscatter, measurement_parameters, status_string, lat_ceilo, lon_ceilo, alt_ceilo,$
                          time_julian=time_julian, base_time=base_time, fill_gaps=fill_gaps, novalue=novalue, verbose=verbose

   i_n = WHERE(vertical_vis GE 0. AND vertical_vis LT 200.)
   IF i_n(0) NE -1 THEN lcb_all(i_n) = 1.

   n_ceilo = N_ELEMENTS(time_c)
   start_ceilo = n_ceilo

  ENDIF

  IF station EQ 'amf' THEN BEGIN
   GET_CEILO_AMF, 0, in_name_ceilo(i_f), time_ceilo, lcb, dts

   n_ceilo = N_ELEMENTS(time_ceilo)

   IF n_ceilo GT 0 THEN BEGIN
    time_c(start_ceilo:start_ceilo+n_ceilo-1) = time_ceilo/3600.
    lcb_all(start_ceilo:start_ceilo+n_ceilo-1) = lcb
    start_ceilo = start_ceilo + n_ceilo
   ENDIF

  ENDIF

ENDFOR

IF start_ceilo GT 1 THEN BEGIN
 lcb_all = lcb_all(0:start_ceilo-1)
 time_c = time_c(0:start_ceilo-1)
ENDIF

IF station EQ 'amm' THEN BEGIN 
 ceilo_path='/data/data_hatpro/amm/other/ceilo/level0b/'+ date_m+'/'
 filename=ceilo_path+date+'_ct25k_AE.CT25K_Od_l0b.nc'
 READ_CT25K_LEVEL0B_NC,filename,$
    time_ceilo, height_ceilo, detection_status, status_flag, first_cbh, second_cbh,$
    third_cbh, laser_pulse_energy,laser_temperature, receiver_sensitivity,$
    window_contamination, tilt_angle, background_light, sum_backscatter,$
    backscatter, measurement_parameters, status_string, lat, lon, alt
 n_ceilo=N_ELEMENTS(time_ceilo)
   
 IF n_ceilo GT 0 THEN BEGIN

    time_c(start_ceilo:start_ceilo+n_ceilo-1) = time_ceilo
    lcb_all(start_ceilo:start_ceilo+n_ceilo-1) = first_cbh
    start_ceilo = start_ceilo + n_ceilo
 ENDIF

IF start_ceilo GT 1 THEN BEGIN
 lcb_all = lcb_all(0:start_ceilo-1)
 time_c = time_c(0:start_ceilo-1)
ENDIF
ENDIF

; get COPS GPS measurements
xx=FLTARR(5)
time_gps=FLTARR(96)
iwv_gps=FLTARR(96)
gpsindex=0
IF station eq 'hor' or station eq 'ach' then begin
   file_gps='/data/data_hatpro/amf/data/gps/'+station+'/gps_'+station+'_'+date+'.dat'
   yes=FILE_TEST(file_gps)
   IF yes eq 1 then begin
   OPENR,unit_gps,file_gps,/get_lun
   gpsindex=1
   i=0
   WHILE not eof(unit_gps) DO BEGIN
     READF,unit_gps,xx
     time_gps(i) = xx(0)/60.
     iwv_gps(i)  = xx(1)
     i=i+1
   ENDWHILE
   ENDIF
ENDIF

;----------------------------------------------------------------------
; Get information about campaign/station
;----------------------------------------------------------------------
GET_STATION_INFO, 0, station, start_time, stop_time, stat_name, latitude,$
                  longitude, altitude, instrument, comment_ret, comment_tech

;----------------------------------------------------------------------
; READ data in netcdf format
;----------------------------------------------------------------------
filename = file_path + date + '_h_'+station+'_l1b.nc'
READ_HATPRO_LEVEL1_NC, filename, algo, commentx,$
                       time, iwv, lwp, flag, temp,pres, relh, az, el

comment = comment_ret

IF station EQ 'sel' OR station EQ 'jue' OR station EQ 'ind' THEN BEGIN
 ang_b = 89.7
 ang_e = 90.3
 i_zen=WHERE(el GT ang_b AND el LT ang_e)
 time = time(i_zen)
 flag = flag(i_zen)
 az = az(i_zen)
 el = el(i_zen)
 iwv = iwv(i_zen)
 lwp =lwp(i_zen)
 temp = temp(i_zen)
 pres = pres(i_zen)
 relh = relh(i_zen)
ENDIF



;----------------------------------------------------------
; plot settings
;----------------------------------------------------------
!p.font =0
!p.multi=[0,2,12]
!P.Charsize=1.6
!P.Charthick=1.5
!P.thick = 2.
SET_PLOT, 'PS'
DEVICE,/color,/landscape,filename=plot_file
TEK_COLOR

x_min = 0.17
x_max = 0.9
xchar = 1.6

;----------------------------------------------------------
; plot integrated water vapor
;----------------------------------------------------------
index = WHERE(iwv GT iwv_min AND iwv LT iwv_max, count)

IF count GT 0 THEN BEGIN
  ymin = Min(iwv(index))
  ymax = Max(iwv(index))
ENDIF ELSE BEGIN
  ymin = dummy - 20
  ymax = dummy + 20
ENDELSE
dy = ymax - ymin
ymin = ymin - dy/10.
ymax = ymax + dy/10.

print,ymin,ymax
;.........................................integrated water vapor
PLOT, time, iwv, yrange=[ymin, ymax], ystyle = 1,$
  xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]',$
  xcharsize = 1.2, ytickv = 2,$ ;ytitle='IWV / kg m!e-2!n',$
  position=[x_min,0.1,x_max,0.26667], xthick=2, xticklen=-0.06, /NODATA

IF index(0) NE -1 THEN BEGIN
 OPLOT, time(index), iwv(index), psym = 3

 m = MEAN(iwv(index))
 str = STRING(m,FORMAT='(F8.2)') + ' +/- '+$
       STRING(STDDEV(iwv(index)),FORMAT='(F8.2)') + ' kg m!e-2!n'
 XYOUTS,0.20,0.245,str,/NORMAL,size=0.8, color = 2
ENDIF

IF gpsindex eq 1 THEN BEGIN
  OPLOT,time_gps,iwv_gps,psym=1,color=3
  XYOUTS,0.75,0.245,'green: IWV GPS',/normal,size=0.8,color=3
ENDIF

XYOUTS,0.1,0.35, 'LWP',alignment=0.5,charsize=1.0,/normal
XYOUTS,0.1,0.32,'[g m!e-2!n]',alignment=0.5,charsize=1.0,/normal
XYOUTS,0.1,0.20,'IWV',alignment=0.5,charsize=1.0,/normal
XYOUTS,0.1,0.17, '[kg m!e-2!n]',alignment=0.5,charsize=1.0,/normal
XYOUTS,0.5,0.93, 'HATPRO level1b: '+stat_name+date,alignment=0.5,charsize=1.3,/NORMAl
XYOUTS,0.05,0.01,comment,/NORMAL,charsize=0.8,color=2

;----------------------------------------------------------
; plot liquid water path
;----------------------------------------------------------
index = WHERE(lwp GT lwp_min AND lwp LT lwp_max, count)
IF index(0) NE -1 THEN lwp(index) = lwp(index)*1000.

IF count GT 0 THEN BEGIN
  ymin = lwp_min
  ymax = lwp_max
ENDIF ELSE BEGIN
  ymin = dummy - 20
  ymax = dummy + 20
ENDELSE
dy = ymax - ymin
ymin = ymin-dy/10.
ymax = ymax + dy/10.

PLOT, time, lwp, yrange=[ymin, ymax], ystyle = 1,$
   xrange = [0., 24.], xstyle=1,$ ; ytitle = 'LWP / g m!e-2!n',$
   xcharsize = 0.00001, ytickv = 2, xticklen=-0.06, $
   position=[x_min,0.2667,x_max,0.43333], /NODATA
OPLOT, [0,24.],[0.,0.], color=2,linestyle=1

IF index(0) NE -1 THEN BEGIN
 OPLOT, time(index), lwp(index), psym = 3
 m = MEAN(lwp(index))
 str = STRING(m,FORMAT='(F8.2)') + ' +/- '+$
       STRING(STDDEV(lwp(index)),FORMAT='(F8.2)') + ' g m!e-2!n'
 XYOUTS,0.20,0.41,str,/NORMAL,size=0.8, color = 2
ENDIF




;----------------------------------------------------------
; plot ceilo
;----------------------------------------------------------
IF start_ceilo GT 1 THEN BEGIN

IF ceilo_type EQ 'dwd' THEN ytitc = 'LCB(DWD) [km]'
IF ceilo_type EQ 'igm' THEN ytitc = 'LCB(IGM) [km]'
IF ceilo_type EQ 'jop' THEN ytitc = 'LCB(JOP) [km]'

 i99 = WHERE(lcb_all GE 0.)
 IF i99(0) NE -1 THEN lcb_all(i99) = lcb_all(i99)/1000.
 PLOT, time_c, lcb_all, yrange=[0., 7.], ystyle = 1,$
       xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]',$
       xcharsize = 0.00001, ytickv = 2, ytitle = ytitc, $
       position=[x_min,0.4333,x_max,0.6], xticklen=-0.06, /NODATA, /NOERASE

 OPLOT, time_c, lcb_all, color = 4, thick = 2, psym = 4, symsize = 0.3
ENDIF


;----------------------------------------------------------
; plot flags
;----------------------------------------------------------
time_res = 60.
xx = 5*time_res/3600.
x = [-(xx), -(xx), (xx), (xx)]
y = [-1, 1, 1,-1]
USERSYM, x, y, /fill

y_min = 0.6
;........................................iwv stddev threshold test
index = WHERE(flag AND 16,count)
print,'Number of samples which did not pass IWV stddev threshold =',count
PLOT,[0,0],[24,1],/NODATA, yrange=[0.8,1.2], ystyle = 1,$
     xrange = [0., 24.], xstyle=1,xcharsize = 0.0000001, ycharsize=0.00001,$
     position=[x_min,y_min+0.08,x_max,y_min+0.1],yticks=1
IF count GT 0 THEN BEGIN
   IF count EQ 1 THEN BEGIN
     count = 2
     index = [index(0),index(0)]
   ENDIF
   h = REPLICATE(1,count)
   OPLOT,time(index),h,psym=8,thick=2,color = 12
ENDIF
XYOUTS,0.1,y_min+0.09,'IWV STD thres',charsize=0.7,alignment=0.5,/Normal

;........................................threshold test
index = WHERE(flag AND 8,count)
Print,'Number of times where IWV and/or LWP threshold is not passed=',count
PLOT,[0,0],[24,1],/NODATA, yrange=[0.8,1.2], ystyle = 1,$
     xrange = [0., 24.], xstyle=1,xcharsize = 0.0000001, ycharsize=0.00001,$
     position=[x_min,y_min+0.06,x_max,y_min+0.08],yticks=1
IF count GT 0 THEN BEGIN
   IF count EQ 1 THEN BEGIN
     count = 2
     index = [index(0),index(0)]
   ENDIF
   h = REPLICATE(1,count)
   OPLOT,time(index),h,psym=8,thick=2,color = 8
ENDIF
XYOUTS,0.1,y_min+0.07,'IWV/LWP thres',charsize=0.7,alignment=0.5,/Normal

;........................................TB thresholds
index = WHERE(flag AND 4,count)
print,'Number of samples where TB threshold not met =',count
PLOT,[0,0],[24,1],/NODATA, yrange=[0.8,1.2], ystyle = 1,$
     xrange = [0., 24.], xstyle=1,xcharsize = 0.0000001, ycharsize=0.00001,$
     position=[x_min,y_min+0.04,x_max,y_min+0.06],yticks=1
IF count GT 0 THEN BEGIN
   IF count EQ 1 THEN BEGIN
     count = 2
     index = [index(0),index(0)]
   ENDIF
   h = REPLICATE(1,count)
   OPLOT,time(index),h,psym=8,thick=2,color = 4
ENDIF
XYOUTS,0.1,y_min+0.05,'TB threshold',alignment=0.5,charsize=0.7,/Normal

;........................................visual inspection
index = WHERE(flag AND 2,count)
print,'Number of samples filtered out by visual inspection =',count
PLOT,[0,0],[24,1],/NODATA, yrange=[0.8,1.2], ystyle = 1,$
     xrange = [0., 24.], xstyle=1,xcharsize = 0.0000001, ycharsize=0.00001,$
     position=[x_min,y_min+0.02,x_max,y_min+0.04],yticks=1
IF count GT 0 THEN BEGIN
   IF count EQ 1 THEN BEGIN
     count = 2
     index = [index(0), index(0)]
   ENDIF
   h = REPLICATE(1,count)
   OPLOT,time(index),h,thick=2,psym=8,color = 2
ENDIF
XYOUTS,0.1,y_min+0.03,'Visual filter',alignment=0.5,charsize=0.7,/Normal

;........................................rain flag

;defekter Regenfilter
;day  = FIX(STRMID(date,4,2))
;IF day GE 20 THEN BEGIN
;XYOUTS,0.5,y_min+0.004,'DEFEKTER REGENFILTER',alignment=0.5,charsize=0.8,color = 3,/Normal
;ENDIF
;defekter Regenfilter
;month= FIX(STRMID(date,2,2))
;if month NE 07 then begin

index = WHERE(flag AND 1,count)
print,'Number of samples with HATPRO rain flag =',count
PLOT,[0,0],[24,1],/NODATA, yrange=[0.8,1.2], ystyle = 1,$
     xrange = [0., 24.], xstyle=1,xcharsize = 0.0000001, ycharsize=0.00001,$
     position=[x_min,y_min,x_max,y_min+0.02],yticks=1, xticklen=0.08
IF count GT 0 THEN BEGIN
IF count EQ 1 THEN BEGIN
     count = 2
     index = [index(0),index(0)]
   ENDIF
   h = REPLICATE(1,count)
   OPLOT,time(index),h,psym=8,thick=2,color = 3
ENDIF
XYOUTS,0.1,y_min+0.01,'HATPRO Rain',alignment=0.5,charsize=0.7,/Normal


;endif
;XYOUTS,0.5,y_min+0.004,'DEFEKTER REGENFILTER',alignment=0.5,charsize=0.8,color = 3,/Normal 
;------------------------------------------------------------------------------------------
; plot environmental temperature and relative humidity
;------------------------------------------------------------------------------------------
index = WHERE(relh NE dummy, count)
ymin = Min(relh(index))
IF ymin LT 0 THEN ymin =0.
ymax = 99.
IF ymin GT ymax THEN ymin = 0.
dy = ymax - ymin
PLOT,time,relh,Yrange=[ymin,ymax],ystyle=1,$
   xcharsize=0.00001,ytickv=1,ytitle='RH / % ',$
   Xrange=[0.,24.], Xstyle=1, ycharsize = 0.7,$
   position=[x_min,0.7,x_max,0.8],xticklen=0.08,/NODATA
OPLOT, time,relh, psym = 3

index = WHERE(temp NE dummy, count)
ymin = Min(temp(index))
ymax = Max(temp(index))
dy = ymax - ymin
ymin = ymin - dy/10.
ymax = ymax + dy/10.
PLOT,time,temp,Yrange=[ymin,ymax],ystyle=1,$
   xcharsize=0.00001,ytickv=2,ytitle='T [K]',ycharsize = 0.7,$
   Xrange=[0.,24.], Xstyle=1, Xtitle='Time (UTC) / h',$
   position=[x_min,0.8,x_max,0.9],xticklen=0.08,/NODATA
OPLOT, time, temp, psym = 3


DEVICE,/CLOSE


END

