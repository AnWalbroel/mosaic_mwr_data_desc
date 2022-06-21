;+
;**********************
PRO PL_LEVEL2A_HATPRO,$
;**********************
;keywords
Date=date,$              ; YYMMDD to be plotted
Station=station,$        ; three letter code
int_algo=int_algo,$              ; IWV/LWP algorithm files without 'IWV_' and 'LWV_' in folders '/home/hatpro/retrievals/IWV/' and '.../LWP/'
q_algo=q_algo,$          ; q algorithm files
T_algo=T_algo,$          ; T algorithm files (elevation scanning!)
IWV_max=IWV_max, $       ; maximum IWV to be plotted, default=80., if IWV_max=0 autmatically adapted to stddev during the day, if IWV_max<0 range=[mean-|IWV_max/2|,mean+|IWV_max/2|]
ceilo_type=ceilo_type, $ ; either 'igm' or 'dwd' currently
rd_file=rd_file,$        ; 'wvl': use WVL/OLC files (default), 'brt': use BRT files
verbose=verbose           
; $Id:$
; Abstract: 
; * plots level 2a HATPRO/WORA data --> only intended for first quicklooks
; * Temp., q&IWV, RH&LWP, HATPRO Az-IWV, Radar&LWP, Ceilometer data ...
; Authors:
; U. Loehnert
; Date:
; 2010-08
; Dependencies:
; - 
; Changes
;-

!Except=2
if ~keyword_set(verbose) then verbose=0
if ~keyword_set(IWV_max) then IWV_max=80.

IF N_ELEMENTS(rd_file) EQ 0 THEN rd_file = 'wvl'
IF N_ELEMENTS(rd_file) NE 0 THEN BEGIN
 IF rd_file NE 'brt' THEN rd_file = 'wvl'
ENDIF

IF N_ELEMENTS(ceilo_type) EQ 0 THEN ceilo_type = 'xxx'

;***set plotting thresholds
t_thres_min = 150.
t_thres_max = 330.
q_thres_min = 0.
q_thres_max = 20.
rh_thres_min = -5.
rh_thres_max = 105.

date_m= STRMID(date, 0, 4)
date_d = STRMID(date, 4, 2)
date_y = STRMID(date, 0, 2)
date_m_x = STRMID(date, 2, 2)

raw_path = '/data/data_hatpro/'+station+'/data/raw/'+date_m+'/'+date+'/'
plot_path = '/data/data_hatpro/'+station+'/plots/level1a/'+date_m+'/'

GET_PREV_DATE, '20'+date, yyyy_prev, mm_prev, dd_prev
date_m_prev = STRMID(yyyy_prev, 2, 2) + mm_prev
date_prev = STRMID(yyyy_prev, 2, 2) + mm_prev + dd_prev
raw_path_prev = '/data/data_hatpro/'+station+'/data/raw/'+date_m_prev+'/'+date_prev+'/'

GET_NEXT_DATE, '20'+date, yyyy_next, mm_next, dd_next
date_d_next = dd_next 
date_y_next = STRMID(yyyy_next, 2, 2)
date_m_x_next = mm_next 
date_m_next = STRMID(yyyy_next, 2, 2) + mm_next
date_next = STRMID(yyyy_next, 2, 2) + mm_next + dd_next

;****brightness temperatures
IF rd_file EQ 'wvl' THEN BEGIN

; search and count WVL files
 in_name_wvl = FILE_SEARCH(raw_path + '*.wvl')
 in_name_wvl_prev = FILE_SEARCH(raw_path_prev + '*.wvl')

 IF in_name_wvl_prev(0) NE '' THEN in_name_wvl = [in_name_wvl_prev, in_name_wvl]
 n_in_wvl = N_ELEMENTS(in_name_wvl)
 IF in_name_wvl(0) EQ '' THEN n_in_wvl = 0 

 if verbose then PRINT,'Total number of WVL files being checked: ',n_in_wvl

ENDIF ELSE BEGIN

 ; search and count BRT files
 in_name_brt = FILE_SEARCH(raw_path + '*.brt')
 in_name_brt_prev = FILE_SEARCH(raw_path_prev + '*.brt')

 IF in_name_brt_prev(0) NE '' THEN in_name_brt = [in_name_brt_prev, in_name_brt]

 n_in_wvl = N_ELEMENTS(in_name_brt)
 IF in_name_brt(0) EQ '' THEN n_in_wvl = 0 

 if verbose then PRINT,'Total number of BRT files being checked: ',n_in_wvl

ENDELSE

;****BLB files
in_name_blb = FILE_SEARCH(raw_path + '*.blb')
in_name_blb_prev = FILE_SEARCH(raw_path_prev + '*.blb')

IF in_name_blb_prev(0) NE '' THEN in_name_blb = [in_name_blb_prev, in_name_blb]

n_in_blb = N_ELEMENTS(in_name_blb)
IF in_name_blb(0) EQ '' THEN n_in_blb = 0 

if verbose then PRINT,'Total number of BLB files being checked: ',n_in_blb

;****ceilometer data
n_in_ceilo = 0
in_name_ceilo = ''

IF station EQ 'amf' THEN BEGIN

 ceilo_path = '/data/data_hatpro/amf_other/ceilometer/'+ date_m+'/'
 spawn,'ls '+ceilo_path + '*' + date + '*', in_name_ceilo, err_msg

  ii = WHERE(in_name_ceilo NE '')
  IF ii(0) EQ -1 THEN BEGIN
   n_in_ceilo = 0
  ENDIF ELSE BEGIN
   in_name_ceilo = in_name_ceilo(ii)
   n_in_ceilo = N_ELEMENTS(in_name_ceilo)
  ENDELSE

 print, 'Total number of CEILO files: ',n_in_ceilo

ENDIF 

IF station EQ 'ufs' and ceilo_type EQ 'dwd' THEN BEGIN
 ceilo_path = '/data/data_hatpro/ufs_ceilo/'+ date_m+'/'
 ceilo_path_next = '/data/data_hatpro/ufs_ceilo/'+ date_m_next+'/'

 IF dn LT 70417. THEN BEGIN
  spawn,'ls '+ceilo_path + '*'+date_y+date_m_x+date_d+'*.T02', in_name_ceilo1
  spawn,'ls '+ceilo_path_next + '*'+date_y_next+date_m_x_next+date_d_next+'*.T02', in_name_ceilo2
  in_name_ceilo = [in_name_ceilo1, in_name_ceilo2]
 ENDIF ELSE IF dn GE 70417. THEN BEGIN
  spawn,'ls '+ceilo_path + '*'+date_y+'-'+date_m_x+'-'+date_d+'*.t02', in_name_ceilo1
  spawn,'ls '+ceilo_path_next + '*'+date_y_next+'-'+date_m_x_next+'-'+date_d_next+'*.t02', in_name_ceilo2
  in_name_ceilo = [in_name_ceilo1, in_name_ceilo2]
 ENDIF

 help = SIZE(in_name_ceilo)
 n_in_ceilo = help(1)
 IF in_name_ceilo(0) EQ '' THEN n_in_ceilo = 0
 PRINT,'Total number of CEILO files: ',n_in_ceilo
ENDIF

IF station EQ 'ufs' and ceilo_type EQ 'igm' THEN BEGIN
 ceilo_path = '/projekt/TOSCA/ceilo_igm/data/'+ date_m+'/'
 spawn,'ls '+ceilo_path + date + '*.nc', in_name_ceilo
 n_in_ceilo = 1
ENDIF

;****met data
in_name_met = FILE_SEARCH(raw_path + '*.met')
in_name_met_prev = FILE_SEARCH(raw_path_prev + '*.met')

IF in_name_met_prev(0) NE '' THEN in_name_met = [in_name_met_prev, in_name_met]
n_in_met = N_ELEMENTS(in_name_met)
IF in_name_met(0) EQ '' THEN n_in_met = 0

if verbose then PRINT,'Total number of MET files being checked: ',n_in_met

;****initialize arrays
n = 200000
n_z = 43
time_w = LONARR(n)
time_woz = LONARR(n)
time_m = LONARR(n)
time_c = FLTARR(n)
time_b = FLTARR(n)

iwv_all = FLTARR(n)
lwp_all = FLTARR(n)
iwv_all_oz = FLTARR(n)
lwp_all_oz = FLTARR(n)
el_all_oz = FLTARR(n)
az_all_oz = FLTARR(n)
scan_all = FLTARR(n)
rain_all_wvl = BYTARR(n)
rain_all_wvl_oz = BYTARR(n)

rain_all = BYTARR(n)
rain_all_blb = BYTARR(n)
met_all  = FLTARR(3, n)
temp_all = FLTARR(n_z, n)
hum_all =  FLTARR(n_z, n)
lcb_all = FLTARR(n)
tb_blb_all = REPLICATE(-999., 14, 6, n)

;****calculate IWV & LWP
start_wvl = 0l
start_wvl_oz = 0l
start_met = 0l
start_ceilo  = 0l
start_blb = 0l

;***process zenith iwv/lwp observations
coeff_dat =  '/home/hatpro/retrievals/IWV/IWV_'+int_algo
GET_COEFF_INT, 0, coeff_dat, lin_qua, angle, f, offsetv, c_lv, c_qv
coeff_dat =  '/home/hatpro/retrievals/LWP/LWP_'+int_algo
GET_COEFF_INT, 0, coeff_dat, lin_qua, angle, f, offsetl, c_ll, c_ql

i_scan = 0l

FOR i_f = 0, n_in_wvl-1 DO BEGIN

 IF rd_file EQ 'wvl' THEN BEGIN
  GET_SPEC, 0, in_name_wvl(i_f), n_wvl, time_wvl, rain_wvl, freq, TB_wvl, az_wvl, el_wvl
 ENDIF ELSE BEGIN
  GET_BRT, 0, in_name_brt(i_f), n_wvl, time_wvl, rain_wvl, freq_brt, TB_brt, az_wvl, el_wvl
  IF n_wvl GT 0 THEN BEGIN
   TB_wvl = TB_brt(0:6, *)
   freq = freq_brt
  ENDIF
 ENDELSE

 IF n_wvl EQ 0 THEN GOTO, NEXTFILE

 i_zen = WHERE(el_wvl GT 89.7 AND el_wvl LT 90.3)
 n_wvl_z = LONG(N_ELEMENTS(i_zen))

 IF i_zen(0) EQ -1 THEN GOTO, NEXTFILE 

 iwv = REPLICATE(offsetv, n_wvl_z)
 lwp = REPLICATE(offsetl, n_wvl_z)
 TB_new = TB_wvl

 nn = 7 ; default
  
;**iwv

 IF N_ELEMENTS(f) EQ 3 AND station EQ 'jue' THEN BEGIN ; 22.24, 23.04, 25.44, 26.24 GHz channels disturbed
  nn = 3
  TB_new = FLTARR(nn, n_wvl)
  TB_new(0, *) = TB_wvl(2, *)      ; 23.84 GHz 
  TB_new(1:2, *) = TB_wvl(5:6, *)  ; 27.84, 31.4 GHz
 ENDIF

 FOR il = 0l, n_wvl_z-1l DO iwv(il) = TOTAL(c_lv*TB_new(0:nn-1, i_zen(il))) + TOTAL(c_qv*(TB_new(0:nn-1, i_zen(il))^2)) + iwv(il)

;**lwp

 IF N_ELEMENTS(f) EQ 3 AND station EQ 'jue' THEN BEGIN ; 22.24, 23.04, 25.44, 26.24 GHz channels disturbed
  nn = 3
  TB_new = FLTARR(nn, n_wvl)
  TB_new(0, *) = TB_wvl(2, *)      ; 23.84 GHz 
  TB_new(1:2, *) = TB_wvl(5:6, *)  ; 27.84, 31.4 GHz
 ENDIF

 FOR il = 0l, n_wvl_z-1l DO lwp(il) = TOTAL(c_ll*TB_new(0:nn-1, i_zen(il))) + TOTAL(c_ql*(TB_new(0:nn-1, i_zen(il))^2)) + lwp(il)

 lwp = lwp*1000.

 IF n_wvl_z GT 0 THEN BEGIN
  time_w(start_wvl:start_wvl+n_wvl_z-1l) = time_wvl(i_zen)
  iwv_all(start_wvl:start_wvl+n_wvl_z-1l) = iwv
  lwp_all(start_wvl:start_wvl+n_wvl_z-1l) = lwp
  rain_all_wvl(start_wvl:start_wvl+n_wvl_z-1l) = rain_wvl(i_zen)
  start_wvl = start_wvl + n_wvl_z
 ENDIF
NEXTFILE:
ENDFOR ; loop over n_in_wvl
;**transform time_w to Julian time
time_w = time_w/86400d
time_w = time_w + JULDAY(1,1,2001,0,0,0)

;***process off-zenith iwv/lwp observations

IF station EQ 'jue' THEN BEGIN


 FOR i_f = 0, n_in_wvl-1 DO BEGIN
  GET_BRT, 0, in_name_brt(i_f), n_wvl, time_wvl, rain_wvl, freq_brt, TB_brt, az_wvl, el_wvl
  IF n_wvl GT 0 THEN BEGIN
   TB_wvl = TB_brt(0:6, *)
   freq = freq_brt
  ENDIF

  IF n_wvl EQ 0 THEN GOTO, NEXTFILE_OZ

  i_off_zen = WHERE(el_wvl LT 89.7)
  n_wvl_oz = LONG(N_ELEMENTS(i_off_zen))

  i_off_zen = WHERE(el_wvl LE 89.7)
  n_wvl_oz = LONG(N_ELEMENTS(i_off_zen))

  IF i_off_zen(0) EQ -1 THEN GOTO, NEXTFILE_OZ 

  IF N_ELEMENTS(f) EQ 3 THEN BEGIN ; 22.24, 23.04, 25.44, 26.24 GHz channels disturbed
   nn = 3
   TB_new = FLTARR(nn, n_wvl_oz)
   TB_new(0, *) = TB_wvl(2, *)      ; 23.84 GHz 
   TB_new(1:2, *) = TB_wvl(5:6, *)  ; 27.84, 31.4 GHz
  ENDIF

  co_path = '/home/hatpro/retrievals/'
  xx = FINDGEN(100)
  ix = WHERE(xx GT 44 AND xx LT 70)
  xx = xx(ix)
  cx = STRING(xx, format = '(i2)')

  coeff_dat = file_search(co_path+'IWV/IWV_db_0112_V'+cx+'.RET')
  RET_IWVLWP_MA, coeff_dat, TB_new(*, i_off_zen), el_wvl(i_off_zen), iwv_oz, freq
  iwv_oz = iwv_oz*sin(!dpi*el_wvl(i_off_zen)/180d) ; airmass correction

  coeff_dat = file_search(co_path+'LWP/LWP_db_0112_V'+cx+'.RET')
  RET_IWVLWP_MA, coeff_dat, TB_new(*, i_off_zen), el_wvl(i_off_zen), lwp_oz, freq
  lwp_oz = lwp_oz*sin(!dpi*el_wvl(i_off_zen)/180d) ; airmass correction

  IF n_wvl_oz GT 0 THEN BEGIN
   time_woz(start_wvl_oz:start_wvl_oz+n_wvl_oz-1l) = time_wvl(i_off_zen)
   iwv_all_oz(start_wvl_oz:start_wvl_oz+n_wvl_oz-1l) = iwv_oz
   lwp_all_oz(start_wvl_oz:start_wvl_oz+n_wvl_oz-1l) = lwp_oz
   el_all_oz(start_wvl_oz:start_wvl_oz+n_wvl_oz-1l) = el_wvl(i_off_zen)
   az_all_oz(start_wvl_oz:start_wvl_oz+n_wvl_oz-1l) = az_wvl(i_off_zen)
   rain_all_wvl_oz(start_wvl_oz:start_wvl_oz+n_wvl_oz-1l) = rain_wvl(i_off_zen)
   scan_all(start_wvl_oz:start_wvl_oz+n_wvl_oz-1l) = i_scan
   start_wvl_oz = start_wvl_oz + n_wvl_oz
   i_scan = i_scan + 1
  ENDIF
 NEXTFILE_OZ:
 ENDFOR ; loop over n_in_wvl_oz

;**transform time_w to Julian time
 time_woz = time_woz/86400d
 time_woz = time_woz + JULDAY(1,1,2001,0,0,0)
ENDIF ; station Juelich

;****met files
FOR i_f = 0,n_in_met-1 DO BEGIN
    GET_MET,in_name_met(i_f),time_met,rain,temp,pres,humi,n_met
    IF n_met GT 1 THEN BEGIN
      n_met= LONG(N_ELEMENTS(time_met))
      time_m(start_met:start_met+n_met-1l) = time_met
      met_all(0,start_met:start_met+n_met-1l) = temp
      met_all(1,start_met:start_met+n_met-1l) = pres
      met_all(2,start_met:start_met+n_met-1l) = humi
      rain_all(start_met:start_met+n_met-1l) = rain
      start_met = start_met + n_met
    ENDIF
ENDFOR

;**transform time_m to Julian time
time_m = time_m/86400d
time_m = time_m + JULDAY(1,1,2001,0,0,0)

;****calculate T&q profiles
FOR i_f = 0, n_in_blb-1 DO BEGIN

 GET_BLB, 0, in_name_blb(i_f), n_blb, time_blb, rain_blb, freq_blb, n_ang, tb_blb, ang_blb
 IF n_blb GT 0 THEN BEGIN
  n_blb= LONG(N_ELEMENTS(time_blb))
  time_b(start_blb:start_blb+n_blb-1l) = time_blb
  tb_blb_all(*, *, start_blb:start_blb+n_blb-1l) = tb_blb   
  rain_all_blb(start_blb:start_blb+n_blb-1l) = rain_blb
  start_blb = start_blb + n_blb
 ENDIF 
ENDFOR

;**transform time_b to Julian time
time_b = time_b/86400d
time_b = time_b + JULDAY(1,1,2001,0,0,0)
n_b = start_blb

;****calculate T-Profile (level2a)

;**read temp. coefficients 
coeff_file_T = '/home/hatpro/retrievals/TProfile_ELE/TEL_'+T_algo
GET_COEFF_TEL_3Z_4BL, 0, coeff_file_T, angles, f, z, offset, coeff_TEL, coeff_Tgr, Tgr_go=0
;**read humidity coefficients
print, 'Warning: all humidity channels considered!!
coeff_file_q =  '/home/hatpro/retrievals/HProfile/HUM_'+q_algo
GET_COEFF_ZEN, 0, coeff_file_q, 1, angle, f, z, offset_q, c_l, c_q

;**HATPRO TB bias correction
IF station EQ 'jue' OR station eq 'sel' THEN BEGIN
 restore, '/data/data_hatpro/jue/tbcor/TB_meas_minus_calc_mm_jue_09.sav'
ENDIF ELSE IF station EQ 'pa0' THEN BEGIN
 restore, '/data/data_hatpro/pa0/tbcor/TB_meas_minus_calc_pa0.sav'
 IF algo EQ 'pay_0112_V03.RET' THEN Tgr_go = 1
 IF Tgr_go EQ 1 THEN BEGIN
  restore, '/home/hatpro/pay_analysis/2m_temp_pay/pay_temp_2006_2009.sav'
 ENDIF
ENDIF ELSE BEGIN
 print, 'warning: no TB-bias correction in V-band implemented!' 
ENDELSE

FOR i = 0l, n_b-1 DO BEGIN

 Tgr = 0.

 IF station EQ 'jue' OR station eq 'sel' THEN BEGIN;für sel gibt es keine tb-Korrekturwerte, deshalb werden die von jue benutzt
;***tb-bias correction
  IF (FLOAT(date_m_x) LE 5. AND FLOAT(date_y) LE 9.) OR (FLOAT(date_y) EQ 8.) then begin ;es gibt nur für 04/09 und 08/09 TB-Korrektur-Werte
   tb_cor = TB_meas_minus_calc_mm(3,*,*)
  ENDIF ELSE IF (date_m_x GE 6. AND date_y GT 9.) then begin
   tb_cor = TB_meas_minus_calc_mm(7,*,*)
  ENDIF ELSE BEGIN
   tb_cor = REPLICATE(0d, 7, n_ang)
   print, 'Warning: no V-Band bias correction implemented!' 
  ENDELSE
  tb_cor=reform(tb_cor)
  TB_algo_T = [tb_blb_all(7:13, 0, i) - tb_cor(*, 0)]
  TB_algo_q = tb_blb_all(0:6, 0, i)

  FOR k = 1, n_ang-1 DO TB_algo_T = [TB_algo_T, tb_blb_all(10:13, k, i) - tb_cor(3:6, k)]
;***retrieval application
  FOR j = 0l, n_z-1 DO temp_all(j, i) = coeff_TEL(j, *)#TB_algo_T + offset(j)
  FOR j = 0l, n_z-1 DO hum_all(j, i) = c_l(j, *)#TB_algo_q + (c_q(j, *)#TB_algo_q^2) + offset_q(j)

 ENDIF ELSE IF station EQ 'pa0' THEN BEGIN ; this is the way the bias correction should be carried out in future!
;***tb-bias correction

  ny = N_ELEMENTS(TB_meas_minus_calc(0, *, 0))
  nz = N_ELEMENTS(TB_meas_minus_calc(0, 0, *))
  tb_cor = REPLICATE(0., ny, nz)

  n_crit_new = N_ELEMENTS(time_crit_new)
  FOR ii = 0, n_crit_new-2 DO BEGIN
   IF time_b(i) GT time_crit_new(ii) AND time_b(i) LE time_crit_new(ii+1) THEN tb_cor = TB_meas_minus_calc(ii, *, *)
  ENDFOR
  tb_cor=reform(tb_cor)
;**with bias correction
  TB_algo_T = [tb_blb_all(7:13, 0, i) - tb_cor(*, 0)]
  FOR k = 1, n_ang-1 DO TB_algo = [TB_algo_T, tb_blb_all(10:13, k, i) - tb_cor(3:6, k)]

;**without bias correction
;  TB_algo_T = [tb_blb_all(7:13, 0, i)]
;  FOR k = 1, n_ang-1 DO TB_algo = [TB_algo_T, tb_blb_all(10:13, k, i)]
   TB_algo_q = tb_blb_all(0:6, 0, i)

;***retrieval application
  icheck = WHERE(TB_algo LE 0.)
  IF icheck(0) EQ -1 THEN BEGIN
   FOR j = 0l, n_z-1 DO temp_all(i, j) = coeff_TEL(j, *)#TB_algo_T + offset(j)
  ENDIF

 ENDIF 

ENDFOR


;****cloud radar data ...

;**** ceilometer data
IF in_name_ceilo(0) NE '' THEN BEGIN
 FOR i_f = 0, n_in_ceilo-1 DO BEGIN

  IF station EQ 'ufs' AND ceilo_type EQ 'dwd' THEN BEGIN
   READ_CEILO_UFS, 0, in_name_ceilo(i_f), time_ceilo, lcb
   n_ceilo = N_ELEMENTS(time_ceilo)

   IF n_ceilo GT 0 THEN BEGIN
    IF i_f EQ 0 THEN BEGIN
     time_c(start_ceilo:start_ceilo+n_ceilo-1) = time_ceilo-1.
     lcb_all(start_ceilo:start_ceilo+n_ceilo-1) = lcb
     start_ceilo = start_ceilo + n_ceilo
    ENDIF
    IF i_f EQ 1 THEN BEGIN
     time_c(start_ceilo:start_ceilo+n_ceilo-1) = time_ceilo-1.+24.
     lcb_all(start_ceilo:start_ceilo+n_ceilo-1) = lcb
     start_ceilo = start_ceilo + n_ceilo
    ENDIF
   ENDIF

  ENDIF

  IF station EQ 'ufs' AND ceilo_type EQ 'igm' THEN BEGIN
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
ENDIF

IF start_wvl EQ 0 THEN start_wvl = 2
IF start_blb EQ 0 THEN start_blb = 2
IF start_met EQ 0 THEN start_met = 2
IF start_wvl_oz EQ 0 THEN start_wvl_oz = 2

time_w = time_w(0l:start_wvl-1l)
IF start_wvl_oz GT 1 THEN time_woz = time_woz(0l:start_wvl_oz-1l)
time_m = time_m(0l:start_met-1l)
time_b = time_b(0l:start_blb-1l)
IF start_ceilo GT 1 THEN time_c = time_c(0l:start_ceilo-1l)

rain_all_wvl = rain_all_wvl(0l:start_wvl-1l)

iwv      = iwv_all(0l:start_wvl-1l)
lwp      = lwp_all(0l:start_wvl-1l)
met_all  = met_all(*,0l:start_met-1l)
rain_all = rain_all(0l:start_met-1l)
IF start_ceilo GT 1 THEN lcb_all = lcb_all(0l:start_ceilo-1l)
temp_all = temp_all(*, 0l:start_blb-1)
hum_all = hum_all(*, 0l:start_blb-1)

;IF start_wvl_oz GT 1 THEN BEGIN
 iwv_oz = iwv_all_oz(0l:start_wvl_oz-1l)
 lwp_oz = lwp_all_oz(0l:start_wvl_oz-1l)
 az_oz = az_all_oz(0l:start_wvl_oz-1l)
 el_oz = el_all_oz(0l:start_wvl_oz-1l) 
 scan = scan_all(0l:start_wvl_oz-1l)
;ENDIF

;***convert TB time to decimal hours
xx_N = n_elements(time_w)
xx_dt = time_w[1:xx_N-1] - time_w[0:xx_N-2]
xx_i = where( xx_dt LT 0 , xx_cnt )
if xx_cnt GT 0 then print,'WARNING: unsorted wvl-time !'

i_sort = SORT(time_w)
time_w = time_w(i_sort)
iwv = iwv(i_sort)
lwp = lwp(i_sort)
rain_all_wvl = rain_all_wvl(i_sort)

CALDAT, time_w, MONTH, DAY, YEAR, hour_wvl, min_wvl, sec_wvl
CALDAT_TO_STRING, day, month, year, yymmdd

i_day = WHERE(yymmdd EQ date)
IF i_day(0) NE -1 THEN BEGIN
 time_w = hour_wvl(i_day) + min_wvl(i_day)/60d + sec_wvl(i_day)/3600d
 rain_all_wvl = rain_all_wvl(i_day)
 iwv = iwv(i_day)
 lwp = lwp(i_day)
ENDIF

ii = WHERE(iwv GT 0.)
iwv_min = -99.
iwv_max = 0.
IF ii(0) NE -1 THEN BEGIN 
 iwv_min = MEAN(iwv(ii))-3.*STDDEV(iwv(ii))
 IF iwv_min LT 0. THEN iwv_min = 0.
 iwv_max = MEAN(iwv(ii))+3.*STDDEV(iwv(ii))
 IF iwv_max GT 50. THEN iwv_max = 50. 
ENDIF

ii = WHERE(lwp GT -50.)
IF ii(0) NE -1 THEN BEGIN
 lwp_min = MEAN(lwp(ii))-3.*STDDEV(lwp(ii))
 IF lwp_min LT -50. THEN lwp_min = -49.
 lwp_max = MEAN(lwp(ii))+3.*STDDEV(lwp(ii))
 IF lwp_max GT 2000. THEN lwp_max = 1999. 
 IF lwp_max LT 300. THEN lwp_max = 300.
ENDIF ELSE BEGIN
 lwp_min = -49.
 lwp_max = 1000.
ENDELSE

;***convert TB off-zenith time to decimal hours
xx_N = n_elements(time_woz)
xx_dt = time_woz[1:xx_N-1] - time_woz[0:xx_N-2]
xx_i = where( xx_dt LT 0 , xx_cnt )
if xx_cnt GT 0 then print,'WARNING: unsorted wvl(woz)-time !'

i_sort = SORT(time_woz)
time_woz = time_woz(i_sort)
iwv_oz = iwv_oz(i_sort)
lwp_oz = lwp_oz(i_sort)
el_oz = el_oz(i_sort)
az_oz = az_oz(i_sort)
scan = scan(i_sort)

;**convert HATPRO azimuth to geographical azimuth
;az_cor = 180. ; radiometer pointing southward
az_cor = 0. ; radiometer pointing northward

IF az_cor EQ 180. THEN BEGIN
 i_180_0 = WHERE(az_oz GE 0. AND az_oz LE 180.)
 i_0_180 = WHERE(az_oz GT 180. AND az_oz LE 360.)
 az_oz_old = az_oz
 IF i_180_0(0) NE -1 THEN az_oz(i_180_0) = (az_oz(i_180_0)+az_cor) - 2*az_oz(i_180_0)
 IF i_0_180(0) NE -1 THEN az_oz(i_0_180) = (az_oz(i_0_180)+az_cor) - 2*az_oz(i_0_180) + 360.
ENDIF ELSE IF az_cor EQ 0. THEN BEGIN
 az_oz = 360. - az_oz
 i_360 = WHERE(az_oz EQ 360.)
 IF i_360(0) NE -1 THEN az_oz(i_360) = 0.
ENDIF

rain_all_wvl_oz = rain_all_wvl_oz(i_sort)

CALDAT, time_woz, MONTH, DAY, YEAR, hour_wvl, min_wvl, sec_wvl
CALDAT_TO_STRING, day, month, year, yymmdd

i_day = WHERE(yymmdd EQ date)  

IF i_day(0) NE -1 THEN BEGIN
 time_woz = hour_wvl(i_day) + min_wvl(i_day)/60d + sec_wvl(i_day)/3600d
 rain_all_wvl_oz = rain_all_wvl(i_day)
 iwv_oz = iwv_oz(i_day)
 lwp_oz = lwp_oz(i_day)
 el_oz = el_oz(i_day)
 az_oz = az_oz(i_day)
 scan = scan(i_day)
ENDIF ELSE BEGIN
 print, 'Warning: no off-zenith data'
 GOTO, SKIP_OZ
ENDELSE

;***convert met time to decimal hours
xx_N = n_elements(time_m)
xx_dt = time_m[1:xx_N-1] - time_m[0:xx_N-2]
xx_i = where( xx_dt LT 0 , xx_cnt )
if xx_cnt GT 0 then print,'WARNING: unsorted met-time !'

i_sort = SORT(time_m)
time_m = time_m(i_sort)
met_all = met_all(*, i_sort)
rain_all = rain_all(i_sort)

CALDAT, time_m, MONTH, DAY, YEAR, hour_met, min_met, sec_met
CALDAT_TO_STRING, day, month, year, yymmdd

i_day = WHERE(yymmdd EQ date)  

IF i_day(0) NE -1 THEN BEGIN
 time_m = hour_met(i_day) + min_met(i_day)/60d + sec_met(i_day)/3600d
 met_all = met_all(*, i_day)
 rain_all = rain_all(i_day)
ENDIF

;***convert blb time to decimal hours
xx_N = N_ELEMENTS(time_b)
xx_dt = time_b[1:xx_N-1] - time_b[0:xx_N-2]
xx_i = where( xx_dt LT 0 , xx_cnt )
if xx_cnt GT 0 then print,'WARNING: unsorted blb-time !'

i_sort = SORT(time_b)
time_b = time_b(i_sort)
temp_all = temp_all(*, i_sort)
hum_all = hum_all(*, i_sort)
rain_all_blb = rain_all_blb(i_sort)

CALDAT, time_b, MONTH, DAY, YEAR, hour_blb, min_blb, sec_blb
CALDAT_TO_STRING, day, month, year, yymmdd

i_day = WHERE(yymmdd EQ date)  

IF i_day(0) NE -1 THEN BEGIN
 time_b = hour_blb(i_day) + min_blb(i_day)/60d + sec_blb(i_day)/3600d
 temp_all = temp_all(*, i_day)
 hum_all = hum_all(*, i_day)
 rain_all_blb = rain_all_blb(i_day)
ENDIF ELSE BEGIN
 print, 'Warning: no BLB files!'
 GOTO, SKIP_TS
ENDELSE

;***bring profiles onto an equidistant time grid from 0 to 24 UTC

time_res = 900. ; (15 minutes)
hs = time_res/3600.
n_equi = 24.*3600./time_res
n_equi = n_equi+1
time_equi = FINDGEN(n_equi)*time_res
time_equi = time_equi/3600.
n_z = N_ELEMENTS(z)
dummy = -999.

T_equi = REPLICATE(dummy, n_equi, n_z)
q_equi = REPLICATE(dummy, n_equi, n_z)
rh_equi = REPLICATE(dummy, n_equi, n_z)
flag_equi = REPLICATE(0. ,n_equi)

FOR i = 0, n_equi-1 DO BEGIN

 i_match = WHERE(time_b GE time_equi(i)-hs AND time_b LT time_equi(i)+hs)
 i_match_rain = WHERE(time_m GE time_equi(i)-hs AND time_m LT time_equi(i)+hs)

 IF i_match_rain(0) NE -1 THEN BEGIN

  i_rain = WHERE(rain_all(i_match_rain) GT 0.)
  IF i_rain(0) NE -1 THEN flag_equi(i) = 1

 ENDIF

 thres_exceed = 0
 IF i_match(0) NE -1 THEN BEGIN
  
  FOR j = 0, n_z-1 DO BEGIN
   T_equi(i, j) = MEAN(temp_all(j, i_match))
   q_equi(i, j) = MEAN(hum_all(j, i_match))
   IF T_equi(i, j) LE t_thres_min OR T_equi(i, j) GT t_thres_max THEN thres_exceed = 1
   IF q_equi(i, j) LE q_thres_min OR q_equi(i, j) GT q_thres_max THEN thres_exceed = 1
  ENDFOR

;***calculate relative humidity profiles from T_equi and q_equi

  FOR j = 0, n_z-1 DO BEGIN
   ABSHUM_TO_RH, T_equi(i, j), q_equi(i, j)/1000., rhx
   IF rhx LE rh_thres_min OR rhx GT rh_thres_max THEN BEGIN
    IF z(j) LE 5000. THEN thres_exceed = 1
   ENDIF
   rh_equi(i, j) = rhx
  ENDFOR

;*bei level2a wir flag entweder auf 0 oder 1 gesetzt, um die Sachen einfach zu halten
  IF thres_exceed EQ 1 THEN flag_equi(i) = 1
 ENDIF
ENDFOR

;***read iwv(radiosonde Essen) (only jue)
iwv_radio=['-99','-99']

if station eq 'jue' then begin
 radiopath = '/data/data_hatpro/jue/data/radio_essen/'
 radiofile = radiopath + date_m + '.txt'
 str = ''
 s = ''
 tt = 0
 search=file_search(radiofile)
 if search(0) ne '' then begin
  OPENR, unit, radiofile, /GET_LUN
  WHILE tt lt 2 && (~ EOF(unit)) DO BEGIN
    READF, unit, s
    if strmid(s,0,2) EQ date_d THEN BEGIN 
     iwv_radio(tt) = STRMID(s, 6, 5)
     tt = tt+1
    endif
  ENDWHILE
  time_radio=[0.,12.]

  CLOSE, unit
  FREE_LUN, unit
 endif
endif

;****get information about campaign/station
GET_STATION_INFO, 0, station, start_time, stop_time, stat_name,$
                  latitude, longitude, altitude, instrument,$
                  comment_ret, comment_tech


;****plot LWP & rh-prof time series
!p.color = 0
!P.charsize = 0.8

psfile = '/data/data_hatpro/'+station+'/plots/level2a/'+date_m+'/'+date+'_h_'+station+'_l2a_rh.ps'
set_plot, 'ps'
!p.multi = [0, 1, 1]
LOADCT,39
device, /color, file = psfile, BITS_PER_PIXEL = 8, /landscape
!P.font = -1

h_min = 0.
h_max = 5.

rh_ret_new = rh_equi
time = time_equi
n_Tz = WHERE(z/1000. LE h_max)
nz = N_ELEMENTS(n_Tz)
;max_rh = MAX(rh_ret_new(*, 0))
;max_rh = max_rh + max_rh*0.002
max_rh = 100.
min_rh = 0.
nnn = N_ELEMENTS(n_Tz)
;iii = WHERE(rh_ret_new(*, n_Tz(nnn-1)) GT 150.)
;IF iii(0) NE -1 THEN min_rh = MIN(rh_ret_new(iii, n_Tz(nnn-1))) - 5.

arr_rh = DBZ_2_BYTE_T(rh_ret_new(*, n_Tz), max_rh, min_rh)
n_new = LONG(N_ELEMENTS(time))

time_start = time(0)
time_end = time(n_new-2l)

FOR i = 0l, n_new-1l DO BEGIN
 i_rh = WHERE(rh_ret_new(i, n_Tz) LT min_rh)
 IF i_rh(0) NE -1 THEN arr_rh(i, i_rh) = 255
ENDFOR

arr_z = (z(1)-z(0))*(FINDGEN((h_max - h_min)*1e3/(z(1)-z(0))) + 1)
arr_rh_new = FLTARR(n_new, N_ELEMENTS(arr_z))
FOR i_arr = 0, n_new-1 DO arr_rh_new(i_arr, *) = INTERPOL(arr_rh(i_arr, n_Tz), z(n_Tz), arr_z)
arr_rh = arr_rh_new

pos1 = [0.1, 0.1, 0.85, 0.65]

tit = 'Liquid clouds & rel. humidity, ' + station +', ' + date

PLOT, time, YRANGE = [h_min, h_max], xrange = [time_start, time_end],$
      color = 0, xtitle = 'Time [UTC]', xticklen = -0.02, yticklen = -0.015,$
      xstyle = 1, ystyle = 1, ytitle = 'Height [km]',$
      Position = pos1, /Nodata, /Normal, /Noerase, font = 3

dx = 5.

dxx = (max_rh - min_rh)/dx
ll = FINDGEN(FIX(dxx))*dx + min_rh
lll = N_ELEMENTS(ll)
clab = REPLICATE(1, lll)
ccol = FLTARR(lll)
FOR ic = 0, lll-1 DO BEGIN
 ccol(ic) = DBZ_2_BYTE_T(FLOAT(ll(ic)), max_rh, min_rh)
ENDFOR

i_valid = WHERE(rh_ret_new(*, 0) GT 0.)
n_valid = N_ELEMENTS(i_valid)
p_cont = FLTARR(n_valid, nz)
t_cont = FLTARR(n_valid)
k = 0

IF i_valid(0) NE -1 THEN BEGIN

 FOR iv = 0, n_valid-2 DO BEGIN

  p_cont(k, *) = rh_ret_new(i_valid(iv), 0:nz-1)
  t_cont(k, *) = time(i_valid(iv))

  IF i_valid(iv)+1 EQ i_valid(iv+1) THEN BEGIN

   k = k + 1

  ENDIF ELSE BEGIN
   IF k GT 0 THEN BEGIN


    p_cont = p_cont(0:k, *)
    t_cont = t_cont(0:k, *)

    CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
            Position = pos1, xrange = [time_start, time_end],$
            Levels = ll, xstyle = 5, ystyle = 5,$
            /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol, c_charsize = 0.8

    CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
            Position = pos1, xrange = [time_start, time_end],$
            Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5, c_charsize = 0.8,$
            /Noerase, /FOLLOW, color = 0, charsize = 0.8

    p_cont = FLTARR(n_valid, nz)
    t_cont = FLTARR(n_valid)
    k = 0
   ENDIF
  ENDELSE
 ENDFOR
 IF k GT 1 THEN BEGIN
  p_cont = p_cont(0:k-1, *)
  t_cont = t_cont(0:k-1, *)

  CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
           Position = pos1, xrange = [time_start, time_end],$
           Levels = ll, xstyle = 5, ystyle = 5, c_charsize = 0.8,$
           /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol

  CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
           Position = pos1, xrange = [time_start, time_end], c_charsize = 0.8,$
           Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5,$
           /Noerase, /FOLLOW, color = 0, charsize = 0.8

 ENDIF

ENDIF

bar_tit = 'RH [%]'
ncolor = 254
BarColor = BINDGEN(ncolor)
BarLabel = [min_rh, max_rh]
!P.Color = 0

COLLABEL2, BarColor, BarLabel, $
    CHARSIZE = 0.5, $
    DIRECTION = 1, FONT = 3, $
    LabFormat = '(I4)',LABOFFSET = 35, $
    POSITION = [0.93, 0.4, 0.98, 0.75]
XYOUTS, 0.93,0.35, bar_tit, color = 0, font = 3, charsize = 0.8, /Normal
XYOUTS, 0.02, 0.01, comment_ret, color = 254, charsize = 0.5, /normal

;***plot flags
xx = 0.6
x = [-(xx), -(xx), (xx), (xx)]
y = [-1, 1, 1,-1]
USERSYM, x, y, /fill

x_min = 0.1
x_max = 0.85
y_min = 0.88
y_max = 0.93

;****rain flag
index = WHERE(flag_equi AND 1,count)
index1 = WHERE(flag_equi EQ 0, count1)

PLOT, [0, 0], [24, 1], /NODATA, yrange=[0.9,1.1], ystyle = 1, title = tit,$
      xrange = [0., 24.], xstyle = 1, xcharsize = 0.0000001, ycharsize = 0.00001,$
      position = [x_min,y_min,x_max,y_max], yticks = 1, xticks = 1, /NOERASE, /NORMAL
IF count GT 0 THEN BEGIN
 IF count EQ 1 THEN BEGIN
  count = 2
  index = [index(0), index(0)]
 ENDIF
 h = REPLICATE(1, count)
 OPLOT, time(index), h, psym=8, thick=5, color=250
ENDIF
IF count1 GT 0 THEN BEGIN
 IF count1 EQ 1 THEN BEGIN
  count1 = 2
  index1 = [index1(0), index1(0)]
 ENDIF
 h = REPLICATE(1, count1)
 OPLOT, time(index1), h, psym=8, thick=5, color=150
ENDIF

XYOUTS,0.05,y_min+0.02,'DATA FLAG',alignment=0.5,charsize=0.6,/Normal

;***plot LWP time series

x_min = 0.1
x_max = 0.85
y_min = 0.65
y_max = 0.88

PLOT, [0, 0], [24, 1], /NODATA, yrange=[lwp_min, lwp_max], ystyle = 1,$
      xrange = [0., 24.], xstyle = 1, xcharsize = 0.0000001,$
      position = [x_min,y_min,x_max,y_max], xticks = 1, ytitle = 'LWP [gm!e-2!N]', /NOERASE, /NORMAL

;OPLOT, time_w, lwp, psym=4, symsize = 0.3, color = 50
OPLOT, time_w, lwp, thick = 3, color = 50
OPLOT, [time_w(0), time_w(N_ELEMENTS(time_w)-1)], [0., 0.], linestyle = 1
DEVICE, /CLOSE

;****plot IWV & q-prof time series
!p.color = 0
!P.charsize = 0.8

psfile = '/data/data_hatpro/'+station+'/plots/level2a/'+date_m+'/'+date+'_h_'+station+'_l2a_q.ps'
set_plot, 'ps'
!p.multi = [0, 1, 1]
LOADCT,39
device, /color, file = psfile, BITS_PER_PIXEL = 8, /landscape
!P.font = -1

h_min = 0.
h_max = 10.

q_ret_new = q_equi
time = time_equi
n_Tz = WHERE(z/1000. LE h_max)
nz = N_ELEMENTS(n_Tz)
;max_q = MAX(q_ret_new(*, 0))
max_q = 15.
max_q = max_q + max_q*0.002
min_q = 0.
nnn = N_ELEMENTS(n_Tz)
iii = WHERE(q_ret_new(*, n_Tz(nnn-1)) GT 150.)
IF iii(0) NE -1 THEN min_q = MIN(q_ret_new(iii, n_Tz(nnn-1))) - 5.

arr_q = DBZ_2_BYTE_Q(q_ret_new(*, n_Tz), max_q, min_q, 255, 30)
n_new = LONG(N_ELEMENTS(time))

time_start = time(0)
time_end = time(n_new-2l)

FOR i = 0l, n_new-1l DO BEGIN
 i_q = WHERE(q_ret_new(i, n_Tz) LT min_q)
 IF i_q(0) NE -1 THEN arr_q(i, i_q) = 255
ENDFOR

arr_z = (z(1)-z(0))*(FINDGEN((h_max - h_min)*1e3/(z(1)-z(0))) + 1)
arr_q_new = FLTARR(n_new, N_ELEMENTS(arr_z))
FOR i_arr = 0, n_new-1 DO arr_q_new(i_arr, *) = INTERPOL(arr_q(i_arr, n_Tz), z(n_Tz), arr_z)
arr_q = arr_q_new

pos1 = [0.1, 0.1, 0.85, 0.65]

tit = 'Water vapor column & abs. humidity, ' + station +', ' + date

PLOT, time, YRANGE = [h_min, h_max], xrange = [time_start, time_end],$
      color = 0, xtitle = 'Time [UTC]', xticklen = -0.02, yticklen = -0.015,$
      xstyle = 1, ystyle = 1, ytitle = 'Height [km]',$
      Position = pos1, /Nodata, /Normal, /Noerase, font = 3

dx = 1.
IF (max_q - min_q) LE 0. THEN BEGIN
 max_q = 15.
 min_q = 0.
ENDIF

dxx = (max_q - min_q)/dx
ll = FINDGEN(FIX(dxx))*dx + min_q
lll = N_ELEMENTS(ll)
clab = REPLICATE(1, lll)
ccol = FLTARR(lll)
FOR ic = 0, lll-1 DO BEGIN
 ccol(ic) = DBZ_2_BYTE_Q(FLOAT(ll(ic)), max_q, min_q, 255, 30)
ENDFOR

i_valid = WHERE(q_ret_new(*, 0) GT 0.)
n_valid = N_ELEMENTS(i_valid)
p_cont = FLTARR(n_valid, nz)
t_cont = FLTARR(n_valid)
k = 0

IF i_valid(0) NE -1 THEN BEGIN

 FOR iv = 0, n_valid-2 DO BEGIN

  p_cont(k, *) = q_ret_new(i_valid(iv), 0:nz-1)
  t_cont(k, *) = time(i_valid(iv))

  IF i_valid(iv)+1 EQ i_valid(iv+1) THEN BEGIN

   k = k + 1

  ENDIF ELSE BEGIN
   IF k GT 0 THEN BEGIN


    p_cont = p_cont(0:k, *)
    t_cont = t_cont(0:k, *)

    CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
            Position = pos1, xrange = [time_start, time_end],$
            Levels = ll, xstyle = 5, ystyle = 5,$
            /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol, c_charsize = 0.8

    CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
            Position = pos1, xrange = [time_start, time_end],$
            Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5, c_charsize = 0.8,$
            /Noerase, /FOLLOW, color = 0, charsize = 0.8

    p_cont = FLTARR(n_valid, nz)
    t_cont = FLTARR(n_valid)
    k = 0
   ENDIF
  ENDELSE
 ENDFOR
 IF k GT 1 THEN BEGIN
  p_cont = p_cont(0:k-1, *)
  t_cont = t_cont(0:k-1, *)

  CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
           Position = pos1, xrange = [time_start, time_end],$
           Levels = ll, xstyle = 5, ystyle = 5, c_charsize = 0.8,$
           /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol

  CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
           Position = pos1, xrange = [time_start, time_end], c_charsize = 0.8,$
           Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5,$
           /Noerase, /FOLLOW, color = 0, charsize = 0.8

 ENDIF

ENDIF

bar_tit = 'q [gm!e-3!N]'
ncolor = 254
BarColor1 = BINDGEN(ncolor)
BarColor=DBZ_2_BYTE_Q(FLOAT(BarColor1), 0., 255., 255., 30.)
BarColor=REVERSE(BarColor)
BarLabel = [min_q, max_q]
!P.Color = 0

COLLABEL2, BarColor, BarLabel, $
    CHARSIZE = 0.5, $
    DIRECTION = 1, FONT = 3, $
    LabFormat = '(I4)',LABOFFSET = 35, $
    POSITION = [0.93, 0.4, 0.98, 0.75]
XYOUTS, 0.93,0.35, bar_tit, color = 0, font = 3, charsize = 0.8, /Normal
XYOUTS, 0.02, 0.01, comment_ret, color = 254, charsize = 0.5, /normal

;***plot flags
xx = 0.6
x = [-(xx), -(xx), (xx), (xx)]
y = [-1, 1, 1,-1]
USERSYM, x, y, /fill

x_min = 0.1
x_max = 0.85
y_min = 0.88
y_max = 0.93

;****rain flag
index = WHERE(flag_equi AND 1,count)
index1 = WHERE(flag_equi EQ 0, count1)

PLOT, [0, 0], [24, 1], /NODATA, yrange=[0.9,1.1], ystyle = 1, title = tit,$
      xrange = [0., 24.], xstyle = 1, xcharsize = 0.0000001, ycharsize = 0.00001,$
      position = [x_min,y_min,x_max,y_max], yticks = 1, xticks = 1, /NOERASE, /NORMAL
IF count GT 0 THEN BEGIN
 IF count EQ 1 THEN BEGIN
  count = 2
  index = [index(0), index(0)]
 ENDIF
 h = REPLICATE(1, count)
 OPLOT, time(index), h, psym=8, thick=5, color=250
ENDIF
IF count1 GT 0 THEN BEGIN
 IF count1 EQ 1 THEN BEGIN
  count1 = 2
  index1 = [index1(0), index1(0)]
 ENDIF
 h = REPLICATE(1, count1)
 OPLOT, time(index1), h, psym=8, thick=5, color=150
ENDIF

XYOUTS,0.05,y_min+0.02,'DATA FLAG',alignment=0.5,charsize=0.6,/Normal

;***plot IWV time series

x_min = 0.1
x_max = 0.85
y_min = 0.65
y_max = 0.88


PLOT, [0, 0], [24, 1], /NODATA, yrange=[iwv_min, iwv_max], ystyle = 1,$
      xrange = [0., 24.], xstyle = 1, xcharsize = 0.0000001,$
      position = [x_min,y_min,x_max,y_max], xticks = 1, ytitle = 'IWV [kgm!e-2!N]', /NOERASE, /NORMAL

;OPLOT, time_w, iwv, psym=4, symsize = 0.3, color = 50
OPLOT, time_w, iwv, thick = 3, color = 50

DEVICE, /CLOSE

;****plot Temp-prof time series
!p.color = 0
!P.charsize = 0.8

psfile = '/data/data_hatpro/'+station+'/plots/level2a/'+date_m+'/'+date+'_h_'+station+'_l2a_t.ps'
set_plot, 'ps'
!p.multi = [0, 1, 1]
LOADCT,39
device, /color, file = psfile, BITS_PER_PIXEL = 8, /landscape
!P.font = -1

h_min = 0.
h_max = 5.

T_ret_new = T_equi
time = time_equi
n_Tz = WHERE(z/1000. LE h_max)
nz = N_ELEMENTS(n_Tz)
max_T = MAX(T_ret_new(*, 0))
max_T = max_T + max_T*0.002
min_T = 0.
nnn = N_ELEMENTS(n_Tz)
iii = WHERE(T_ret_new(*, n_Tz(nnn-1)) GT 150.)
IF iii(0) NE -1 THEN min_T = MIN(T_ret_new(iii, n_Tz(nnn-1))) - 5.

arr_T = DBZ_2_BYTE_T(T_ret_new(*, n_Tz), max_T, min_T)
n_new = LONG(N_ELEMENTS(time))

time_start = time(0)
time_end = time(n_new-2l)

FOR i = 0l, n_new-1l DO BEGIN
 i_T = WHERE(T_ret_new(i, n_Tz) LT min_T)
 IF i_T(0) NE -1 THEN arr_T(i, i_T) = 255
ENDFOR

arr_z = (z(1)-z(0))*(FINDGEN((h_max - h_min)*1e3/(z(1)-z(0))) + 1)
arr_T_new = FLTARR(n_new, N_ELEMENTS(arr_z))
FOR i_arr = 0, n_new-1 DO arr_T_new(i_arr, *) = INTERPOL(arr_T(i_arr, n_Tz), z(n_Tz), arr_z)
arr_T = arr_T_new

pos1 = [0.1, 0.1, 0.85, 0.80]

tit = 'Temperature, ' + station +', ' + date

PLOT, time, YRANGE = [h_min, h_max], xrange = [time_start, time_end],$
      color = 0, xtitle = 'Time [UTC]', xticklen = -0.02, yticklen = -0.015,$
      xstyle = 1, ystyle = 1, ytitle = 'Height [km]',$
      Position = pos1, /Nodata, /Normal, /Noerase, font = 3

dx = 2.
IF (max_T - min_T) LE 0. THEN BEGIN
 max_T = 300.
 min_T = 200.
ENDIF

dxx = (max_T - min_T)/dx
ll = FINDGEN(FIX(dxx))*dx + min_T
lll = N_ELEMENTS(ll)
clab = REPLICATE(1, lll)
ccol = FLTARR(lll)
FOR ic = 0, lll-1 DO BEGIN
 ccol(ic) = DBZ_2_BYTE_T(FLOAT(ll(ic)), max_T, min_T)
ENDFOR

i_valid = WHERE(T_ret_new(*, 0) GT 0.)
n_valid = N_ELEMENTS(i_valid)
p_cont = FLTARR(n_valid, nz)
t_cont = FLTARR(n_valid)
k = 0

IF i_valid(0) NE -1 THEN BEGIN

 FOR iv = 0, n_valid-2 DO BEGIN

  p_cont(k, *) = T_ret_new(i_valid(iv), 0:nz-1)
  t_cont(k, *) = time(i_valid(iv))

  IF i_valid(iv)+1 EQ i_valid(iv+1) THEN BEGIN

   k = k + 1

  ENDIF ELSE BEGIN
   IF k GT 0 THEN BEGIN


    p_cont = p_cont(0:k, *)
    t_cont = t_cont(0:k, *)

    CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
            Position = pos1, xrange = [time_start, time_end],$
            Levels = ll, xstyle = 5, ystyle = 5,$
            /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol, c_charsize = 0.8

    CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
            Position = pos1, xrange = [time_start, time_end],$
            Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5, c_charsize = 0.8,$
            /Noerase, /FOLLOW, color = 0, charsize = 0.8

    p_cont = FLTARR(n_valid, nz)
    t_cont = FLTARR(n_valid)
    k = 0
   ENDIF
  ENDELSE
 ENDFOR
 IF k GT 1 THEN BEGIN
  p_cont = p_cont(0:k-1, *)
  t_cont = t_cont(0:k-1, *)

  CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
           Position = pos1, xrange = [time_start, time_end],$
           Levels = ll, xstyle = 5, ystyle = 5, c_charsize = 0.8,$
           /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol

  CONTOUR, p_cont, t_cont, z(n_Tz)/1000.,$
           Position = pos1, xrange = [time_start, time_end], c_charsize = 0.8,$
           Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5,$
           /Noerase, /FOLLOW, color = 0, charsize = 0.8

 ENDIF

ENDIF

bar_tit = 'T [K]'
ncolor = 254
BarColor = BINDGEN(ncolor)
BarLabel = [min_T, max_T]
!P.Color = 0

COLLABEL2, BarColor, BarLabel, $
    CHARSIZE = 0.5, $
    DIRECTION = 1, FONT = 3, $
    LabFormat = '(I4)',LABOFFSET = 35, $
    POSITION = [0.93, 0.4, 0.98, 0.75]
XYOUTS, 0.93,0.35, bar_tit, color = 0, font = 3, charsize = 0.8, /Normal
XYOUTS, 0.02, 0.01, comment_ret, color = 254, charsize = 0.5, /normal

;***plot flags
xx = 0.6
x = [-(xx), -(xx), (xx), (xx)]
y = [-1, 1, 1,-1]
USERSYM, x, y, /fill

x_min = 0.1
x_max = 0.85
y_min = 0.8
y_max = 0.9

;****rain flag
index = WHERE(flag_equi AND 1,count)
index1 = WHERE(flag_equi EQ 0, count1)

PLOT, [0, 0], [24, 1], /NODATA, yrange=[0.9,1.1], ystyle = 1, title = tit,$
      xrange = [0., 24.], xstyle = 1, xcharsize = 0.0000001, ycharsize = 0.00001,$
      position = [x_min,y_min,x_max,y_max], yticks = 1, xticks = 1, /NOERASE, /NORMAL
IF count GT 0 THEN BEGIN
 IF count EQ 1 THEN BEGIN
  count = 2
  index = [index(0), index(0)]
 ENDIF
 h = REPLICATE(1, count)
 OPLOT, time(index), h, psym=8, thick=5, color=250
ENDIF
IF count1 GT 0 THEN BEGIN
 IF count1 EQ 1 THEN BEGIN
  count1 = 2
  index1 = [index1(0), index1(0)]
 ENDIF
 h = REPLICATE(1, count1)
 OPLOT, time(index1), h, psym=8, thick=5, color=150
ENDIF

XYOUTS,0.05,y_min+0.075,'DATA FLAG',alignment=0.5,charsize=0.6,/Normal

DEVICE, /CLOSE

SKIP_TS:

;****plot IWV azimuth-time series
!P.charsize = 0.8

IF start_wvl GT 2 THEN BEGIN 

 psfile = '/data/data_hatpro/'+station+'/plots/level2a/'+date_m+'/'+date+'_h_'+station+'_l2a_iwv_az.ps'
 set_plot, 'ps'
 !p.multi = [0, 2, 1]
 LOADCT, 5
 device, /color, file = psfile, BITS_PER_PIXEL = 8, /landscape
 !P.font = 3

 el_plot = [45., 30.]
 n_el = N_ELEMENTS(el_plot)

 plot_off_xx = 0.5
;r_min = FIX(MIN(iwv_oz))-1.
;r_max = FIX(MAX(iwv_oz))+1

 r_min = FLOOR(iwv_min)
 r_max = CEIL(iwv_max)
 IF r_min LT 0. THEN r_min = 0.
 IF r_max GT 50. THEN r_max = 50.

 FOR i_el = 0, n_el-1 DO BEGIN
  plot_off_x = FLOAT(i_el)*plot_off_xx
  el_c = STRING(el_plot(i_el), format = '(i2)')

  tit = 'IWV (az_vs_time)@'+el_c+' deg, ' + station +', ' + date

  PLOT, [0,360], [0,24], /NODATA, Title=tit,$
  xrange = [0, 360], xstyle=1, xcharsize = 1.2, ycharsize = 1.2,$
  xticklen=0.02, xticks=4, xtickname = ['N','E','S','W','N'],$
  yrange = [0,24], ystyle=1, ytitle = 'Time of day [UTC]',$
  yticklen = -0.02, yticks=8, position=[0.1+plot_off_x,0.15,0.4+plot_off_x,0.9]

  i_plot = WHERE(el_oz EQ el_plot(i_el))
  IF i_plot(0) NE -1 THEN BEGIN
   n_plot = N_ELEMENTS(i_plot)
   time_plot = time_woz(i_plot)
   az_plot = az_oz(i_plot)
   iwv_plot = iwv_oz(i_plot)
   scan_plot = scan(i_plot)
   max_scan = MAX(scan_plot(i_plot))
   min_scan = MIN(scan_plot(i_plot))

   da = (az_plot(1)-az_plot(2))
   dt = 0.28

   FOR i_s = min_scan, max_scan-1 DO BEGIN
    i_scan = WHERE(scan_plot EQ i_s)  
    IF i_s LT max_scan THEN BEGIN 
     i_scan_next = WHERE(scan_plot EQ i_s+1)
     dt = time_plot(i_scan_next(0))-time_plot(i_scan(0))
    ENDIF 

    n_scan = N_ELEMENTS(i_scan)
 
    IF i_scan(0) NE -1 THEN BEGIN
     tim = time_plot(i_scan(0))
     FOR i_a = 0, n_scan-1 DO BEGIN
      x = [az_plot(i_scan(i_a)), az_plot(i_scan(i_a)), az_plot(i_scan(i_a))+da, az_plot(i_scan(i_a))+da]
      y = [tim, tim+dt, tim+dt, tim]
      col = iwv_plot(i_scan(i_a))
      IF ABS(col-MEAN(iwv_plot(i_scan))) GT 2. AND STDDEV(iwv_plot(i_scan)) GT 1.0 THEN col = r_max
      IF col GT r_max THEN col = r_max
      IF col LT r_min THEN col = r_min
      POLYFILL, x, y, color=(r_max-col)*255./(r_max-r_min)   

     ENDFOR
    ENDIF  
   ENDFOR

   PLOT, [0,360], [0,24], /NODATA,$
   xrange = [0, 360], xstyle=1, xcharsize = 1e-5, ycharsize = 1e-5,$
   xticklen=0.02, xticks=4, xtickname = ['N','E','S','W','N'],$
   yrange = [0,24], ystyle=1, ytitle = 'Time of day [UTC]',$
  yticklen = -0.02, yticks=8, position=[0.1+plot_off_x,0.15,0.4+plot_off_x,0.9], /NOERASE

  ENDIF
 
 ENDFOR

 bar = 255-INDGEN(255)
 n_ticks =5
 print,'Min',r_min,r_max

 FOR i=0,n_ticks DO BEGIN
   z = r_min + (r_max-r_min)*i/n_ticks
   x = 0.2 + i*0.6/n_ticks
   XYOUTS,x,0.08,'|',alignment=0.5,size=0.8,/NORMAL
   XYOUTS,x,0.05,STRING(z,FORMAT='(F4.1)'),charsize=0.8,alignment=0.5,/NORMAL
 ENDFOR
 XYOUTS,0.5,0.015,'IWV [kg m!e-2!n] (airmass-corrected)',alignment=0.5, charsize=1.2,/NORMAL

 TV,bar,0.2,0.07,xsize=0.6,ysize=0.02,/NORMAL

 DEVICE, /CLOSE

;****plot LWP azimuth-time series
 !P.charsize = 0.8

 psfile = '/data/data_hatpro/'+station+'/plots/level2a/'+date_m+'/'+date+'_h_'+station+'_l2a_lwp_az.ps'
 set_plot, 'ps'
 !p.multi = [0, 2, 1]
 LOADCT, 1
 device, /color, file = psfile, BITS_PER_PIXEL = 8, /landscape
 !P.font = 3

 el_plot = [45., 30.]
 n_el = N_ELEMENTS(el_plot)

 plot_off_xx = 0.5
;r_min = 0.
;r_max = FIX(MAX(lwp_oz*1000.))

;IF r_min LT 0. THEN r_min = 0.
;IF r_max GT 2000. THEN r_max = 2000.

 r_min = 0.
 r_max = CEIL(lwp_max)
 
 FOR i_el = 0, n_el-1 DO BEGIN
  plot_off_x = FLOAT(i_el)*plot_off_xx
  el_c = STRING(el_plot(i_el), format = '(i2)')
 
  tit = 'LWP (az_vs_time)@'+el_c+' deg, ' + station +', ' + date

  PLOT, [0,360], [0,24], /NODATA, Title=tit,$
  xrange = [0, 360], xstyle=1, xcharsize = 1.2, ycharsize = 1.2,$
  xticklen=0.02, xticks=4, xtickname = ['N','E','S','W','N'],$
  yrange = [0,24], ystyle=1, ytitle = 'Time of day [UTC]',$
  yticklen = -0.02, yticks=8, position=[0.1+plot_off_x,0.15,0.4+plot_off_x,0.9]

  i_plot = WHERE(el_oz EQ el_plot(i_el))
  IF i_plot(0) NE -1 THEN BEGIN
   n_plot = N_ELEMENTS(i_plot)
   time_plot = time_woz(i_plot)
   az_plot = az_oz(i_plot)
   lwp_plot = lwp_oz(i_plot)*1000.
   scan_plot = scan(i_plot)
   max_scan = MAX(scan_plot(i_plot))
   min_scan = MIN(scan_plot(i_plot))

   da = (az_plot(1)-az_plot(2))
   dt = 0.28

   FOR i_s = min_scan, max_scan-1 DO BEGIN
    i_scan = WHERE(scan_plot EQ i_s)  
    IF i_s LT max_scan THEN BEGIN 
     i_scan_next = WHERE(scan_plot EQ i_s+1)
     dt = time_plot(i_scan_next(0))-time_plot(i_scan(0))
    ENDIF 

    n_scan = N_ELEMENTS(i_scan)

    IF i_scan(0) NE -1 THEN BEGIN
     tim = time_plot(i_scan(0))
     FOR i_a = 0, n_scan-1 DO BEGIN
      x = [az_plot(i_scan(i_a)), az_plot(i_scan(i_a)), az_plot(i_scan(i_a))+da, az_plot(i_scan(i_a))+da]
      y = [tim, tim+dt, tim+dt, tim]      
      col = lwp_plot(i_scan(i_a))
;     IF col LT lwp_min THEN col = lwp_min
      IF col LT 0. THEN col = 0.
      IF col GT lwp_max THEN col = lwp_max
;     IF ABS(col-MEAN(lwp_plot(i_scan))) GT 2. AND STDDEV(lwp_plot(i_scan)) GT 1.0 THEN col = r_max
      POLYFILL, x, y, color=(r_max-col)*255./(r_max-r_min)   

     ENDFOR
    ENDIF  
   ENDFOR

   PLOT, [0,360], [0,24], /NODATA,$
   xrange = [0, 360], xstyle=1, xcharsize = 1e-5, ycharsize = 1e-5,$
   xticklen=0.02, xticks=4, xtickname = ['N','E','S','W','N'],$
   yrange = [0,24], ystyle=1, ytitle = 'Time of day [UTC]',$
  yticklen = -0.02, yticks=8, position=[0.1+plot_off_x,0.15,0.4+plot_off_x,0.9], /NOERASE

  ENDIF

 ENDFOR

 bar = 255-INDGEN(255)
 n_ticks =5
 print,'Min',r_min,r_max

 FOR i=0,n_ticks DO BEGIN
   z = r_min + (r_max-r_min)*i/n_ticks
   x = 0.2 + i*0.6/n_ticks
   XYOUTS,x,0.08,'|',alignment=0.5,size=0.8,/NORMAL
   XYOUTS,x,0.05,STRING(z,FORMAT='(F6.1)'),charsize=0.8,alignment=0.5,/NORMAL
 ENDFOR
 XYOUTS,0.5,0.015,'LWP [g m!e-2!n] (airmass-corrected)',alignment=0.5, charsize=1.2,/NORMAL

 TV,bar,0.2,0.07,xsize=0.6,ysize=0.02,/NORMAL
 
 DEVICE, /CLOSE

ENDIF; start_wvl GT 2

SKIP_OZ:

;****plot radar&LWP time series
!P.charsize = 1.8

psfile = '/data/data_hatpro/'+station+'/plots/level2a/'+date_m+'/'+date+'_h_'+station+'_l2a_wora.ps'
set_plot, 'ps'
!p.multi = [0, 1, 4]
LOADCT, 1
device, /color, file = psfile, BITS_PER_PIXEL = 8, /landscape
!P.font = 3

path_cr = '/data/data_hatpro/jue/data/mmclx/incoming/'
res = FILE_SEARCH(path_cr+'20'+date_m+date_d+'*.mmclx.gz')

IF res(0) NE '' THEN BEGIN
print, 'unzipping radar file: ', res
SPAWN, 'gunzip -f '+res
 
res = FILE_SEARCH(path_cr+'20'+date_m+date_d+'*.mmclx')

IF res(0) NE '' THEN BEGIN
 print, 'reading cloud radar'
;*ze(range, time)
 READ_CLDRAD_NC, res(0), time_cr, range, ze, ldr, vd, sigma
 time_cr = time_cr/86400d + JULDAY(1,1,1970,0,0,0)
 CALDAT, time_cr, month_cr, day_cr, year_cr, hour_cr, minute_cr, second_cr
 time_cr = hour_cr + minute_cr/60. + second_cr/3600.
ENDIF

print, 'zipping radar file: ', res
SPAWN, 'gzip -f '+res

;***plot ze, vd, sigma time series
pos_lw = [0.1, 0.78, 0.9, 0.9]
pos_ze = [0.1, 0.55, 0.9, 0.78]
pos_vd = [0.1, 0.32, 0.9, 0.55]
pos_si = [0.1, 0.09, 0.9, 0.32]

;***convert to dbz
ze = 10d*ALOG10(ze)

;**parameter ranges
ze_min  = -60.
ze_max  =  15.
ldr_min = -35.
ldr_max =  -15.
vd_max =   2.
vd_min =  -10.
sigma_min =  0.
sigma_max =  2.

range_max =  9999.
range_min = MIN(range)
range_ind = WHERE(range GE range_min and range LE range_max)

ydim = N_ELEMENTS(time_cr)
xdim = N_ELEMENTS(range_ind)
arr_ze  = BYTARR(xdim, ydim)
arr_ldr = BYTARR(xdim, ydim)
arr_vd = BYTARR(xdim, ydim)
arr_sigma = BYTARR(xdim, ydim)

;** convert ze to 0...255
arr_ze = (ABS(ze_min) + ze(range_ind,*)) / ((ze_max-ze_min)/255)
arr_ldr = (ABS(ldr_min) + ldr(range_ind,*)) / ((ldr_max-ldr_min)/255)
arr_vd = (ABS(vd_min) + vd(range_ind,*)) / ((vd_max-vd_min)/255);*(-1)+255
arr_sigma = (ABS(sigma_min) + sigma(range_ind,*)) / ((sigma_max-sigma_min)/255)
;** fill NANs
FOR i = 0, ydim-1 DO BEGIN
 ii = WHERE(~FINITE(arr_ze(*, i)) OR arr_ze(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_ze(ii, i) = 255b
 ii = WHERE(~FINITE(arr_vd(*, i)) OR arr_vd(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_vd(ii, i) = 255b
 ii = WHERE(~FINITE(arr_sigma(*, i)) OR arr_sigma(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_sigma(ii, i) = 255b
ENDFOR

;**plot ze
LOADCT, 39
TV, TRANSPOSE(arr_ze), pos_ze(0), pos_ze(1), xsize=pos_ze(2)-pos_ze(0), ysize=pos_ze(3)-pos_ze(1), /NORMAL
LOADCT, 0
PLOT, time_cr, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [m AGL]',$
      xrange=[time_cr(0), time_cr(ydim-1)], xstyle=1, ystyle=1,$
      xcharsize=1e-5, /NOERASE, position = pos_ze, yticklen = -0.009, xticklen = 0.04

n_ticks = 5
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = ze_min + (ze_max-ze_min)*i/n_ticks
ENDFOR

LOADCT, 39
COLORBAR, position = [pos_ze[2]+0.02, pos_ze[1], pos_ze[2]+0.04, pos_ze[3]-(pos_ze[3]-pos_ze[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(I4)'),$
          title = 'Ze [dBZe]', charsize=1.5, color=0 , font=0, right=1, /vertical

;*** plot Doppler velocity
LOADCT, 39
TV, TRANSPOSE(arr_vd), pos_vd(0), pos_vd(1), xsize=pos_vd(2)-pos_vd(0), ysize=pos_vd(3)-pos_vd(1), /NORMAL
LOADCT, 0
PLOT, time_cr, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [m AGL]',$
      xrange=[time_cr(0), time_cr(ydim-1)], xstyle=1, ystyle=1,$
      xcharsize=1e-5, /NOERASE, position = pos_vd, yticklen = -0.009, xticklen = 0.04

n_ticks = 5
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = vd_min + (vd_max-vd_min)*i/n_ticks
ENDFOR

LOADCT, 39
COLORBAR, position = [pos_vd[2]+0.02, pos_vd[1], pos_vd[2]+0.04, pos_vd[3]-(pos_vd[3]-pos_vd[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(f5.1)'),$
          title = 'Dop. vel. [m/s]', charsize=1.5, color=0 , font=0, right=1, /vertical

;*** plot sigma
LOADCT, 39
TV, TRANSPOSE(arr_sigma), pos_si(0), pos_si(1), xsize=pos_si(2)-pos_si(0), ysize=pos_si(3)-pos_si(1), /NORMAL
LOADCT, 0
PLOT, time_cr, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [m AGL]',$
      xrange=[time_cr(0), time_cr(ydim-1)], xstyle=1, ystyle=1,$
      /NOERASE, position = pos_si, yticklen = -0.009, xticklen = 0.04,$
      xtitle = 'Time [UTC]

n_ticks = 5
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = sigma_min + (sigma_max-sigma_min)*i/n_ticks
ENDFOR

LOADCT, 39
COLORBAR, position = [pos_si[2]+0.02, pos_si[1], pos_si[2]+0.04, pos_si[3]-(pos_si[3]-pos_si[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(F4.1)'),$
          title = 'Spectral width [m/s]', charsize=1.5, color=0 , font=0, right=1, /vertical

;***plot LWP time series

IF start_wvl GT 2 THEN BEGIN

 ii = WHERE(lwp GT -50.)
 tit = 'LWP & Cloud radar, ' + station +', ' + date


 PLOT, [0, 0], [24, 1], /NODATA, yrange=[lwp_min, lwp_max], ystyle = 1,$
       xrange = [time_w(0), time_w(N_ELEMENTS(time_w)-1)], xstyle = 1, xcharsize = 0.0000001, color=0, title=tit,$
       position = pos_lw, xticks = 1, ytitle = 'LWP [gm!e-2!N]', /NOERASE, /NORMAL
     

 OPLOT, time_w, lwp, thick = 3, color = 50
 OPLOT, [time_w(0), time_w(N_ELEMENTS(time_w)-1)], [0., 0.], linestyle = 1

ENDIF

DEVICE, /CLOSE

ENDIF ELSE BEGIN

 print, 'no radar files found for ', date

ENDELSE

END
