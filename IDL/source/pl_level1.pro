;+
;**********************
PRO PL_LEVEL1A_HATPRO,$
;**********************
;keywords
Date=date,$              ; YYMMDD to be plotted
Station=station,$        ; three letter code
Algo=algo,$              ; name of algorithm files without 'IWV_' and 'LWV_' in folders '/home/hatpro/retrievals/IWV/' and '.../LWP/'
IWV_max=IWV_max, $       ; maximum IWV to be plotted, default=80., if IWV_max=0 autmatically adapted to stddev during the day, if IWV_max<0 range=[mean-|IWV_max/2|,mean+|IWV_max/2|]
ceilo_type=ceilo_type, $ ; either 'igm' or 'dwd' currently
rd_file=rd_file,$        ; 'wvl': use WVL/OLC files (default), 'brt': use BRT files
trk = trk,$              ; ='non': do nothing (default), ='sat': write line of sight satellite IWV retrievals into ascii file
verbose=verbose           
; $Id: pl_level1a_hatpro.pro,v 1.10 2010/02/09 09:54:41 loehnert Exp hatpro $
; Abstract: 
; * plots level 1a HATPRO data - IWV & LWP, only intended for first quicklooks 
; Authors:
; S. Crewell, U. Loehnert
; Date:
; 2008-08-21
; Dependencies:
; - 
; Changes
; 2009-11-09
;   added rd_file keyword and subsequent processing
; 2009-10-08:
;   added absolute humidity (JHS)
; 2009-10-01: 
;   added water vapor mixing ratio in humidity plot (JHS)
; 2008-10-27: included ceilo_type keyword (UL)
;-
!Except=2
if ~keyword_set(verbose) then verbose=0
if ~keyword_set(IWV_max) then IWV_max=80.

IF N_ELEMENTS(rd_file) EQ 0 THEN rd_file = 'wvl'
IF N_ELEMENTS(rd_file) NE 0 THEN BEGIN
 IF rd_file NE 'brt' THEN rd_file = 'wvl'
ENDIF

IF N_ELEMENTS(ceilo_type) EQ 0 THEN ceilo_type = 'xxx'
IF N_ELEMENTS(trk) EQ 0 THEN trk = 'non'

date_m= STRMID(date, 0, 4)
date_d = STRMID(date, 4, 2)
date_y = STRMID(date, 0, 2)
date_m_x = STRMID(date, 2, 2)

raw_path = '/data/data_hatpro/'+station+'/data/raw/'+date_m+'/'+date+'/'
plot_path = '/data/data_hatpro/'+station+'/plots/level1a/'+date_m+'/'
sat_path = '/data/data_hatpro/'+station+'/data/sat/'

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
 in_name_brt_today = FILE_SEARCH(raw_path + '*.brt')
 in_name_brt_prev = FILE_SEARCH(raw_path_prev + '*.brt')

 in_name_brt = [in_name_brt_prev, in_name_brt_today]
 IF in_name_brt_prev(0) EQ '' THEN in_name_brt = [in_name_brt_today]
 IF in_name_brt_today(0) EQ '' THEN in_name_brt = [in_name_brt_prev]

 n_in_wvl = N_ELEMENTS(in_name_brt)
 IF in_name_brt(0) EQ '' THEN n_in_wvl = 0 

 if verbose then PRINT,'Total number of BRT files being checked: ',n_in_wvl

ENDELSE

n_in_brt_sat = 0
n_in_trk_sat = 0
IF trk EQ 'sat' THEN BEGIN

 ; search and count sat BRT files
 in_name_brt_sat = FILE_SEARCH(raw_path + '*sat.brt')
 in_name_brt_sat_prev = FILE_SEARCH(raw_path_prev + '*sat.brt')  

 IF in_name_brt_sat_prev(0) NE '' THEN in_name_brt_sat = [in_name_brt_sat_prev, in_name_brt_sat]

 n_in_brt_sat = N_ELEMENTS(in_name_brt_sat)
 IF in_name_brt_sat(0) EQ '' THEN n_in_brt_sat = 0 

 if verbose then PRINT,'Total number of SAT BRT files being checked: ',n_in_brt_sat

 ; search and count sat files
 in_name_trk_sat = FILE_SEARCH(raw_path + '*sat.trk')
 in_name_trk_sat_prev = FILE_SEARCH(raw_path_prev + '*sat.trk')  

 IF in_name_trk_sat_prev(0) NE '' THEN in_name_trk_sat = [in_name_trk_sat_prev, in_name_trk_sat]

 n_in_trk_sat = N_ELEMENTS(in_name_trk_sat)
 IF in_name_trk_sat(0) EQ '' THEN n_in_trk_sat = 0 

 if verbose then PRINT,'Total number of SAT TRK files being checked: ',n_in_trk_sat

;***read .trk files of this day
 FOR i = 0, n_in_trk_sat-1 DO BEGIN

  READ_HATPRO_TRK, in_name_trk_sat(i), sat_type, sat_nr, yys, mms, dds, hhs, mis, sss, sat_elev, sat_azi, err
  IF err EQ 1 THEN GOTO, SKIP_THIS_TRK

  n_sat = N_ELEMENTS(sss)
  jds = REPLICATE(0d, n_sat)

  FOR j = 0, n_sat-1 DO BEGIN
   yyc = STRING(yys(j), format = '(i2)')
   IF yys(j) LT 10. THEN yyc = '0' + STRING(yys(j), format = '(i1)')
   mmc = STRING(mms(j), format = '(i2)')
   IF mms(j) LT 10. THEN mmc = '0' + STRING(mms(j), format = '(i1)')
   ddc = STRING(dds(j), format = '(i2)')
   IF dds(j) LT 10. THEN ddc = '0' + STRING(dds(j), format = '(i1)')
   dates = '20'+yyc+mmc+ddc

   times = hhs(j)+mis(j)/60.+sss(j)/3600.
   
   DATE_TIME_TO_JULDAT, dates, times, jd
   jds(j) = jd
  ENDFOR

  IF i EQ 0 THEN BEGIN
   jds_a = jds
   sat_type_a = sat_type
   sat_nr_a = sat_nr
   sat_elev_a = sat_elev
   sat_azi_a = sat_azi        
  ENDIF ELSE BEGIN
   jds_a = [jds_a, jds]
   sat_type_a = [sat_type_a, sat_type]
   sat_nr_a = [sat_nr_a, sat_nr]
   sat_elev_a = [sat_elev_a, sat_elev]
   sat_azi_a = [sat_azi_a, sat_azi]
  ENDELSE  

SKIP_THIS_TRK:

 ENDFOR

ENDIF

; ceilometer data
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

in_name_met = FILE_SEARCH(raw_path + '*.met')
in_name_met_prev = FILE_SEARCH(raw_path_prev + '*.met')

IF in_name_met_prev(0) NE '' THEN in_name_met = [in_name_met_prev, in_name_met]
n_in_met = N_ELEMENTS(in_name_met)
IF in_name_met(0) EQ '' THEN n_in_met = 0

if verbose then PRINT,'Total number of MET files being checked: ',n_in_met

n = 200000
n_z = 35
time_w = LONARR(n)
time_m = LONARR(n)
time_c = FLTARR(n)

iwv_all = FLTARR(n)
lwp_all = FLTARR(n)
rain_all_wvl = BYTARR(n)
rain_all = BYTARR(n)
met_all  = FLTARR(3, n)
temp_all = FLTARR(n_z, n)
hum_all =  FLTARR(n_z, n)
lcb_all = FLTARR(n)

start_wvl = 0l
start_met = 0l
start_ceilo  = 0l

coeff_dat =  '/home/hatpro/retrievals/IWV/IWV_'+algo
GET_COEFF_INT, 0, coeff_dat, lin_qua, angle, f, offsetv, c_lv, c_qv
coeff_dat =  '/home/hatpro/retrievals/LWP/LWP_'+algo
GET_COEFF_INT, 0, coeff_dat, lin_qua, angle, f, offsetl, c_ll, c_ql

FOR i_f = 0,n_in_wvl-1 DO BEGIN

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

 IF N_ELEMENTS(f) EQ 6 AND station EQ 'amm' THEN BEGIN ; 23.84 GHz channel disturbed
  nn = 6
  TB_new = FLTARR(nn, n_wvl)
  TB_new(0:1, *) = TB_wvl(0:1, *)   ; brightness temperatures without 23.84 GHz
  TB_new(2:5, *) = TB_wvl(3:6, *)   ; brightness temperatures without 23.84 GHz
 ENDIF

 IF N_ELEMENTS(f) EQ 3 AND station EQ 'jue' THEN BEGIN ; 22.24, 23.04, 25.44, 26.24 GHz channels disturbed
  nn = 3
  TB_new = FLTARR(nn, n_wvl)
  TB_new(0, *) = TB_wvl(2, *)      ; 23.84 GHz 
  TB_new(1:2, *) = TB_wvl(5:6, *)  ; 27.84, 31.4 GHz
 ENDIF

 FOR il = 0l, n_wvl_z-1l DO iwv(il) = TOTAL(c_lv*TB_new(0:nn-1, i_zen(il))) + TOTAL(c_qv*(TB_new(0:nn-1, i_zen(il))^2)) + iwv(il)
;**lwp

 IF N_ELEMENTS(f) EQ 6 AND station EQ 'amm' THEN BEGIN ;23.84 GHz channel disturbed
  nn = 6
  TB_new = FLTARR(nn, n_wvl)
  TB_new(0:1, *) = TB_wvl(0:1, *)   ; brightness temperatures without 23.84 GHz
  TB_new(2:5, *) = TB_wvl(3:6, *)   ; brightness temperatures without 23.84 GHz
 ENDIF

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


IF trk EQ 'sat' AND n_in_brt_sat GE 1 AND n_in_trk_sat GE 1 THEN BEGIN

 go = 0
 FOR i_s = 0, n_in_brt_sat-1 DO BEGIN

  if verbose GE 5 then print,'reading "', in_name_brt_sat(i_s),'"'

  GET_BRT, 0, in_name_brt_sat(i_s), n_wvls, time_wvls, rain_wvls, freq_brt, TB_brts, az_wvls, el_wvls

;**azimuth correction to geo-system
  IF n_wvls GT 1 THEN BEGIN
   az_wvls_cor = 0.
   i_180_0 = WHERE(az_wvls GE 0. AND az_wvls LE 180.)
   i_0_180 = WHERE(az_wvls GT 180. AND az_wvls LE 360.)
   IF i_180_0(0) NE -1 THEN az_wvls(i_180_0) = (az_wvls(i_180_0)+az_wvls_cor) - 2*az_wvls(i_180_0)
   IF i_0_180(0) NE -1 THEN az_wvls(i_0_180) = (az_wvls(i_0_180)+az_wvls_cor) - 2*az_wvls(i_0_180) + 360.
   ikk = WHERE(az_wvls LT 0.)
   IF ikk(0) NE -1 THEN az_wvls(ikk) = az_wvls(ikk)+360.
 
   IF n_wvls GT 0 THEN BEGIN
    TB_wvls = TB_brts(0:6, *)
    freqs = freq_brt
   ENDIF

   IF station EQ 'jue' THEN BEGIN ; 22.24, 23.04, 25.44, 26.24 GHz channels disturbed
    nn = 3
    TB_news = FLTARR(nn, n_wvls)
    TB_news(0, *) = TB_wvls(2, *)      ; 23.84 GHz 
    TB_news(1:2, *) = TB_wvls(5:6, *)   ; 27.84, 31.4 GHz
   ENDIF

   coeff_dat = file_search('/home/loehnert/retrieval/stat_algo/RET_files/IWV_GPS/IWV_db_0112_V*.RET')
   RET_IWVLWP_MA, coeff_dat, TB_news, el_wvls, iwvs, f

   IF go EQ 0 THEN BEGIN

    times_all = time_wvls
    rains_all = rain_wvls
    azs_all = az_wvls
    els_all = el_wvls
    iwvs_all = iwvs
    go = 1

   ENDIF ELSE BEGIN

    times_all = [times_all, time_wvls]
    rains_all = [rains_all, rain_wvls]
    azs_all = [azs_all, az_wvls]
    els_all = [els_all, el_wvls]
    iwvs_all = [iwvs_all, iwvs]
    
   ENDELSE 

  ENDIF

 ENDFOR

 times_all = times_all/86400d
 times_all = times_all + JULDAY(1,1,2001,0,0,0)

 i_sort = SORT(times_all)
 times_all = times_all(i_sort)
 iwvs_all = iwvs_all(i_sort)
 rains_all = rains_all(i_sort) 
 azs_all = azs_all(i_sort)
 els_all = els_all(i_sort)

 n_all = N_ELEMENTS(times_all)
 sat_type_all = REPLICATE('dummy', n_all)
 sat_nr_all = REPLICATE(-99., n_all)
 sat_elev_all = REPLICATE(-99., n_all)
 sat_azi_all = REPLICATE(-99., n_all) 

;***find closest .trk time
 FOR i = 0, n_all-1 DO BEGIN
;  ind = WHERE(ABS(times_all(i)-jds_a) LT 1d/86400d AND ABS(azs_all(i)-sat_azi_a) LT 0.01 AND ABS(els_all(i)-sat_elev_a) LT 0.5) 
   ind = WHERE(ABS(times_all(i)-jds_a) LT .5/86400d AND ABS(els_all(i)-sat_elev_a) LT 0.5)
;  CLOSEST, times_all(i), jds_a, dmin, ind
;  IF dmin LT 1d/86400d AND ABS(azs_all(i)-sat_azi_a(ind)) LT 0.01 THEN BEGIN

;IF ind(0) NE -1 then print, i, ind, ABS(times_all(i)-jds_a(ind)), 1d/86400d, ABS(azs_all(i)-sat_azi_a(ind)), 0.01, ABS(els_all(i)-sat_elev_a(ind)), 0.5
;aa=''
;read,aa
  IF N_ELEMENTS(ind) EQ 1 AND ind(0) NE -1 THEN BEGIN
   sat_type_all(i) = sat_type_a(ind)
   sat_nr_all(i) = sat_nr_a(ind) 
   sat_elev_all(i) = sat_elev_a(ind)
   sat_azi_all(i) = sat_azi_a(ind)  
  ENDIF
 ENDFOR

 CALDAT, times_all, MONTH, days, YEAR, hour_wvls, min_wvls, sec_wvls
 i_day = WHERE(days EQ date_d)

 IF i_day(0) NE -1 THEN BEGIN
  times_all = hour_wvls(i_day) + min_wvls(i_day)/60d + sec_wvls(i_day)/3600d
  iwvs_all = iwvs_all(i_day)
  rains_all = rains_all(i_day)
  azs_all = azs_all(i_day)
  els_all = els_all(i_day)
  sat_type_all = sat_type_all(i_day) 
  sat_nr_all = sat_nr_all(i_day)
  sat_elev_all = sat_elev_all(i_day)
  sat_azi_all = sat_azi_all(i_day)

  n_all = N_ELEMENTS(times_all)

  filename = sat_path +  date+'_h_'+station+'_sat.sav' 
  filename_asc = sat_path +  date+'_h_'+station+'_sat.txt'
; SAVE, filename = filename, times_all, iwvs_all, rains_all, azs_all, els_all, sat_nr_all, sat_elev_all, sat_azi_all

  comment = 'sat_type_all and sat_nr_all are set to dummy values when no unambiguous assignment between .TRK and .BRT file is possible'

;**save to .SAV file
  SAVE, filename = filename, date, times_all, iwvs_all, rains_all, azs_all, els_all, sat_type_all, sat_nr_all, comment

;***save to ascii file
  OPENW, unit_asc, filename_asc, /GET_LUN
  PRINTF, unit_asc, 'Date: ', date, format = '(a6, a6)'
  PRINTF, unit_asc, '|time, UTC|iwv, kgm-2|rain flag|azimuth angle|elevation angle|sat number|'
  FOR k = 0, n_all-1 DO BEGIN
   PRINTF, unit_asc, times_all(k), iwvs_all(k), rains_all(k), azs_all(k), els_all(k), sat_nr_all(k),$
           format = '(f6.2, 1x, f6.2, 1x, i1, 1x, f6.2, 1x, f6.2, 1x, i4)'
  ENDFOR  
  FREE_LUN, unit_asc

 ENDIF ; data of actual day available?

ENDIF ; trk option

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

;----------------------------------------------------------------------
; READ ceilometer data
;----------------------------------------------------------------------

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

time_w = time_w(0l:start_wvl-1l)
time_m = time_m(0l:start_met-1l)

IF start_ceilo GT 1 THEN time_c = time_c(0l:start_ceilo-1l)

rain_all_wvl = rain_all_wvl(0l:start_wvl-1l)

iwv      = iwv_all(0l:start_wvl-1l)
lwp      = lwp_all(0l:start_wvl-1l)
met_all  = met_all(*,0l:start_met-1l)
rain_all = rain_all(0l:start_met-1l)
IF start_ceilo GT 1 THEN lcb_all = lcb_all(0l:start_ceilo-1l)

;.........................................convert time water vapour channels
time_w = time_w/86400d
time_w = time_w + JULDAY(1,1,2001,0,0,0)

xx_N = n_elements(time_w)
xx_dt = time_w[1:xx_N-1] - time_w[0:xx_N-2]
xx_i = where( xx_dt LT 0 , xx_cnt )
if xx_cnt GT 0 then print,'WARNING: unsorted wvl-time !'

i_sort = SORT(time_w)
time_w = time_w(i_sort)
iwv = iwv(i_sort)
lwp = lwp(i_sort)

CALDAT, time_w, MONTH, DAY, YEAR, hour_wvl, min_wvl, sec_wvl
i_day = WHERE(day EQ date_d)

time_w = hour_wvl(i_day) + min_wvl(i_day)/60d + sec_wvl(i_day)/3600d

rain_all_wvl = rain_all_wvl(i_day)
iwv = iwv(i_day)
lwp = lwp(i_day)

;.........................................convert time met channels
time_m = time_m/86400d
time_m = time_m + JULDAY(1,1,2001,0,0,0)

xx_N = n_elements(time_m)
xx_dt = time_m[1:xx_N-1] - time_m[0:xx_N-2]
xx_i = where( xx_dt LT 0 , xx_cnt )
if xx_cnt GT 0 then print,'WARNING: unsorted met-time !'

i_sort = SORT(time_m)
time_m = time_m(i_sort)
met_all = met_all(*, i_sort)
rain_all = rain_all(i_sort)

CALDAT, time_m, MONTH, DAY, YEAR, hour_met, min_met, sec_met
i_day = WHERE(day EQ date_d)

time_m = hour_met(i_day) + min_met(i_day)/60d + sec_met(i_day)/3600d
met_all = met_all(*, i_day)
rain_all = rain_all(i_day)



;read iwv(radiosonde Essen) (only jue)
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


;.........................................Plot IWV/LWP time series
!p.font =0
!p.multi=[0,2,4]
!P.Charsize=1.6
!P.thick = 1.5
SET_PLOT, 'PS'
if verbose ge 1 then print,'creating PS file :'+date+'_h_'+station+'_l1a.ps'
DEVICE, /color, /landscape, filename = plot_path +  date+'_h_'+station+'_l1a.ps' 

; load a 32-color colortable similar to the default Tektronix printer colortable
; The first 9 colors are: Index 0=black, 1=white, 2=red, 3=green, 4=blue, 5=cyan, 6=magenta, 8=orange.
TEK_COLOR

; x_min = [0.17,0.57]
; x_max = [0.50,0.9]
x_min = [0.17 ,0.585]
x_max = [0.485,0.90 ]
;y_min = 0.1 
;y_max = y_min + 0.5
y_min = [0.1, 0.27, 0.44]+0.1
y_max = [0.27, 0.44, 0.58]+0.1

xchar = 1.4

; determine range of iwv plot:

if IWV_max GT 0 THEN BEGIN
  ; old method: absolute scaling
  iwvx = iwv
  ; exclude all negative or those larger than IWV_max 
  iok = WHERE(iwv LT IWV_max AND iwv GT 0. , cnt )
  IF cnt gt 0 THEN iwvx = iwv(iok)
  ; determine range 
  dy = MAX(iwvx) - MIN(iwvx)
  ; plot range is min..max +2*10%
  ymin = Min(iwvx)-dy/10.
  ymax = MAX(iwvx)+dy/10.
ENDIF ELSE BEGIN ; IWV_max LE 0 => semi autmatic scaling
  no_rain = WHERE(rain_all_wvl EQ 0, cnt )
  IF cnt gt 0 THEN iwvx = iwv[no_rain] $
              ELSE iwvx = iwv
  iwv_mean = mean(iwvx)
  iwv_stddev = stddev(iwvx)
  IF IWV_max NE 0 THEN BEGIN ; IW_max < 0 => use range given with |iwv_max|
    ymin = iwv_mean-abs(iwv_max/2.)
    ymax = iwv_mean+abs(iwv_max/2.)
  ENDIF ELSE BEGIN ; IWV_max=0 ; IW_max < 0 => use stddev for the range
    ymin = iwv_mean-2*iwv_stddev
    ymax = iwv_mean+2*iwv_stddev
  ENDELSE
  ; no negative values
  IF ymin < 0 then BEGIN
    ymax = ymax-ymin
    ymin = 0 
  ENDIF ; ymin < 0
  ; offset to plot the rain event bars
  dy = 0.8*(ymax-ymin)
ENDELSE ; IWV_max LE 0
if verbose ge 1 then print,'Range of IWV plot will be:[',ymin,',',ymax,']'

; define range etc. of plot
PLOT,time_w,iwv,Yrange=[ymin,ymax],ystyle=1,$
     Xrange=[0., 24.], Xstyle=1, Xtitle='Time (UTC) [h]',xthick=2, xticklen=-0.08,$
     xcharsize=xchar,ytickv=2,$; title = 'Integrated Water Vapor',$
     position=[x_min(0),y_min(0),x_max(1),y_max(0)],/NODATA

; plot rain periods
index = WHERE(rain_all_wvl EQ 1,count)
IF count GT 0 THEN BEGIN
 FOR j =0l,LONG(count-1l) DO BEGIN
  IF time_w(index(j)) GT 0.05 AND time_w(index(j)) LT 23.95 THEN BEGIN
   x = [time_w(index(j)),time_w(index(j))]
   OPLOT,x,[ymin+dy/90.,ymax-dy/90.],color=3,thick=3
  ENDIF
 ENDFOR
ENDIF

; plot iwv data
OPLOT, time_w, iwv, psym = 3

; triangles for positive outlyers
idx = where( iwv GT ymax , cnt )
IF cnt GT 0 THEN BEGIN
  OPLOT, time_w(idx), replicate(ymax-0.001*(ymax-ymin),cnt), psym = 5
ENDIF 

; plot iwv(radiosonde essen) (only jue)
if station eq 'jue' AND iwv_radio(0) ne -99 then begin
 OPLOT, time_radio, iwv_radio, psym = 2, color = 6
endif



ymin = -50.
ymax = 1000.
dy = ymax - ymin
ymin = ymin-dy/10.
ymax = ymax + dy/10.

PLOT,time_w,lwp,Yrange=[ymin,ymax],ystyle=1,$
     ytickv=2,$;title='Integrated Liquid Water',$
     Xrange=[0.,24.], Xstyle=1, Xtitle='Time (UTC) [h]',$
     position=[x_min(0),y_min(1),x_max(1),y_max(1)],xcharsize = 0.000001, /NODATA

index = WHERE(rain_all_wvl EQ 1,count)
IF count GT 0 THEN BEGIN
 FOR j =0l,LONG(count-1l) DO BEGIN
  IF time_w(index(j)) GT 0.05 AND time_w(index(j)) LT 23.95 THEN BEGIN
   x = [time_w(index(j)),time_w(index(j))]
   OPLOT,x,[ymin+dy/90.,ymax-dy/90.],color=3,thick=3
  ENDIF
 ENDFOR
ENDIF

OPLOT, time_w, lwp, psym = 3
OPLOT, [0., 90000.], [0., 0.], linestyle = 1, thick = 1

ymin = 20.
ymax = 100.
dy = ymax-ymin

; define axis for humidity plot
PLOT,time_m,met_all(2,*),Yrange=[ymin, ymax],ystyle=1,$
     ytickv=2,ytitle='RH [%]',title='Humidity',$
     Xrange=[0.,24.], Xstyle=1, Xtitle='Time (UTC) [h]',$
     position=[x_min(0),0.75,x_max(0),0.95],/NODATA

; plot rain flag as green fat lines so that they give green blocks...
index = WHERE(rain_all EQ 1,count)
IF count GT 0 THEN BEGIN
 FOR i =0l,LONG(count-1l) DO BEGIN
  IF time_m(index(i)) GT 0.05 AND time_m(index(i)) LT 23.95 THEN BEGIN
   x = [time_m(index(i)),time_m(index(i))]
   OPLOT,x,[ymin+dy/90.,ymax-dy/90.],color=3,thick=3
  ENDIF
 ENDFOR
ENDIF

; plot relative humidity
OPLOT, time_m, met_all(2,*), psym = 3

; plot specific humidity or water vapor mixing ratio
; q = mu * e/P = mu * E_sat(T)*RH/100 * / P
; mu = 622 g/kg (mole mass ratio water/air)
; E_sat water vapor pressrure at saturation
; E_sat( T ) = E0*exp(c1*T/(c2+T))   [Magnus formula]
; T in deg Celsius
; ... and absolute water vapor (=vapor density in g/m^3)
; a = q*ro = q * P/(R_air * T_abs)
R_air = 287.04 ; J/(kg K) = 1/100 hPa/(kg K)
E_sat_0 = 6.1078
E_sat_c1 = 17.0809
E_sat_c2 = 234.175
T_null = 273.15
; calc wv saturation pressure
e_sat = E_sat_0 * exp( E_sat_c1 * (met_all(0,*)-T_null) / (E_sat_c2 + (met_all(0,*)-T_null)) )
; specific hmidity
q_spec = 622.* (e_sat * met_all(2,*)/100.) / met_all(1,*)
q_abs = q_spec * met_all(1,*) / ( R_air*0.01 * met_all(0,*) )
;-)
if verbose then begin
  print,'T in ',min(met_all(0,*)),' -- ',max(met_all(0,*)),' K'
  print,'P in ',min(met_all(1,*)),' -- ',max(met_all(1,*)),' hPa'
  print,'RH in',min(met_all(2,*)),' -- ',max(met_all(2,*)),' %'
  print,'esat in',min(e_sat),' -- ',max(e_sat),' hPa'
  print,'q_spec in ',min(q_spec),' -- ',max(q_spec),' g/kg'
  print,'q_abs in ',min(q_abs),' -- ',max(q_abs),' g/m^3'
endif
; define axis
AXIS, yaxis=1, yrange=[0.,20.], ytitle='q [g/kg], a [g/m!e3!N]', yticklen=-!P.ticklen, color=4, /save
; plot q_spec and q_abs ...
OPLOT, time_m, q_spec, color=4
OPLOT, time_m, q_abs , color=4, linestyle=2

; plot temperature
; determine range +2*10%
dy   = MAX(met_all(0,*)) - MIN(met_all(0,*))
ymin = MIN(met_all(0,*)) - dy/10.
ymax = MAX(met_all(0,*)) + dy/10.
; define axes
PLOT,time_m,met_all(0,*),Yrange=[ymin, ymax],ystyle=1,$
     ytickv=2,ytitle='T [K]',title='Air Temperature',$
     Xrange=[0.,24.], Xstyle=1, Xtitle='Time (UTC) [h]',$
     position=[x_min(1),0.75,x_max(1),0.95],/NODATA
; mark rain events
IF count GT 0 THEN BEGIN
 FOR i =0l,LONG(count-1) DO BEGIN
  IF time_m(index(i)) GT 0.05 AND time_m(index(i)) LT 23.95 THEN BEGIN
   x = [time_m(index(i)),time_m(index(i))]
   OPLOT,x,[ymin+dy/90.,ymax-dy/90.],color=3,thick=3
  ENDIF
 ENDFOR
ENDIF
; do the temperature plot
OPLOT, time_m, met_all(0,*), psym = 3
; an axis with deg C tickmarks on the right side
AXIS, yaxis=1, yrange=[ymin-T_null , ymax-T_null], ystyle=1, ytitle='T [!Uo!NC]', yticklen=-!P.ticklen, /save


IF ceilo_type EQ 'dwd' THEN ytitc = 'LCB(DWD) [km]'
IF ceilo_type EQ 'igm' THEN ytitc = 'LCB(IGM) [km]'

IF start_ceilo GT 1 THEN BEGIN
 i99 = WHERE(lcb_all GE 0.)
 IF i99(0) NE -1 THEN lcb_all(i99) = lcb_all(i99)/1000.
 PLOT, time_c, lcb_all, yrange=[0., 7.], ystyle = 1,$
       xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]',$
       xcharsize = 0.00001, ytickv = 2, ytitle = ytitc, $
       position=[x_min(0),y_min(2),x_max(1),y_max(2)], /NODATA, /NOERASE

 OPLOT, time_c, lcb_all, color = 4, thick = 2, psym = 4, symsize = 0.3
ENDIF


XYOUTS,0.532,0.054,station+', '+date,alignment=0.5,charsize=1.2,/NORMAl
XYOUTS,0.532,0.008,'Rain events - blocked green',alignment=0.5, color=3,charsize=0.7,/NORMAl
XYOUTS,0.1,0.25,'IWV [kg m!e-2!n]',charsize = 0.8, alignment=0.5,orientation=90.,/NORMAl
XYOUTS,0.1,0.5,'LWP [g m!e-2!n]',charsize = 0.8, alignment=0.5,orientation=90.,/NORMAl
IF station eq 'jue' then XYOUTS, 0.532, 0.028,'IWV Radiosonde Essen',alignment=0.5, color=6 , charsize=0.7,/NORMAL
DEVICE,/CLOSE

END
