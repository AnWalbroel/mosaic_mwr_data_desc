;+
;****************************
PRO TB_OFFSET_CORRECTION_NDS;,$
;****************************
;INPUT:
;OUPUT:
; $Id: $
; Abstract:
; * calculation of TB offset correction using COSMO-DE, ICON, AROME or Cloudnet output
; * This program uses the parameter file par_tb_offset_correction.pro.
; Authors:
; U. Loehnert
; Date:
; 2012-09-06
; Dependencies:
; -
; Changes:
; 2017-12-12 (UL):
; update to include ICON or AROME output for offset correction
; 2018-09-03 (UL):
; include new directory structure
;-

@par_tb_offset_correction

FOR i_off = 0, n_cal-1 DO BEGIN ; loop over calibration periods

 IF plot_only EQ 1 THEN BEGIN
  RESTORE, file = file_save(i_off)
  GOTO, PL_ONLY
 ENDIF

;*make n_cs_user hour time array
 t_cs_user = FINDGEN((24./n_cs_user)+1)*n_cs_user
 n_t_cs_user = N_ELEMENTS(t_cs_user)

 two_min = 2./60.

 big = 5000.
 iii = 0

;****define parameters to compare
 TB_meas_mean = REPLICATE(-999., n_freq, n_ang, big)
 TB_meas_std = REPLICATE(-999., n_freq, n_ang, big)
 TB_calc = REPLICATE(-999., n_freq, n_ang, big)
 TB_diff = REPLICATE(-999., n_freq, n_ang, big)
 jd_diff = REPLICATE(-999., big)
 iwv_rs = REPLICATE(-999., big)
 time_rs = REPLICATE(-999., big)
 iwv_gps = REPLICATE(-999., big)
 time_gps = REPLICATE(-999., big)
 iwv_rad = REPLICATE(-999., big)
 iwv_mod = REPLICATE(-999., big)

;****find level2a clear sky data within calibration periods

;**number of days in calibration period

 n_days = LONG(e_date_jd_offset_new(i_off+1) - s_date_jd_offset_new(i_off)) + 1l

;**create array of dates
 jd_dates = LINDGEN(n_days) + LONG(s_date_jd_offset_new(i_off)) + 0.5d
 CALDAT, jd_dates, mm, dd, yy, hh, min, ss 
 yyyy_c_a = STRING(yy, format = '(i4)')
 mm_c_a = REPLICATE('', n_days)
 dd_c_a = REPLICATE('', n_days)

 FOR i = 0, n_days-1 DO BEGIN
  IF dd(i) LT 10 THEN BEGIN
   dd_c_a(i) = '0'+STRING(dd(i), format = '(i1)')
  ENDIF ELSE BEGIN
   dd_c_a(i) = STRING(dd(i), format = '(i2)') 
  ENDELSE
  IF mm(i) LT 10 THEN BEGIN
   mm_c_a(i) = '0'+STRING(mm(i), format = '(i1)')
  ENDIF ELSE BEGIN
   mm_c_a(i) = STRING(mm(i), format = '(i2)') 
  ENDELSE
 ENDFOR

 months_a = STRMID(yyyy_c_a, 2, 2) + mm_c_a
 dates_a = STRMID(yyyy_c_a, 2, 2) + mm_c_a + dd_c_a

;**create array of file descriptors
 dirs = yyyy_c_a + '/' + mm_c_a + '/' + dd_c_a + '/'
 file_names_clwvi_a = FILE_SEARCH(path_l2 + dirs + '*clwvi*.nc')
 n_file_names_a = N_ELEMENTS(file_names_clwvi_a)

 n_dirs = N_ELEMENTS(dirs)
 IF file_names_clwvi_a(0) EQ '' THEN BEGIN
  print, '> no lwp file found from '
  print, dirs(0)
  print, '> to '
  print, dirs(n_dirs-1) 
  GOTO, SKIP_CAL
 ENDIF 

;***read lwp files and check for clear sky using lwp_thres

 FOR i = 0, n_days-1 DO BEGIN

  yyyy_c = yyyy_c_a(i)
  mm_c = mm_c_a(i)
  dd_c = dd_c_a(i)
  
  dates = yyyy_c_a(i) + mm_c_a(i) + dd_c_a(i)

  print, 'working on: ', dates

  i_exclude = WHERE(date_exclude EQ dates)

  IF i_exclude(0) NE -1 THEN BEGIN
   print, 'excluding ' + dates + ' manually...'
   GOTO, SKIP_DAY
  ENDIF

;**check if TB data is available in the first place

  IF data_type EQ 'mwrBL' THEN BEGIN
   file_mwrBL = FILE_SEARCH(path_l1 + dirs(i) + '*mwrBL*.nc')
  
   IF file_mwrBL(0) EQ '' THEN BEGIN
    print, 'no matching level mwrBL data found!'
    GOTO, SKIP_DAY 
   ENDIF
  ENDIF ELSE IF data_type EQ 'mwr' THEN BEGIN
   file_mwr = FILE_SEARCH(path_l1 + months_a(i) + dirs(i) + '*mwr*.nc')
  
   IF file_mwr(0) EQ '' THEN BEGIN
    print, 'no matching level1 mwr data found!'
    GOTO, SKIP_DAY 
   ENDIF
  ENDIF 

  file_name_clwvi = FILE_SEARCH(path_l2 + dirs(i) + '*clwvi*.nc')   
  file_name_prw = FILE_SEARCH(path_l2 + dirs(i) + '*prw*.nc')

  IF file_name_clwvi(0) EQ '' THEN BEGIN 
   print, '> no lwp file found for ' + dates
   GOTO, SKIP_DAY
  ENDIF

;  READ_LEVEL2A_NC, file_name, algo, comment, time_2a, ele, azi, iwv, lwp,$
;                   flag_2a, temp, pres, relh

  READ_CLWVI_NC, file_name_clwvi, time_clwvi, ele, azi, lwp, lwp_org, lwp_oz, flag_2a_clwvi
  SEC2DATE1970, time_clwvi, x, time_2a_clwvi

  READ_PRW_NC, file_name_prw, time_prw, ele, azi, prw, prw_oz, flag_2a_prw
  SEC2DATE1970, time_prw, x, time_2a_prw


;**loop over n_cs_user intervals and find out if they are clear or cloudy
  opened = 0
  FOR ii = 0, n_t_cs_user-1 DO BEGIN

;*loop over 2-minute intervals within each t_cs_user interval and calculate lwp std 
   clear = 0
   cloudy = 0

   FOR j = 0, (n_cs_user*30)-2 DO BEGIN
    i_two = WHERE(time_2a_clwvi GE t_cs_user(ii)+FLOAT(j)*two_min AND time_2a_clwvi LT t_cs_user(ii)+FLOAT(j+1)*two_min AND ele GT 89.4 AND ele LT 90.6 AND lwp NE -999., nn)
    i_flag = WHERE(time_2a_clwvi GE t_cs_user(ii)+FLOAT(j)*two_min AND time_2a_clwvi LT t_cs_user(ii)+FLOAT(j+1)*two_min AND (flag_2a_clwvi GT 0. AND flag_2a_clwvi LT 1024.)) 
    IF i_flag(0) NE -1 THEN cloudy = 1

    IF nn GE 30 THEN BEGIN
     lwp_std = STDDEV(lwp(i_two))
     IF lwp_std GE lwp_std_thres/1000. THEN BEGIN
      cloudy = 1
     ENDIF ELSE BEGIN
      clear = 1
;print, nn, t_cs_user(ii)+FLOAT(j)*two_min, lwp_std*1000.

     ENDELSE

    ENDIF

   ENDFOR ; loop over 2-minute intervals 

;print, file_names_clwvi(i), t_cs_user(ii), clear, cloudy

;**search for model output in case of clear-sky cases
   IF clear EQ 1 AND cloudy EQ 0 THEN BEGIN
    print, 'clear sky case detected on: ', yyyy_c + mm_c + dd_c
    date_d = DOUBLE(yyyy_c) + DOUBLE(mm_c) + DOUBLE(dd_c)

    print, 'at time (UTC): ', t_cs_user(ii)

    file_name_tar =''

    IF model_choice EQ 'cosmo' THEN BEGIN
     IF date_d LE 20140529. THEN BEGIN   
      file_name_tar = model_path + yyyy_c + '/' + yyyy_c + mm_c + dd_c + '_gop9_lmk_cols.tar'
      file_name_m  = yyyy_c + mm_c + dd_c + '_gop9_lmk_' + station_str + '.nc' 

      tag = 0
      check = FILE_SEARCH(file_name_m)
      IF check(0) EQ '' THEN BEGIN
       check = FILE_SEARCH(file_name_tar)
       IF check(0) EQ '' THEN BEGIN
        print, 'No model profile available for ' + yyyy_c + mm_c + dd_c
        GOTO, SKIP_DAY
       ENDIF ELSE BEGIN
        print, 'copying file ' + file_name_tar
        SPAWN, 'cp ' + file_name_tar + ' .'
        file_name_tar = yyyy_c + mm_c + dd_c + '_gop9_lmk_cols.tar'
        SPAWN, 'tar -xvf ' + file_name_tar + ' ' + file_name_m + '.gz' 
        SPAWN, 'gunzip ' + file_name_m + '.gz'
        tag = 1
       ENDELSE
      ENDIF

      IF tag EQ 1 THEN BEGIN
       SPAWN, 'rm *.tar' 
       SPAWN, 'rm *.nc'
      ENDIF
  
     ENDIF ELSE BEGIN ; model_choice = 'cosmo' and date GT 20140529

      IF site EQ 'joy' OR site EQ 'rao' THEN BEGIN
       file_name_m = model_path + yyyy_c + '/sitesdomain_7x7/' + yyyy_c + mm_c + dd_c + '_lmk_cols7x7_ana.nc'
       i_stat = 24  
      ENDIF ELSE IF site EQ 'cab' OR site EQ 'ham' or site EQ 'pay' OR site EQ 'ufs' THEN BEGIN 
       file_name_m = model_path + yyyy_c + '/sitesdomain_3x3/' + yyyy_c + mm_c + dd_c + '_lmk_cols3x3_ana.nc' 
       i_stat = 4
      ENDIF ELSE BEGIN
       file_name_m = model_path + yyyy_c + '/sitesdomain_1x1/' + yyyy_c + mm_c + dd_c + '_lmk_cols1x1_ana.nc' 
       i_stat = 0
      ENDELSE
     
     ENDELSE     
  
    ENDIF

    IF model_choice EQ 'arome' THEN BEGIN
     file_name_m = model_path + '*' + yyyy_c + mm_c + dd_c + '*.nc'
     print, 'Sorry, AROME option for tb offset correction not yet implemented'
     STOP
    ENDIF

    IF model_choice EQ 'icon' THEN BEGIN
     file_name_m = model_path + 'yyyy_c' + '/' + mm_c + '/' + yyyy_c + mm_c + dd_c + '*.nc'
     print, 'Sorry, ICON option for tb offset correction not yet implemented'
     STOP
    ENDIF

    IF model_choice EQ 'cloudnet' THEN BEGIN
     file_name_m = model_path + yyyy_c + '/' + yyyy_c + mm_c + dd_c + '_' + stat_cloudnet + '_categorize.nc'     
    ENDIF

;STOP fuer COSMO Sauelen - noch nicht fertig!

;**read model output
    IF opened EQ 0 AND model_choice EQ 'cosmo' AND date_d LE 20140529 THEN BEGIN
     infid = NCDF_OPEN(file_name_m, /NOWRITE)
 
     time = NCDF_VARID(infid, 'time1h')
     longitude = NCDF_VARID(infid, 'longitude')
     latitude = NCDF_VARID(infid, 'latitude')
     z = NCDF_VARID(infid, 'hfl')
     T = NCDF_VARID(infid, 'temperature')
     p = NCDF_VARID(infid, 'p')
     q = NCDF_VARID(infid, 'qv')

     NCDF_VARGET, infid, time, time
     NCDF_VARGET, infid, longitude, longitude
     NCDF_VARGET, infid, latitude, latitude
     NCDF_VARGET, infid, z, z
     NCDF_VARGET, infid, T, T
     NCDF_VARGET, infid, p, p
     NCDF_VARGET, infid, q, q

     SEC2DATE1970, time, x, time_mod
 
     opened = 1

    ENDIF ELSE IF opened EQ 0 AND model_choice EQ 'cosmo' AND date_d GE 20140530 THEN BEGIN 

     infid = NCDF_OPEN(file_name_m, /NOWRITE)
 
     time = NCDF_VARID(infid, 'TIME')
     coordinates = NCDF_VARID(infid, 'STATION_COORD_GEO')
     T = NCDF_VARID(infid, 'T') ; K
     p = NCDF_VARID(infid, 'P') ; Pa
     q = NCDF_VARID(infid, 'QV') ; kg/kg

     NCDF_VARGET, infid, time, time
     NCDF_VARGET, infid, coordinates, coordinates
     NCDF_VARGET, infid, T, T
     NCDF_VARGET, infid, p, p
     NCDF_VARGET, infid, q, q

     SEC2DATE1970, time, x, time_mod

;*find matching station in model output
     ii_lat = WHERE(coordinates(*, 0) EQ site_lat)
     ii_lon = WHERE(coorindates(0, *) EQ site_lon) 

;**match model times to measurement times

     imatch = WHERE(time_mod GE t_cs_user(ii) AND time_mod LT t_cs_user(ii+1))
 
;     IF imatch(0) NE -1 THEN BEGIN

;*derive height array above ground z in m
;      nz = N_ELEMENTS(p(imatch, *, 0, ii))
;      z = REPLICATE(0., nz+1) 
;      Rl = 287. ; J / kg * K
;      g = 9.81 ; m / s * s

;      FOR j = 0, nz-2 DO BEGIN            
;       Tm = (T(j+1)+T(j))/2.
;       z(j+1) = z(j) -1*ALOGN(p(j+1)/p(j))*(Rl)
;      ENDFOR 
 
;      f = 0 ;time since model forecast
;      c = [19, 20, 25, 26] ;averaged columns

;      z1 = AVERAGE(z(*, c), 2)
;      T1 = REFORM(AVERAGE(T(*, f, *, c), 4))
;      p1 = REFORM(AVERAGE(p(*, f, *, c), 4))
;      q1 = REFORM(AVERAGE(q(*, f, *, c), 4))


;      z2 = REVERSE(z1)
;      T2 = REVERSE(AVERAGE(T1(imatch, *), 1))
;      p2 = REVERSE(AVERAGE(p1(imatch, *), 1)) 
;      q2 = REVERSE(AVERAGE(q1(imatch, *), 1))

;      ENDIF

;*derive height array above ground z in m
;     ii = WHERE(coordinates        

;     opened = 1

;* derive height above ground z in m
;

    ENDIF ELSE IF opened EQ 0 AND model_choice EQ 'arome' THEN BEGIN

;HIER WEITER falls AROME Daten einfliessen sollen

    ENDIF ELSE IF opened EQ 0 AND model_choice EQ 'icon' THEN BEGIN

;HIER WEITER falls ICON Daten einfliessen 

    ENDIF ELSE IF opened EQ 0 AND model_choice EQ 'cloudnet' THEN BEGIN 

     file_exist = FILE_SEARCH(file_name_m, COUNT=n_exist)
     IF n_exist EQ 0 THEN BEGIN
      print, '> no cloudnet data found'
      GOTO, SKIP_DAY
     ENDIF

     infid = NCDF_OPEN(file_name_m, /NOWRITE)
     time_mod = NCDF_VARID(infid, 'time')
     z_msl = NCDF_VARID(infid, 'altitude')
     z = NCDF_VARID(infid, 'model_height') ; m above MSL
     T = NCDF_VARID(infid, 'temperature') ; K
     p = NCDF_VARID(infid, 'pressure') ; Pa
     q = NCDF_VARID(infid, 'specific_humidity') ; 1

     NCDF_VARGET, infid, time_mod, time_mod
     NCDF_VARGET, infid, z_msl, z_msl
     NCDF_VARGET, infid, z, z
;     NCDF_ATTGET, infid, z_at, att_name, att_value
     NCDF_VARGET, infid, T, T
     NCDF_VARGET, infid, p, p
     NCDF_VARGET, infid, q, q

;**match model times to measurement times

     imatch = WHERE(time_mod GE t_cs_user(ii) AND time_mod LT t_cs_user(ii+1))

     z2 = z
     nz = N_ELEMENTS(z2)
     T2 = REPLICATE(-999., nz)
     p2 = REPLICATE(-999., nz)
     q2 = REPLICATE(-999., nz)
     FOR k = 0, nz-1 DO BEGIN
      T2(k) = MEAN(T(k, imatch))
      p2(k) = MEAN(p(k, imatch))
      q2(k) = MEAN(q(k, imatch))
     ENDFOR  
    
    ENDIF; cloudnet data use

;*extend model profile to 30 km
    mmx = FIX(mm_c)
    z_model_top = z2(N_ELEMENTS(z2)-1)
    i_low = WHERE(z_model_top LT (z_final_m+z2(0)-z_sup))

    IF i_low(0) NE -1 THEN BEGIN
     z2 = [z2, (z_final_m(i_low)+z2(0)-z_sup)]
     n_low = N_ELEMENTS(i_low)
     a_null = REPLICATE(0., n_low)
     T2 = [T2, T_mean(i_low, mmx-1)]
     p2 = [p2, p_mean(i_low, mmx-1)]
     q2 = [q2, a_null]
    ENDIF

    e2 = (q2*p2)/0.622
    a2 = e2/(462.*T2)

    IWV_CALC, a2, z2, iwv2
    iwv_mod(iii) = iwv2

    i_match_iwv = WHERE(time_2a_prw GT t_cs_user(ii) AND time_2a_prw LT t_cs_user(ii+1) AND prw GT 0. AND flag_2a_prw EQ 0 AND ele GT 89.4 AND ele LT 90.6)
    IF i_match_iwv(0) NE -1 THEN BEGIN
     iwv_rad(iii) = MEAN(prw(i_match_iwv))
     IF ABS(iwv_rad(iii)-iwv2) GE iwv_thres THEN BEGIN      
      print, 'Model minus MWR IWV (in kgm-2) larger than threshold ', iwv_thres
      print, '= ', iwv2-iwv_match
      print, 'Aborting period'
      GOTO, ABORT
     ENDIF
    ENDIF

;*find matching radiosonde
    rs_files = FILE_SEARCH(rs_path + yyyy_c + '/' + mm_c + '/20' + dates + '*.html')
    n_rs_files = N_ELEMENTS(rs_files)
 
    IF rs_files(0) NE '' THEN BEGIN      
     FOR irs = 0, n_rs_files-1 DO BEGIN
 
      READ_RS_WYOMING, rs_files(irs), time_rs_jd, iwv_rs_x
      CALDAT, time_rs_jd, month_rs, day_rs, year_rs, hour_rs, sec_rs

      imatch = WHERE(hour_rs GE t_cs_user(ii) AND hour_rs LT t_cs_user(ii+1) AND iwv_rs_x GT 0.)

      IF imatch(0) NE -1 THEN BEGIN
       iwv_rs(iii) = iwv_rs_x
       time_rs(iii) = time_rs_jd      
       print, 'found clear RS for ' + dates
       GOTO, SKIP_RS
      ENDIF
 
     ENDFOR
    ENDIF

SKIP_RS:

;*find matching GPS
    gps_files = FILE_SEARCH(gps_path + yyyy_c + '/' + '*' + dates + '*')
    n_gps_files = N_ELEMENTS(gps_files)

    IF gps_files(0) NE '' THEN BEGIN
     FOR igps = 0, n_gps_files-1 DO BEGIN

      infid = NCDF_OPEN(gps_files(igps), /NOWRITE)

      time_gps_ = NCDF_VARID(infid, 'time')
      iwv_gps_ = NCDF_VARID(infid, 'pwv')

      NCDF_VARGET, infid, time_gps_, time_gps_x
      NCDF_VARGET, infid, iwv_gps_, iwv_gps_x

      NCDF_CLOSE, infid

      imatch = WHERE(time_gps_x GE t_cs_user(ii) AND time_gps_x LT t_cs_user(ii+1) AND iwv_gps_x GT 0.)

      IF imatch(0) NE -1 THEN BEGIN
       iwv_gps(iii) = MEAN(iwv_gps_x(imatch))
       time_gps(iii) = MEAN(time_gps_x(imatch))
       print, 'GPS data found for ' + dates
      ENDIF

     ENDFOR
    ENDIF

;**derive offsets from mwr or mwrBL data

    IF data_type EQ 'mwrBL' THEN BEGIN

;*average MWR data

     READ_LEVEL1_MWRBL_NC, file_mwrBL(0), commment, lat, lon, alt, angles, freq, time_1c, tb_1c, flag_1c, temp, pres, relh
     SEC2DATE1970, time_1c, x, time_1c

     imatch = WHERE(time_1c GE t_cs_user(ii) AND time_1c LT t_cs_user(ii+1), n_match) 
     IF n_match GE 2 THEN BEGIN
      FOR k = 0, n_freq-1 DO BEGIN
       FOR kk = 0, n_ang-1 DO BEGIN
        TB_meas_mean(k, kk, iii) = MEAN(tb_1c(k, kk, imatch))   
        TB_meas_std(k, kk, iii) = STDEV(tb_1c(k, kk, imatch))
       ENDFOR
      ENDFOR
     ENDIF

; HIER WEITER
stop
;*RT
     jk = WHERE(FINITE(z2) AND FINITE(T2) AND FINITE(p2) AND FINITE(a2))

     IF N_ELEMENTS(jk) EQ -1 THEN GOTO, ABORT  
     IF N_ELEMENTS(angles) EQ 0 THEN GOTO, ABORT

     print, 'calculating STP for ', dates
     FOR k = 0, n_ang-1 DO BEGIN
      STP, z2(jk), T2(jk), p2(jk), a2(jk), angles(k)-90., freq,$
           TB_calcx, tau, tau_wv, tau_o2, tau_liq, T_mr,$
           air_corr = 'rueeger_avai_02', z_site = z_site          
 
      icheck = WHERE(TB_meas_mean(*, k, iii) LT 0. OR TB_calcx LT 0.)
      IF icheck(0) EQ -1 THEN BEGIN
       TB_calc(*, k, iii) = TB_calcx
       TB_diff(*, k, iii) = TB_meas_mean(*, k, iii) - TB_calcx
       DATE_TIME_TO_JULDAT, yyyy_c + mm_c + dd_c, t_cs_user(ii), jd  
       jd_diff(iii) = jd 
      ENDIF 

     ENDFOR
    
    ENDIF ELSE IF data_type EQ 'mwr' THEN BEGIN

;*average MWR data
     FOR kk = 0, n_ang-1 DO BEGIN
      READ_LEVEL1_MWR_NC, file_mwr(0), comment, time_1b, freq, tb_1b, wavel_ir, tb_ir, ele_ir, flag_1b, temp, pres, relh, az, el
      SEC2DATE1970, time_1b, x, time_1b

      imatch = WHERE(time_1b GE t_cs_user(ii) AND time_1b LT t_cs_user(ii+1)$
                     AND el LT angles(kk)+0.6 AND el GT angles(kk)-0.6, n_match) 

      IF n_match GE 2 THEN BEGIN
       FOR k = 0, n_freq-1 DO BEGIN
        TB_meas_mean(k, kk, iii) = MEAN(tb_1b(k, imatch))   
        TB_meas_std(k, kk, iii) = STDEV(tb_1b(k, imatch))
       ENDFOR
      ENDIF
     ENDFOR

;*RT
     jk = WHERE(FINITE(z2) AND FINITE(T2) AND FINITE(p2) AND FINITE(a2))

     IF N_ELEMENTS(jk) EQ -1 THEN GOTO, ABORT  
     IF N_ELEMENTS(angles) EQ 0 THEN GOTO, ABORT

     print, 'calculating STP for ', dates
     FOR k = 0, n_ang-1 DO BEGIN
      STP, z2(jk), T2(jk), p2(jk), a2(jk), angles(k)-90., freq,$
           TB_calcx, tau, tau_wv, tau_o2, tau_liq, T_mr,$
           air_corr = 'rueeger_avai_02', z_site = z_site          
 
      icheck = WHERE(TB_meas_mean(*, k, iii) LT 0. OR TB_calcx LT 0.)
      IF icheck(0) EQ -1 THEN BEGIN
       TB_calc(*, k, iii) = TB_calcx
       TB_diff(*, k, iii) = TB_meas_mean(*, k, iii) - TB_calcx
       DATE_TIME_TO_JULDAT, yyyy_c + mm_c + dd_c, t_cs_user(ii), jd  
       jd_diff(iii) = jd 
      ENDIF 

     ENDFOR
    
    ENDIF ELSE BEGIN; brt OR blb

     print, 'You must specifiy either brt or blb as data_type!' 
     stop

    ENDELSE 

   ENDIF; case clear/cloudy

ABORT:
   DATE_TIME_TO_JULDAT, yyyy_c + mm_c + dd_c, t_cs_user(ii), jd  
   jd_diff(iii) = jd 
   iii = iii + 1

;HIER WEITER: test was am 5.7. 21 bzw. am 29.7. 0 los ist - überhaupt noch drin nach erneutem Krähenflaggen am 5.2.15??

  ENDFOR ; loop over n_ts_user intervals ; ii

SKIP_DAY:

 ENDFOR ; loop over daily files ; i

 aa = WHERE(TB_calc(0, 0, *) GT 0. AND TB_meas_mean(0, 0, *) GT 0., n_aa)

 IF aa(0) NE -1 THEN BEGIN
  TB_calc = TB_calc(*, *, aa)
  TB_meas_mean = TB_meas_mean(*, *, aa)
  TB_meas_std = TB_meas_std(*, *, aa)
  TB_diff = TB_diff(*, *, aa)
  jd_diff = jd_diff(aa)
  iwv_rs = iwv_rs(aa)
  time_rs = time_rs(aa) 
  iwv_gps = iwv_gps(aa)
  time_gps = time_gps(aa) 
  iwv_rad = iwv_rad(aa)
  iwv_mod = iwv_mod(aa)

  IF N_ELEMENTS(aa) GT 1 THEN BEGIN
   TB_calc_m = AVERAGE(TB_calc(*, *, *), 3)
   TB_meas_mean_m = AVERAGE(TB_meas_mean(*, *, *), 3)
   TB_meas_std_m = AVERAGE(TB_meas_std(*, *, *), 3)
   TB_diff_m = AVERAGE(TB_diff(*, *, *), 3) 
  ENDIF ELSE BEGIN
   TB_calc_m = REFORM(TB_calc(*, *, 0:0))
   TB_meas_mean_m = REFORM(TB_meas_mean(*, *, 0:0))
   TB_meas_std_m = REFORM(TB_meas_std(*, *, 0:0))
   TB_diff_m = REFORM(TB_diff(*, *, 0:0))
  ENDELSE

;****save offset data to file
  DATE_TIME_TO_JULDAT, start_date(i_off), start_time(i_off), jd_start_date
  DATE_TIME_TO_JULDAT, end_date(i_off), end_time(i_off), jd_end_date

  jd_start_date_oc = jd_start_date
  jd_end_date_oc = jd_end_date
  freq_oc = freq
  angles_oc = angles
  lwp_std_thres_oc = lwp_std_thres
  TB_diff_oc = TB_diff
  jd_diff_oc = jd_diff
  TB_diff_m_oc = TB_diff_m
  iwv_rs_oc = iwv_rs
  iwv_gps_oc = iwv_gps
  iwv_mod_oc = iwv_mod
  iwv_rad_oc = iwv_rad
  TB_calc_oc = Tb_calc_m
  TB_meas_mean_oc = Tb_meas_mean_m
  TB_meas_std_oc = Tb_meas_std_m

  SAVE, jd_start_date_oc, jd_end_date_oc, freq_oc, angles_oc, lwp_std_thres_oc, TB_calc_oc, TB_meas_mean_oc,$ 
        TB_meas_std_oc, TB_diff_oc, jd_diff_oc, TB_diff_m_oc,$
        iwv_rs_oc, iwv_gps_oc_oc, iwv_rad_oc, iwv_mod_oc, filename=file_save(i_off)

PL_ONLY:

  TB_diff_std = REPLICATE(-999., n_freq, n_ang)
  n_aa = N_ELEMENTS(iwv_rs_oc)

;****plot OFFSETs as a function of Frequency
  ik = WHERE(freq_oc GT 20 AND freq_oc LT 35)
  iv = WHERE(freq_oc GT 50 AND freq_oc LT 60)

  SET_PLOT, 'ps'
  LOADCOL, 'col1'
  DEVICE, /color, file=ps_file(i_off)
  !p.multi = [0, 1, 1]
 
;***K-Band (means)
  max_TB_k = -999.
  min_TB_k = 999.

  FOR i = 0, N_ELEMENTS(TB_diff_m_oc(ik, 0))-1 DO BEGIN
   FOR j = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
    IF TB_diff_m_oc(ik(i), j) GT max_TB_k THEN max_TB_k = TB_diff_m_oc(ik(i), j)  
    IF TB_diff_m_oc(ik(i), j) LT min_TB_k THEN min_TB_k = TB_diff_m_oc(ik(i), j)
   ENDFOR
  ENDFOR

  PLOT, freq_oc(ik), TB_diff_m_oc(ik, 0), yrange = [min_TB_k, max_TB_k],$
        xtitle='GHz', ytitle='TB_meas - TB_calc [K]', title=tit(i_off), /NODATA

  FOR i = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
   OPLOT, freq_oc(ik), TB_diff_m_oc(ik, i), color = i, thick = 3, psym = -2 
   XYOUTS, 0.15, 0.9-0.07*FLOAT(i), STRING(angles_oc(i), format = '(f4.1)'), /normal, color = i
  ENDFOR
  OPLOT, [0., 100.], [0., 0.],linestyle = 1
  XYOUTS, 0.15, 0.2, 'N = ' + STRING(n_aa, format = '(i3)'), /normal

;***K-Band (variances)
  IF n_aa GT 2 THEN BEGIN
   max_TB_k = -999.
   min_TB_k = 999.

   FOR i = 0, N_ELEMENTS(TB_diff_m_oc(ik, 0))-1 DO BEGIN
    FOR j = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
     IF STDEV(TB_diff_oc(ik(i), j, *)) GT max_TB_k THEN max_TB_k = STDEV(TB_diff_oc(ik(i), j, *))  
     IF STDEV(TB_diff_oc(ik(i), j, *)) LT min_TB_k THEN min_TB_k = STDEV(TB_diff_oc(ik(i), j, *))
     TB_diff_std(ik(i), j) = STDEV(TB_diff_oc(ik(i), j, *))
    ENDFOR
   ENDFOR

   PLOT, freq_oc(ik), TB_diff_m_oc(ik, 0), yrange = [min_TB_k, max_TB_k],$
         xtitle='GHz', ytitle='STD(TB_meas-TB_calc) [K]', title=tit(i_off), /NODATA

   FOR i = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
    OPLOT, freq_oc(ik), TB_diff_std(ik, i), color = i, thick = 3, psym = -2
    XYOUTS, 0.15, 0.9-0.07*FLOAT(i), STRING(angles_oc(i), format = '(f4.1)'), /normal, color = i
   ENDFOR
   OPLOT, [0., 100.], [0., 0.],linestyle = 1
   XYOUTS, 0.15, 0.2, 'N = ' + STRING(n_aa, format = '(i3)'), /normal
  ENDIF

;***V-Band (means)
  max_TB_v = -999.
  min_TB_v = 999.

  FOR i = 0, N_ELEMENTS(TB_diff_m_oc(iv, 0))-1 DO BEGIN
   FOR j = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
    IF TB_diff_m_oc(iv(i), j) GT max_TB_v THEN max_TB_v = TB_diff_m_oc(iv(i), j)  
    IF TB_diff_m_oc(iv(i), j) LT min_TB_v THEN min_TB_v = TB_diff_m_oc(iv(i), j)
   ENDFOR
  ENDFOR

  PLOT, freq_oc(iv), TB_diff_m_oc(iv, 0), yrange = [min_TB_v, max_TB_v],$
        xtitle='GHz', ytitle='TB_meas - TB_calc [K]', title=tit(i_off), /NODATA

  FOR i = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
   OPLOT, freq_oc(iv), TB_diff_m_oc(iv, i), color = i, thick = 3, psym = -2
   XYOUTS, 0.15, 0.9-0.07*FLOAT(i), STRING(angles_oc(i), format = '(f4.1)'), /normal, color = i
  ENDFOR
  OPLOT, [0., 100.], [0., 0.],linestyle = 1
  XYOUTS, 0.15, 0.2, 'N = ' + STRING(n_aa, format = '(i3)'), /normal

;***V-Band (variances)
  IF n_aa GT 2 THEN BEGIN
   max_TB_v = -999.
   min_TB_v = 999.

   FOR i = 0, N_ELEMENTS(TB_diff_m_oc(iv, 0))-1 DO BEGIN
    FOR j = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
     IF STDEV(TB_diff_oc(iv(i), j, *)) GT max_TB_v THEN max_TB_v = STDEV(TB_diff_oc(iv(i), j, *))  
     IF STDEV(TB_diff_oc(iv(i), j, *)) LT min_TB_v THEN min_TB_v = STDEV(TB_diff_oc(iv(i), j, *))
     TB_diff_std(iv(i), j) = STDEV(TB_diff_oc(iv(i), j, *))
    ENDFOR
   ENDFOR

   PLOT, freq_oc(iv), TB_diff_m_oc(iv, 0), yrange = [min_TB_v, max_TB_v],$
         xtitle='GHz', ytitle='STD(TB_meas-TB_calc) [K]', title=tit(i_off), /NODATA
 
   FOR i = 0, N_ELEMENTS(TB_diff_m_oc(0, *))-1 DO BEGIN
    OPLOT, freq_oc(iv), TB_diff_std(iv, i), color = i, thick = 3, psym = -2
    XYOUTS, 0.15, 0.9-0.07*FLOAT(i), STRING(angles_oc(i), format = '(f4.1)'), /normal, color = i
   ENDFOR
   OPLOT, [0., 100.], [0., 0.],linestyle = 1
   XYOUTS, 0.15, 0.2, 'N = ' + STRING(n_aa, format = '(i3)'), /normal
  ENDIF

;***TB offset time series at 90deg elevation
  jj = WHERE(TB_diff_oc(0, 0, *) NE -999.) 

  dummyy = LABEL_DATE(DATE_FORMAT=['%D-%M','%Y'])
 
  min1 = 999.
  max1 = -999.
  FOR k = 0, n_freq-1 DO BEGIN
   IF MIN(TB_diff_oc(k, 0, jj)) LT min1 THEN min1 = MIN(TB_diff_oc(k, 0, jj))
   IF MAX(TB_diff_oc(k, 0, jj)) GT max1 THEN max1 = MAX(TB_diff_oc(k, 0, jj))
  ENDFOR  

  PLOT, jd_diff_oc(jj), TB_diff_oc(0, 0, jj), xtitle = 'Date', ytitle = 'TB_meas-TB_calc [K]',$
        title = tit(i_off) +', (90 deg elev,)',$
        yrange = [min1, max1], /NODATA,$
        XTICKFORMAT='LABEL_DATE', XTICKUNITS = ['Time', 'Time'], XTICKS=6 

  kk = [2, 6, 7, 8, 9]
  n_kk = N_ELEMENTS(kk)
  char = STRING(freq_oc, format = '(f5.2)')

  FOR i = 0, n_kk-1 DO BEGIN
   OPLOT, jd_diff_oc(jj), TB_diff_oc(kk(i), 0, jj), color = i+1, psym = sym(1), thick = 2
   XYOUTS, 0.15, 0.9-0.03*FLOAT(i), char(kk(i)) + ' GHz', /normal, color = i+1
  ENDFOR

  OPLOT, [jd_diff_oc(0), jd_diff_oc(N_ELEMENTS(jd_diff_oc)-1)], [0., 0.], linestyle = 2 

;***IWV time series at 90deg elevation
  dummyy = LABEL_DATE(DATE_FORMAT=['%D-%M','%Y'])

  min1 = 0.
  max1 = MAX(iwv_rad_oc(jj)) + MAX(iwv_rad_oc(jj))/10.
 
  PLOT, jd_diff_oc(jj), iwv_rad_oc(jj), xtitle = 'Date', ytitle = 'IWV kg m!e-2!N',$
        title = tit1(i_off),$
        yrange = [min1, max1], /NODATA,$
        XTICKFORMAT='LABEL_DATE', XTICKUNITS = ['Time', 'Time'], XTICKS=6 

  char = ['IWV MWR', 'IWV MODEL', 'IWV RS', 'IWV GPS'] 
  n_char = N_ELEMENTS(char) 

  FOR i = 0, n_char-1 DO BEGIN
   IF i EQ 0 THEN OPLOT, jd_diff_oc(jj), iwv_rad_oc(jj), color = i+1, thick = 2, psym = sym(1)
   IF i EQ 1 THEN OPLOT, jd_diff_oc(jj), iwv_mod_oc(jj), color = i+1, thick = 2, psym = sym(1)
   IF i EQ 2 THEN OPLOT, jd_diff_oc(jj), iwv_rs_oc(jj), color = i+1, thick = 2, psym = sym(1)
   IF i EQ 3 THEN OPLOT, jd_diff_oc(jj), iwv_gps_oc(jj), color = i+1, thick = 2, psym = sym(1)
   XYOUTS, 0.15, 0.15+0.03*FLOAT(i), char(i), /normal, color = i+1
  ENDFOR

;***Scatterplots of different IWV sources

  !p.multi = [0, 2, 2]
  !p.charsize = 0.6
 
  i_jj = WHERE(iwv_rad_oc(jj) GT 0. AND iwv_mod_oc(jj) GT 0., n1)
  IF i_jj(0) NE -1 THEN BEGIN
   BIAS, iwv_rad_oc(i_jj), iwv_mod_oc(i_jj), b1
   r1 = STDDEV(iwv_rad_oc(i_jj)-iwv_mod_oc(i_jj))
   PLOT, iwv_rad_oc(i_jj), iwv_mod_oc(i_jj), psym = 3, xtitle = 'IWV MWR [kgm!e-2!N]', ytitle = 'IWV COSMO [kgm!e-2!N]',$
         xrange = [0., 35.], yrange = [0., 35.], xstyle = 1, ystyle = 1
   OPLOT, [0., 35.], [0., 35.], linestyle = 0
   XYOUTS, 2., 32., 'N = ' + STRING(n1, format = '(i3)') 
   XYOUTS, 2., 30., 'BIAS = ' + STRING(b1, format = '(f6.2)')
   XYOUTS, 2., 28., 'RMS = ' + STRING(r1, format = '(f6.2)')
  ENDIF

  i_jj = WHERE(iwv_rad_oc(jj) GT 0. AND iwv_rs_oc(jj) GT 0., n2)
  IF i_jj(0) NE -1 THEN BEGIN
   BIAS, iwv_rad_oc(i_jj), iwv_rs_oc(i_jj), b2
   r2 = STDDEV(iwv_rad_oc(i_jj)-iwv_rs_oc(i_jj))
   PLOT, iwv_rad_oc(i_jj), iwv_rs_oc(i_jj), psym = 3, xtitle = 'IWV MWR [kgm!e-2!N]', ytitle = 'IWV RS Essen [kgm!e-2!N]',$
         xrange = [0., 35.], yrange = [0., 35.], xstyle = 1, ystyle = 1
   OPLOT, [0., 35.], [0., 35.], linestyle = 0
   XYOUTS, 2., 32., 'N = ' + STRING(n2, format = '(i3)') 
   XYOUTS, 2., 30., 'BIAS = ' + STRING(b2, format = '(f6.2)')
   XYOUTS, 2., 28., 'RMS = ' + STRING(r2, format = '(f6.2)')
  ENDIF

  i_jj = WHERE(iwv_rad_oc(jj) GT 0. AND iwv_gps_oc(jj) GT 0., n3)
  IF i_jj(0) NE -1 THEN BEGIN
   BIAS, iwv_rad_oc(i_jj), iwv_gps_oc(i_jj), b3
   r3 = STDDEV(iwv_rad_oc(i_jj)-iwv_gps_oc(i_jj))
   PLOT, iwv_rad_oc(i_jj), iwv_gps_oc(i_jj), psym = 3, xtitle = 'IWV MWR [kgm!e-2!N]', ytitle = 'IWV GPS [kgm!e-2!N]',$
         xrange = [0., 35.], yrange = [0., 35.], xstyle = 1, ystyle = 1
   OPLOT, [0., 35.], [0., 35.], linestyle = 0
   XYOUTS, 2., 32., 'N = ' + STRING(n3, format = '(i3)') 
   XYOUTS, 2., 30., 'BIAS = ' + STRING(b3, format = '(f6.2)')
   XYOUTS, 2., 28., 'RMS = ' + STRING(r3, format = '(f6.2)')
  ENDIF

  i_jj = WHERE(iwv_rs_oc(jj) GT 0. AND iwv_mod_oc(jj) GT 0., n4)
  IF i_jj(0) NE -1 THEN BEGIN
   BIAS, iwv_rs_oc(i_jj), iwv_mod_oc(i_jj), b4
   r4 = STDDEV(iwv_rs_oc(i_jj)-iwv_mod_oc(i_jj))
   PLOT, iwv_rs_oc(i_jj), iwv_mod_oc(i_jj), psym = 3, xtitle = 'IWV RS Essen [kgm!e-2!N]', ytitle = 'IWV COSMO [kgm!e-2!N]',$
         xrange = [0., 35.], yrange = [0., 35.], xstyle = 1, ystyle = 1
   OPLOT, [0., 35.], [0., 35.], linestyle = 0
   XYOUTS, 2., 32., 'N = '+ STRING(n4, format = '(i3)') 
   XYOUTS, 2., 30., 'BIAS = ' + STRING(b4, format = '(f6.2)')
   XYOUTS, 2., 28., 'RMS = ' + STRING(r4, format = '(f6.2)')
  ENDIF

  DEVICE, /CLOSE

 ENDIF

SKIP_CAL:

ENDFOR; loop over calibration periods



END
