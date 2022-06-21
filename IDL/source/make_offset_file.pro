;+
;*******************
PRO MAKE_OFFSET_FILE
;*******************
;KEYWORDS
; Abstract:
;* read TB offset files created by
;* tb_offset_correction.pro
;* and store them into the format need by pl_mk.pro
; Author:
; U. Loehnert
; Date:
; 2012-09-11
; Dependencies:
; -
; Changes:
; XXXX-XX-XX: ???
;-

@par_tb_offset_correction

;***PL_MK compatible file
file = 'tb_offset_'+station+'.sav'

;***plot all results
ps_file = 'tb_offset_'+station+'.ps'
tit = 'TB offsets '+station

files = FILE_SEARCH('offsets/OC_'+station+'*.sav')
nfiles = N_ELEMENTS(files)
date_offset = REPLICATE(-99d, nfiles)

IF files(0) NE '' THEN BEGIN

 n_diff_s = 0
 FOR i = 0, nfiles-1 DO BEGIN

  RESTORE, files(i)
  IF i EQ 0 THEN BEGIN 
   freq_offset = freq_oc
   nf = N_ELEMENTS(freq_offset)
   ang_offset = angles_oc
   na = N_ELEMENTS(ang_offset)
   tb_offset = REPLICATE(-99d, nfiles, nf, na)
   tb_diff_a = REPLICATE(-99d, nf, na, 10000) 
   jd_diff_a = REPLICATE(-99d, 10000)
  ENDIF

  date_offset(i) = jd_start_date_oc
  tb_offset(i, *, *) = TB_diff_m_oc   
  n_diff = N_ELEMENTS(TB_diff_oc(0, 0, *))
  FOR i_f = 0, nf-1 DO BEGIN
   FOR i_a = 0, na-1 DO BEGIN
    tb_diff_a(i_f, i_a, n_diff_s:n_diff_s+n_diff-1) = TB_diff_oc(i_f, i_a, *)  
   ENDFOR
  ENDFOR 
  jd_diff_a(n_diff_s:n_diff_s+n_diff-1) = jd_diff_oc
  n_diff_s = n_diff_s + n_diff 

 ENDFOR 
  
ENDIF

jd_diff_a = jd_diff_a(0:n_diff_s-1)
tb_diff_a = tb_diff_a(*, *, 0:n_diff_s-1)

SAVE, filename = file, date_offset, freq_offset, ang_offset, tb_offset

;***plot TB offset time series at 90deg elevation

SET_PLOT, 'ps'
LOADCOL, 'col1'
DEVICE, /color, file=ps_file

!p.multi = [0, 1, 1]

jj = WHERE(tb_diff_a(0, 0, *) NE -99d, n_jj)

IF jj(0) NE -1 THEN BEGIN
 CALDAT, jd_diff_a(jj), mo_jj, da_jj, ye_jj, ho_jj, mi_jj, se_jj

 FOR i = 0, n_jj-1 DO BEGIN

  print, 'time (yy/mm/dd/hh): ', ye_jj(i), mo_jj(i), da_jj(i), ho_jj(i)  
  print, 'tb_diff (23.8, 31.4, 51.3 GHz): ', tb_diff_a(2, 0, jj(i)), tb_diff_a(6, 0, jj(i)), tb_diff_a(7, 0, jj(i))

 ENDFOR
ENDIF

dummyy = LABEL_DATE(DATE_FORMAT=['%D-%M','%Y'])

min1 = 99.
max1 = -99.

FOR k = 0, nf-1 DO BEGIN
 IF MIN(tb_diff_a(k, 0, jj)) LT min1 THEN min1 = MIN(tb_diff_a(k, 0, jj))
 IF MAX(tb_diff_a(k, 0, jj)) GT max1 THEN max1 = MAX(tb_diff_a(k, 0, jj))
ENDFOR

cc = 'lwp_std_thres: '+STRING(lwp_std_thres, format = '(f3.1)')

PLOT, jd_diff_a(jj), tb_diff_a(0, 0, jj), xtitle = 'Date', ytitle = 'TB_meas-TB_calc [K]',$
      title = tit +', (90 deg elev., '+cc+')', yrange = [min1, max1], /NODATA,$
      XTICKFORMAT='LABEL_DATE', XTICKUNITS = ['Time', 'Time'], XTICKS=6

kk = [2, 6, 7, 8, 9]
n_kk = N_ELEMENTS(kk)   
char = STRING(freq_oc, format = '(f5.2)')

FOR i = 0, n_kk-1 DO BEGIN
 OPLOT, jd_diff_a(jj), tb_diff_a(kk(i), 0, jj), color = i+1, psym = sym(1), thick = 2
 XYOUTS, 0.15, 0.9-0.03*FLOAT(i), char(kk(i)) + ' GHz', /normal, color = i+1
ENDFOR

OPLOT, [jd_diff_a(0), jd_diff_a(N_ELEMENTS(jd_diff_a)-1)], [0., 0.], linestyle = 2

n_offset = N_ELEMENTS(date_offset)

FOR i = 0, n_offset-1 DO BEGIN
 OPLOT, [date_offset(i), date_offset(i)], [min1, max1], thick = 5, color = 7
ENDFOR

DEVICE, /close

stop

END
