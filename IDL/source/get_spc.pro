;+
;************
PRO GET_SPC,$
;************
;INPUT:
filename,$                  ;*.spc file generated with a software version 7.4 or higher
;OUPUT:
n,$                         ;number of measurements
time,$                      ;seconds since 1.1.1970
rain,$                      ;rain flag (0/1)
f,$                         ;array of frequencies
tb,$                        ;TB arry (f x time)
az,$                        ;array of azimuth angles
el,$                        ;array of elevation angles  
;KEYWORDS
verbose = verbose,$
angle_calc = angle_calc     ;default is 0, 1 means calculate angles in old manner (see below), only necessary for SPC files with file codes 666666
; Abstract:
; * reads binary *.spc files from software versions 7.4 and higher
; Authors:
; U. Loehnert
; Date:
; 2018-08-23
; Dependencies:
; -
; Changes:
; 2014-10-24 (UL)
; - bux fix to to read azimuth and elevation in a consistent manner, not sure when
;   RPG changed there angle coding, so option is built in to switch back to the old 
;   type of angle calculations.
; 2018-02-26 (UL)
; - adjusted to read new angle format introduced with file code = 666000
; 2018-08-23 (UL)
; - adapted the get_brt.pro file to read the RPG spc files
; changed program ...
;-

ON_ERROR, 2

OPENR, unit, filename, /GET_LUN, ERROR=err

IF err NE 0 THEN BEGIN
  n = 0
  GOTO, END_SPC
ENDIF
ON_ERROR, 0

code = 1l
READU, unit, code

IF verbose THEN print, 'File code: ', code

IF code NE 666667 AND code NE 667000 THEN BEGIN 
 print, 'Error in SPC file code'    ; changed in SW version 7.4
 GOTO, END_SPC 
ENDIF

IF N_ELEMENTS(angle_calc) EQ 0 THEN angle_calc = 0

n = 0l
READU, unit, n
IF verbose THEN print, 'Number of samples=', n
IF n LT 1 THEN GOTO, SPC_END

; Time Reference - 1 : UTC , 0 : Local Time

time_ref = 0l
READU, unit, time_ref
IF verbose THEN BEGIN
CASE time_ref of
  1: print, 'UTC'
  0: print, 'Local Time'
  ELSE: print, 'No valid time'
  ENDCASE
ENDIF

n_f = 0l
READU, unit, n_f
IF verbose THEN print, 'Number of frequencies=', n_f

f = FLTARR(n_f)
READU, unit, f
IF verbose THEN print,'Frequencies= ', f

; array of xmin , xmax instead of a total value for min and max - SW V7.4
xmin = MAKE_ARRAY( n_f, /FLOAT, VALUE = -9e+33)
xmax = MAKE_ARRAY( n_f, /FLOAT, VALUE = -9e+33)
READU, unit, xmin
READU, unit, xmax

time = LONARR(n)
rain = BYTARR(n)
tb = FLTARR(n_f, n)
el  = FLTARR(n)
az  = FLTARR(n)

x1 = 0l
x2 = 0b
help = FLTARR(n_f)
IF code EQ 666667 THEN x = 0.
IF code EQ 667000 THEN x = 0l

FOR i = 0l, LONG(n-1l) DO BEGIN

  READU, unit, x1, x2, help, x
  time(i) = x1
  rain(i) = x2
  tb(*, i) = help

  IF code EQ 666667 THEN BEGIN 

;*the following lines (from ## to ##) were created on 24.10.14 to correct the angle calculation from RPG .SPC files
;##
   IF angle_calc EQ 0 THEN BEGIN
    els = x - 100.*FIX(x/100.)
    azs = (x-els)/1000.
  
    IF azs LE 360. THEN BEGIN ; first quadrant
     el(i) = els 
     az(i) = azs
    ENDIF ELSE IF azs GT 360. AND azs LT 1000. THEN BEGIN
     print, 'Inconsistency in angle calculation!'
     print, 'ABORT in get_spc.pro'
     stop
    ENDIF ELSE IF azs GT 1000. THEN BEGIN
     az(i) = azs - 1000.
     el(i) = 100. + els
    ENDIF  
   ENDIF
;##
;*the following lines (from ### to ###) state the angle calculation until 24.10.14
;###
   IF angle_calc EQ 1 THEN BEGIN
    sign = 1
    IF x LT 0. THEN sign = -1.
;    az(i)   = sign*FLOOR(x/100.)/10.
    az(i)   = sign*FIX(x/100.)/10.
    el(i)   = x - (sign*az(i)*1000.)
   ENDIF 
;###
  ENDIF

  IF code EQ 667000 THEN BEGIN 
   c = STRCOMPRESS(x, /REMOVE_ALL)
   c_az = STRMID(c, strlen(c)-5, strlen(c)) 
   c_el = STRMID(c, 0, strlen(c)-5)  
   IF verbose THEN print, x, ' ', c_az, ' ', c_el
   az(i) = FLOAT(c_az)/100.
   el(i) = FLOAT(c_el)/100.
  ENDIF
  IF verbose THEN print, i, time(i), rain(i), tb(1,i), x, az(i), el(i)
ENDFOR

SPC_END:

CLOSE,unit
FREE_LUN,unit

END_SPC:


END
