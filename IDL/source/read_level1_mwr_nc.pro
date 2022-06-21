;+
;**************************
PRO READ_LEVEL1_MWR_NC,$
;**************************
;INPUT:
filename,$                      ;name of netcdf file to be read
;OUTPUT:
comment,$                       ;any specific comments
time,$                          ;array of decimal hours
freq_hat,$			;Hatpro frequencies [GHz]
tb_hat,$			;Hatpro TBs [K]
wavel_ir,$			;Center wavelengths of IRT Brightness Temperatures (microns)
tb_ir,$				;broadband IR brightness temperature(s), azimuth angle same as MWR
ele_ir,$			;elevation angle IR radiometer (degree)
flag,$                          ;data flagging
temp,$                          ;surface temperature in K
pres,$                          ;surface pressure in hPa
relh,$                          ;surface relative humidity in %
az,$                            ;array of azimuth angles (if available)
el, $                           ;array of elevation angles
; geographical coordinates
latitude = latitude, $          ; latitude of radiometer position (float)
longitude = longitude, $        ; longitude of radiometer position (float)
altitude = altitude,$           ; altitude of radiometer position
; infos from global attributes
;Date_of_file_generation = Date_of_file_generation, $ ; ...
;Conventions = Conventions, $    ; version of the netcdf fileformat
;location = location, $          ; ...
;system = system, $              ; type and name of the radiometer
;day = day,$                     ; day of dataset
;month = month,$                 ; month of dataset
;year = year, $                  ; year of dataset
;institution = institution, $    ; ...
; additional info
N_data = N_data                 ; number of elements in time
;
; $Id:  $
; Abstract:
; * read level1 mwr  netcdf output files
; Authors:
; U. Loehnert
; Date:
; 2019-02-14
; Dependencies:
;   none
; Changes:
;-

infid  = NCDF_OPEN(filename, /NOWRITE)
tid    = NCDF_VARID(infid,'time')
latid  = NCDF_VARID(infid, 'lat')
longid = NCDF_VARID(infid, 'lon')
altid  = NCDF_VARID(infid, 'zsl')
tid    = NCDF_VARID(infid, 'time')
fid    = NCDF_VARID(infid, 'freq_sb')
tbid   = NCDF_VARID(infid, 'tb')
wlirid = NCDF_VARID(infid, 'wl_irp')
tbirid = NCDF_VARID(infid, 'tb_irp')
elirid = NCDF_VARID(infid, 'ele_irp')
rid    = NCDF_VARID(infid, 'flag')
teid   = NCDF_VARID(infid, 'ta')
pid    = NCDF_VARID(infid, 'pa')
rhid   = NCDF_VARID(infid, 'hur')
azid   = NCDF_VARID(infid, 'azi')
elid   = NCDF_VARID(infid, 'ele')  

NCDF_VARGET, infid, latid,  latitude
NCDF_VARGET, infid, longid, longitude
NCDF_VARGET, infid, altid,  altitude
NCDF_VARGET, infid, tid,    time
NCDF_VARGET, infid, tLXid,  time_Linux
NCDF_VARGET, infid, fid,    freq_hat
NCDF_VARGET, infid, tbid,   tb_hat
NCDF_VARGET, infid, wlirid, wavel_ir
NCDF_VARGET, infid, tbirid, tb_ir
NCDF_VARGET, infid, elirid, ele_ir
NCDF_VARGET, infid, rid,    flag
NCDF_VARGET, infid, teid,   temp
NCDF_VARGET, infid, pid,    pres
NCDF_VARGET, infid, rhid,   relh
NCDF_VARGET, infid, azid,   az
NCDF_VARGET, infid, elid,   el


; read global attributes 
; they are strings which are read as array of byte by idl
; they must be converted with string(byte(...))

; get 1st attribute name = NCDF_ATTNAME(infid,/GLOBAL,1)
;NCDF_ATTGET,infid,/GLOBAL,'comments',com
;comment = STRING(BYTE(com))

;NCDF_ATTGET,infid,/GLOBAL,'Date_of_file_generation',Date_of_file_generation
;Date_of_file_generation = string(byte(Date_of_file_generation))

,NCDF_ATTGET,infid,/GLOBAL,'system',system
;system = string(byte(system))

;NCDF_ATTGET,infid,/GLOBAL,'day',   day
;day = string(byte(day))

;NCDF_ATTGET,infid,/GLOBAL,'month', month
;month = string(byte(month))

;NCDF_ATTGET,infid,/GLOBAL,'year',  year
;year = string(byte(year))

;NCDF_ATTGET,infid,/GLOBAL,'institution',  institution
;institution = string(byte(institution))

NCDF_CLOSE, infid

N_data = N_ELEMENTS(time)

END
