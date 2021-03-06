;+
;***************************
PRO READ_LEVEL1_MWRBL_NC,$
;***************************
;INPUT:
filename,$                      ;name of netcdf file to be read
;OUTPUT:
comment,$                       ;any specific comments
latitude,longitude,altitude,$   ;three FLOAT
angles,$                        ;array of elevation angles
freq,$                          ;array of HATPRO frequencies
time,$                          ;array of decimal hours
tb,$                            ;3D TB array (freq x angles x time)
flag,$                          ;data flagging
temp,$                          ;surface temperature in K
pres,$                          ;surface pressure in hPa
relh                            ;surface relative humidity in %
; $Id: read_hatpro_level0c_nc.pro,v 1.2 2010/12/27 17:35:32 loehnert Exp $
; Abstract:
; * read HATPRO level1 BL netcdf files
; Authors:
; U. Loehnert
; Date:
; 2008-09-19
; Dependencies:
; -
; Changes:
; 2019-08-16
; updated to HD(CP)2 output format
;-

infid  = NCDF_OPEN(filename, /NOWRITE)
tid    = NCDF_VARID(infid,'time')
fid    = NCDF_VARID(infid,'freq_sb')
tbid   = NCDF_VARID(infid,'tb')
latid  = NCDF_VARID(infid, 'lat')
longid = NCDF_VARID(infid, 'lon')
altid  = NCDF_VARID(infid, 'zsl')
eleid  = NCDF_VARID(infid, 'ele')
aziid  = NCDF_VARID(infid, 'azi') 
rid    = NCDF_VARID(infid, 'flag')
teid   = NCDF_VARID(infid, 'ta')
pid    = NCDF_VARID(infid, 'pa')
rhid   = NCDF_VARID(infid, 'hur')
  
NCDF_VARGET, infid, latid,  latitude
NCDF_VARGET, infid, longid, longitude
NCDF_VARGET, infid, altid,  altitude
NCDF_VARGET, infid, anglid, angles
NCDF_VARGET, infid, tid,    time
NCDF_VARGET, infid, fid,    freq
NCDF_VARGET, infid, tbid,   tb
NCDF_VARGET, infid, rid,    flag
NCDF_VARGET, infid, teid,   temp
NCDF_VARGET, infid, pid,    pres
NCDF_VARGET, infid, rhid,   relh

; get 1st attribute name = NCDF_ATTNAME(infid,/GLOBAL,1)
;NCDF_ATTGET, infid, /GLOBAL, 'comments', com
;comment = STRING(BYTE(com))

NCDF_CLOSE, infid

n = N_ELEMENTS(time)

END
