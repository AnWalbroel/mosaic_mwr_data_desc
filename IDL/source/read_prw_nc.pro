;+
;********************
PRO READ_PRW_NC,$
;********************
;INPUT:
filename,$                      ;name of netcdf file to be read
;OUTPUT:
time,$                          ;array of decimal hours
ele,$                           ;array of elevation angles 
azi,$                           ;array of azimuth angles
prw,$                           ;array of uncorrected pwr
prw_oz,$                        ;array of uncorrected off-zenith integrated water vapor 
flag                            ;data flagging

; $Id: $
; Abstract:
; * read HATPRO level2a netcdf files
; Authors:
; U. Loehnert
; Date:
; 2012-09-07
; Dependencies:
; -
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

infid = NCDF_OPEN(filename(0), /NOWRITE)
tid  = NCDF_VARID(infid,'time')
latid = NCDF_VARID(infid, 'lat')
lonid = NCDF_VARID(infid, 'lon')
altid = NCDF_VARID(infid, 'zsl')
elid  = NCDF_VARID(infid, 'ele')
azid  = NCDF_VARID(infid, 'azi')
prwid  = NCDF_VARID(infid, 'prw')
prwozid = NCDF_VARID(infid, 'prw_off_zenith')
rid = NCDF_VARID(infid, 'flag')
  
NCDF_VARGET, infid, tid, time
NCDF_VARGET, infid, latid, lat
NCDF_VARGET, infid, lonid, lon
NCDF_VARGET, infid, elid, ele
NCDF_VARGET, infid, azid, azi
NCDF_VARGET, infid, prwid, prw
NCDF_VARGET, infid, prwozid, prw_oz
NCDF_VARGET, infid, rid, flag

END
