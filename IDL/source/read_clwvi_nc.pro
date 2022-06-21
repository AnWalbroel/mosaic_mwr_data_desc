;+
;********************
PRO READ_CLWVI_NC,$
;********************
;INPUT:
filename,$                      ;name of netcdf file to be read
;OUTPUT:
time,$                          ;array of decimal hours
ele,$                           ;array of elevation angles 
azi,$                           ;array of azimuth angles
lwp_oc,$                        ;array of liquid water path, with zeroing
lwp_org,$                       ;array of uncorrected liquid water path
lwp_oz,$                        ;array of uncorrected off-zenith liquid water path
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
lwporgid  = NCDF_VARID(infid, 'clwvi')
lwpozid = NCDF_VARID(infid, 'clwvi_off_zenith')
lwpocid = NCDF_VARID(infid, 'clwvi_offset_zeroing')
rid = NCDF_VARID(infid, 'flag')
  
NCDF_VARGET, infid, tid, time
NCDF_VARGET, infid, latid, lat
NCDF_VARGET, infid, lonid, lon
NCDF_VARGET, infid, elid, ele
NCDF_VARGET, infid, azid, azi
NCDF_VARGET, infid, lwporgid, lwp_org
NCDF_VARGET, infid, lwpozid, lwp_oz
NCDF_VARGET, infid, lwpocid, lwp_oc
NCDF_VARGET, infid, rid, flag

END
