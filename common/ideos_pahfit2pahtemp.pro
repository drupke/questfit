; docformat = 'rst'
;
;+
;
; Convert output from PAHFIT into a pure PAH spectrum.
; Compute and print to the screen a few PAH line ratios.
; Print the spectrum to a file, if desired.
; Print the PAH parameters to a file, if desired.
;
; :Categories:
;    IDEOS
;
; :Returns:
;
; :Params:
;    datafile: in, required, type=string
;      File of best-fit parameters, output by IDEOS_BESTFITPARTAB
;
; :Keywords:
;    outdir: in, required, type=string, default=./
;      Directory for output files
;
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2007sep27  DSNR  created
;      2015dec11, DSNR, documented
;
; :Copyright:
;    Copyright (C) 2015 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
function ideos_pahfit2pahtemp,fit,file=file,parfile=parfile

  li=1d
  lf=40d
  ldelt=0.01d

  ; populate wavelength array
  l=dblarr((lf-li)/ldelt)
  f=l
  l[0]=li
  for i=1,n_elements(l)-1 do l[i]=l[i-1]+ldelt

  if (keyword_set(parfile)) then openw,unit0,parfile,/get_lun

  f6=0d
  f7=0d
  f11=0d
  f17=0d
  for i=0,n_elements(fit.dust_features)-1 do begin
     for j=0,n_elements(l)-1 do $
        f[j]=f[j]+fit.dust_features[i].central_inten*$
              fit.dust_features[i].fwhm^2d / $
              ((l[j]/fit.dust_features[i].wavelength - $
                fit.dust_features[i].wavelength/l[j])^2d + $
               fit.dust_features[i].fwhm^2d)
     if (i eq 2) then f6=f6+fit.dust_features[i].int_strength
     if (i eq 4 OR i eq 5 OR i eq 6) then f7=f7+fit.dust_features[i].int_strength
     if (i eq 10 OR i eq 11) then f11=f11+fit.dust_features[i].int_strength
     if (i eq 19 OR i eq 20 OR i eq 21 OR i eq 22) then f17=f17+fit.dust_features[i].int_strength
     if (keyword_set(parfile)) then $
        printf,unit0,fit.dust_features[i].wavelength,$
               fit.dust_features[i].fwhm,$
               fit.dust_features[i].central_inten,$
               fit.dust_features[i].int_strength
  endfor

; normalize to unity
  f=f/max(f)

; print 11/17 ratio to screen
  print,"The 11/17 ratio is ",f11/f17
  print,"The  6/ 7 ratio is ",f6/f7
  print,"The  7/11 ratio is ",f7/f11

; print PAH template to file
  if (keyword_set(file)) then begin
     openw,unit1,file,/get_lun
     for i=1,n_elements(l)-1 do begin
        printf,unit1,l[i],f[i]
     endfor
     close,unit1
   endif

   if (keyword_set(parfile)) then close,unit0

; return the pure PAH spectrum as an array of arrays
   return,[[l],[f]]

end
