; docformat = 'rst'
;
;+
;
; :Categories:
;    IDEOS
;
; :Returns:
;
; :Params:
;
; :Keywords:
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
;      2015dec15, DSNR, created
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
pro ideos_testpahtemp

   rootdir='/Volumes/fingolfin/spitzer/ideos/'

   restore,rootdir+'docs/datafiles/smith_nftemp3.xdr'
   lx = double(memstor.aar.data.wave)
   fx = double(memstor.aar.data.flux)
   readcol300,rootdir+'docs/datafiles/smith_nftemp3.pah.ext.dat',$
              ld,fd,/silent,format='(D,D)'

   npts = n_elements(ld)
   
   readcol300,rootdir+'docs/datafiles/smith_nftemp3.pahpars.ext.dat',$
              pahwave,fwhm,peakint,totint,/silent,format='(D,D,D,D)'
   f = dblarr(npts)
   f6=0d
   f7=0d
   f11=0d
   f17=0d
   for i=0,n_elements(pahwave)-1 do begin
      f += peakint[i]*fwhm[i]^2d / $
           ((ld/pahwave[i] - pahwave[i]/ld)^2d + fwhm[i]^2d)
      if (i eq 2) then f6+=totint[i]
      if (i eq 4 OR i eq 5 OR i eq 6) then f7+=totint[i]
      if (i eq 10 OR i eq 11) then f11+=totint[i]
      if (i eq 19 OR i eq 20 OR i eq 21 OR i eq 22) then f17+=totint[i]
   endfor

;   f/=max(f)

   ; print 11/17 ratio to screen
   print,"The 11/17 ratio is ",f11/f17
   print,"The  6/ 7 ratio is ",f6/f7
   print,"The  7/11 ratio is ",f7/f11

   cgps_open,rootdir+'docs/datafiles/smith_nftemp3_test.eps',charsize=1,/encap,$
             /inches,xs=7.5d,ys=7.5d,/qui
   cgplot,ld,f/fd,xran=[5,32.5],yran=[0.9,1.1],/xsty,/ysty
;   cgoplot,lx,fd/fx,color='Red'
   cgps_close

end
