; docformat = 'rst'
;
;+
;
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
;      2015dec01, DSNR, created
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
pro ideos_plot

   rootdir='/Volumes/fingolfin/spitzer/ideos/'

   readcol,rootdir+'docs/ideos_pahweq.txt',id,pah6,pah7,pah11,$
           /silent,skip=7,format='(A,X,X,X,X,D,D,D)'
   readcol,rootdir+'docs/ideos_tau.txt',id,siltau,$
           /silent,skip=5,format='(A,X,D)'

   cgps_open,rootdir+'plots/ideos_pah6_v_siltau.eps',charsize=1,/encap,$
             /inches,xs=7.5,ys=7.5,/qui
   cgplot,[0],/xsty,/ysty,yran=[-0.6,1.3],xran=[-4,1]
   cgoplot,pah6,alog10(siltau),psym=16
   cgps_close

   xtit='log[ EW(7.7$\mu$m PAH) / $\mu$m]'
   ytit='$\tau$$\down9.7$$\upeffective$'
   cgps_open,rootdir+'plots/ideos_pah7_v_siltau.eps',charsize=1.1,/encap,$
             /inches,xs=7.5,ys=7.5,/qui
   cgplot,[0],/xsty,/ysty,yran=[-0.99,13.99],xran=[-1.99,1.5],xtit=xtit,$
          ytit=ytit
   cgoplot,pah7,siltau,psym=16
   cgoplot,[-2,0.8],[12,-1],linesty=2
   cgoplot,[-0.3,1.6],[4,8.5],linesty=2
   cgps_close
   
end