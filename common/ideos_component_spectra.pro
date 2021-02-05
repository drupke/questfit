; docformat = 'rst'
;
;+
;
; Extracts component spectra from each SAV file output by QUESTFIT.
;
; :Categories:
;    IDEOS
;
; :Returns:
;    XDR files with component spectra.
;
; :Params:
;    idlist: in, required, type=string
;      File with list of fits to include, given by IDEOS ID. One per line.
;    outdir: in, required, type=string
;      Location of output XDR files.
;
; :Keywords:
;    savdir: in, required, type=string
;      Directory where .sav files are located. Default is current directory.
; 
; :Author:
;    Vincent Viola::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2014--2015, VV, created
;      2015nov20, DSNR, copied from fit_feature_spectra.pro
;      2015dec08, DSNR, added ability to remove ice absorption
;      2016aug25, DSNR, now autodetect # of rows in input list
;      2016aug31, DSNR, changed call to INTERPOLATE to SAP_REBIN to match
;                       the procedure in QUESTFIT for interpolating
;                       extinction and ice curves
;      2016oct25, DSNR, added logic for case of 4 BB components
;    
; :Copyright:
;    Copyright (C) 2015--2016 Vincent Viola, David S. N. Rupke
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
pro ideos_component_spectra, idlist, outdir, savdir=savdir, icefile=icefile

   if not keyword_set(savdir) then savdir='./'

;  if requested, get ice absorption
   readcol,icefile,icewl,icetau,/silent,skipline=1

;  Galaxy fit results table
   nrows = file_lines(idlist)
   rows=[2,nrows]
   idarr = read_csvcol(idlist,1,type='str',sep=',',$
                       rows=rows)
   clustarr = read_csvcol(idlist,2,type='str',sep=',',$
                          rows=rows)
   num=n_elements(idarr)
   for ind=0,num-1 do begin

      label=strcompress(idarr[ind]+'_'+clustarr[ind],/remove_all)
      restore,savdir+label+'.sav'

;     Get number of blackbody components
      sizebbmem = size(bbmem)
      nbb = sizebbmem[1]

      nice = n_elements(icetau)
      bb1_noice = bbmem[0,*]
      bb2_noice = bbmem[1,*]
      bb3_noice = bbmem[2,*]
      if keyword_set(icefile) then begin
         status_rb=sap_rebin(icewl,icetau,dblarr(nice),$
                             icewl_int,icetau_int,iceerr_int,$
                             ref=sourcewave[indsource])
         bb1_noice *= exp(result[1]*icetau_int)
         bb2_noice *= exp(result[5]*icetau_int)
         bb3_noice *= exp(result[9]*icetau_int)
      endif
      
      features=create_struct('Wavelength',sourcewave[indsource],$
                             'Blackbody1',bbmem[0,*],$
                             'Blackbody2',bbmem[1,*],$
                             'Blackbody3',bbmem[2,*],$
                             'Blackbody1_noice',bb1_noice,$
                             'Blackbody2_noice',bb2_noice,$
                             'Blackbody3_noice',bb3_noice,$
                             'PAH1',tempmem[0,*],$
                             'PAH2',tempmem[1,*])
      if nbb eq 4 then begin
         bb4_noice = bbmem[3,*]
         if keyword_set(icefile) then $
            bb4_noice *= exp(result[13]*icetau_int)
         features = create_struct(features,$
                                  'Blackbody4',bbmem[3,*],$
                                  'Blackbody4_noice',bb4_noice)
      endif

      save,features,filename=outdir+label+'.xdr',/xdr   ;save file

   endfor

end
