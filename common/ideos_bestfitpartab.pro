; docformat = 'rst'
;
;+
;
; Creates a table of best-fit parameters for all fits. Power-law parameters 
; are ignored. For order of parameters, see *.dat files output by QUESTFIT.
; 
;
; :Categories:
;    IDEOS
;
; :Returns:
;    Data file with PAH fluxes.
;
; :Params:
;    idlist: in, required, type=string
;      File with list of fits to include, given by IDEOS ID. One per line.
;      First column is ID, second is "cluster" number, third column is 
;      best fit (1, 2, or 3). 
;    outfile: in, required, type=string
;      File with output table.
;
; :Keywords:
;    savdir: in, optional, type=string, default='./'
;      Directory where .sav files are located.
;    skipline: in, optional, type=integer, default=0
;      Number of lines to skip in input file.
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
;      2013--2015, VV, created
;      2015nov17, DSNR, copied from readFitResults.pro; added docs
;      2015dec01, DSNR, added printing of best fit and sorting by ID;
;                       uses READCOL instead of READ_CSVCOL
;      2016oct25, DSNR, added logic for case of 4 BB components
;      2021dec08, DSNR, fixed treatment of IDs with leading zeros
;    
; :Copyright:
;    Copyright (C) 2015--2021 Vincent Viola, David S. N. Rupke
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
pro ideos_bestfitpartab,idlist,outfile,savdir=savdir,skipline=skipline

   if not keyword_set(savdir) then savdir='./'

;  Galaxy fit results table
   if not keyword_set(skipline) then skipline=0
   readcol,idlist,idarr,clustarr,bestfit,format='(A,I,I)',skipline=skipline,$
           delim=',',/silent,nlines=num

;  Sort by IDEOS ID
   isort_idarr = sort(long(idarr))
   s_idarr = idarr[isort_idarr]
   s_clustarr = clustarr[isort_idarr]
   s_bestfit = bestfit[isort_idarr]

   openw,lun,outfile,/get_lun
;  These indices only apply for case of 4 components (including absorption);
;  can be either BB or PL
   igoodpah1=[dindgen(2)+16]
   igoodpah2=[dindgen(2)+19]
   for ind=0,num-1 do begin
      label=string(s_idarr[ind],'_',s_clustarr[ind],format='(A0,A0,I0)')
      printf,lun,label,s_bestfit[ind],format='(A-15,I3,$)'
      restore,savdir+label+'.sav'
;     This bit of logic only works for numberofBB = 3 or 4
      igoodbb=[dindgen(numberofBB*4)]
      printf,lun,result[igoodbb],format='('+$
             string(numberofBB*4,format='(I0)')+'E12.4,$)'
;     Set Temp to something other than 0; otherwise IDEOS_PAH doesn't process
;     (absent) BB#4 correctly
      if numberofBB eq 3 then $
         printf,lun,0d,0d,9999d,0d,format='(4E12.4,$)'
      printf,lun,result[igoodpah1],format='(2E12.4,$)'
      printf,lun,tempmemnorm[0],format='(E12.4,$)'
      printf,lun,result[igoodpah2],format='(2E12.4,$)'
      printf,lun,tempmemnorm[1],format='(E12.4,$)'
      printf,lun,minsourcewave,maxsourcewave,format='(2D8.4)'
   endfor

   free_lun,lun

end
