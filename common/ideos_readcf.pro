; docformat = 'rst'
;
;+
;
; Reads control file for QUESTFIT.
;
; :Categories:
;    IDEOS
;
; :Returns:
;    A whole bunch of variables ...
;    
; :Params:
;
; :Keywords:
;     pathin: in, optional, type=string, default='./'
; 
; :Author:
;    Mario Schweitzer::
;      MPE
;      Garching, Germany
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
;      2015nov20, DSNR, copied from readcf.pro; changed path specification
;      2016aug24, DSNR, changed NUMLINES call to FILE_LINES; small change to
;                       error message
;    
; :Copyright:
;    Copyright (C) 2015 Mario Schweitzer, Vincent Viola, David S. N. Rupke
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
pro ideos_readcf,controlfile=controlfile,numberofBB,numberofPL,numberoftemplates,$
                 numberofext,sourcefiles,range,templatefiles,templatefilesb,$
                 normtemp,fixfreetemp,extcurvetemp,extvaluetemp,normbb,fixfreebb,$
                 fixfreebbtemp,extcurvebb,extvaluebb,tempbb,normpl,fixfreepl,$
                 fixfreeplindex,extcurvepl,extvaluepl,indexpl,templatewave,$
                 templateflux,sourcewave,sourceflux,sourcestdev,maxwavelengthticks,$
                 sizetemp,sizeext,extwave,extvalue,extvaluefixfreetemp,$
                 extvaluefixfreebb,extvaluefixfreepl,extcurvesynonym,$
                 screenmixedBB,screenmixedpl,screenmixedtemplate,taupeaktemp,$
                 taupeakbb,taupeakpl,taupeaktempfixfree,taupeakbbfixfree,$
                 taupeakplfixfree, numberofabstemp, numberofabsbb, numberofabspl,$
                 absorbtempwave,absorbbbwave,absorbplwave,absorbtemptau,$
                 absorbbbtau,absorbpltau,abstempsize,absbbsize,absplsize,z,$
                 title,ID,cluster,pathin=pathin

if not keyword_set(pathin) then path='./' else path=pathin

readcol,pathin+controlfile,f='A,A,F,F,A,F,F,A,F,F',dtype,dname,norm,fixfree,$
        extcurve,extvalue,num,model,temp,nume,/silent
  
;dtype is data type, dname is the file name, norm is the normalization
;factor, fixfree is the fix or free parameter, extcurve is the synonym
;of the extinction curve, extvalue is average extinction value, temp
;is black body temp, will investigate what num, model and nume are

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NOTE: MEMSTOR is the variable restored after restoring any XDR file
;necessary for this program and that is where most of your important
;data is read from when you see memstor followed by an extention
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 
  numberofline=file_lines(pathin+controlfile)
  numberoftemplates=0
  numberofBB=0
  numberofpl=0
  numberofext=0
  numberofabs=0
  numberofsources=0

 
  for i=0,numberofline-1 do begin                    ;counting how much of each data type

      if dtype(i) eq 'extinction' then numberofext+=1
      if dtype(i) eq 'template' then numberoftemplates+=1
      if dtype(i) eq 'blackbody' then numberofbb+=1
      if dtype(i) eq 'powerlaw' then numberofpl+=1
      if dtype(i) eq 'absorption' then numberofabs+=1
      if dtype(i) eq 'source' then numberofsources+=1

  endfor

  sourcefiles=strarr(numberofsources)

  templatefiles=strarr(numberoftemplates)                     ;creates an array containing the number of template files
  templatefilesb=strarr(numberoftemplates)                    ;and then creates arrays for the normalization, fix free values,
  normtemp=fltarr(numberoftemplates)                          ;etc. based on the number of template files
  fixfreetemp=fltarr(numberoftemplates)
  extcurvetemp=strarr(numberoftemplates)
  extvaluetemp=fltarr(numberoftemplates)
  extvaluefixfreetemp=fltarr(numberoftemplates)
  sizetemp=fltarr(numberoftemplates)
  screenmixedtemplate=strarr(numberoftemplates)
  taupeaktemp=fltarr(numberoftemplates,numberofabs)
  taupeaktempfixfree=fltarr(numberoftemplates,numberofabs)
  numberofabstemp=fltarr(numberoftemplates)
  absorbtempwave=fltarr(numberoftemplates,numberofabs,1.e4)
  absorbtemptau=fltarr(numberoftemplates,numberofabs,1.e4) 
  abstempsize=fltarr(numberoftemplates,numberofabs)
  

  blackbody=strarr(numberofBB)                     ;repeat last comment block for blackbodies
  normbb=fltarr(numberofBB)                        
  fixfreebb=fltarr(numberofBB)
  fixfreebbtemp=fltarr(numberofBB)
  extcurvebb=strarr(numberofBB)
  extvaluebb=fltarr(numberofBB)
  extvaluefixfreebb=fltarr(numberofBB)
  tempbb=fltarr(numberofBB)
  screenmixedBB=strarr(numberofBB)
  taupeakBB=fltarr(numberofBB,numberofabs)
  taupeakBBfixfree=fltarr(numberofBB,numberofabs)
  numberofabsBB=fltarr(numberofBB)
  absorbBBwave=fltarr(numberofBB,numberofabs,1.e4) ; maximum 1.e4 wavelengthticks for absfeatures !
  absorbBBtau=fltarr(numberofBB,numberofabs,1.e4) 
  absbbsize=fltarr(numberofBB,numberofabs) 

  if numberofpl gt 0 then begin 
    powerlaw=strarr(numberofpl)                        ;repeat last comment block for powerlaws
    normpl=fltarr(numberofpl)
    fixfreepl=fltarr(numberofpl)
    fixfreeplindex=fltarr(numberofpl)
    extcurvepl=strarr(numberofpl)
    extvaluepl=fltarr(numberofpl)
    extvaluefixfreepl=fltarr(numberofpl)
    indexpl=fltarr(numberofpl)
    screenmixedpl=strarr(numberofpl)
    screenmixedpl=strarr(numberofpl)
    taupeakpl=fltarr(numberofpl,numberofabs)
    taupeakplfixfree=fltarr(numberofpl,numberofabs)
    numberofabspl=fltarr(numberofpl)
    absorbplwave=fltarr(numberofpl,numberofabs,1.e4)
    absorbpltau=fltarr(numberofpl,numberofabs,1.e4) 
    absplsize=fltarr(numberofpl,numberofabs)
  endif else begin
    numberofabspl=0
  endelse

  extinctionfiles=strarr(numberofext)               ;creates an array for the extinction files and extinction synonyms
  extcurvesynonym=strarr(numberofext)
  extfiles=strarr(numberofext)
  sizeext=fltarr(numberofext)

  absorptionfiles=strarr(numberofabs)              ;same for absorption files
  absorpsynonym=strarr(numberofabs)
  absfiles=strarr(numberofabs)
  sizeabs=fltarr(numberofabs)

  counta=0 ;count of absorption component
  countb=0 ;count of BB
  counte=0 ;count of extinction
  countp=0 ;count of PL
  countt=0 ;count of templates
  counts=0 ;count of sources

  countabsorb=0

  for i=0,numberofline-1 do begin

      if dtype(i) eq 'source' then begin   ;sets variable 'sourcefile' equal to the source filename
          sourcefiles(counts)=dname(i)              ;and the range at the from the lower wavelength limit to the
          range=[norm(i),fixfree(i)]       ;upper wavlength limit, which are in the norm and fixfree columns respectively
          counts=counts+1
      endif
      ;;;;;;;;;;;;;;;;
      if dtype(i) eq 'template' then begin
          templatefiles(countt)=dname(i)      ;sets current template name to a position in the templatefiles array
          templatefilesb(countt)=dname(i)     ;and sets it's corresponding parameterization values to it
          normtemp(countt)=norm(i)
          fixfreetemp(countt)=fixfree(i)
          extcurvetemp(countt)=extcurve(i)
          extvaluetemp(countt)=extvalue(i)
          extvaluefixfreetemp(countt)=num(i)
          screenmixedtemplate(countt)=model(i)

          ii=i+1                                   ;following while loop checks for how many absoprtion files are attributed to the current template
          countabsorb=0                            ;and sets the absorption curve and tau peak to a position in the absorbtempwave array

          while dtype(ii) eq 'absorption' do begin

              taupeaktemp(countt,countabsorb)=norm(ii)
              taupeaktempfixfree(countt,countabsorb)=fixfree(ii)              
              restore,path+dname(ii)
              abstempsize(countt,countabsorb)=n_elements(memstor.aar.data.wave)
              absorbtempwave(countt,countabsorb,0:abstempsize(countt,countabsorb)-1)=memstor.aar.data.wave
              absorbtemptau(countt,countabsorb,0:abstempsize(countt,countabsorb)-1)=memstor.aar.data.flux

              countabsorb=countabsorb+1
              ii=ii+1

          endwhile

          numberofabstemp(countt)=countabsorb
         
       	  countt=countt+1

      endif
      
      if dtype(i) eq 'blackbody' then begin  ;this if statement and following while loop repeats the same
          blackbody(countb)=dname(i)         ;set of procedures for blackbodies as was completed with templates
          normbb(countb)=norm(i)
          fixfreebb(countb)=fixfree(i)
          fixfreebbtemp(countb)=nume(i)
          extcurvebb(countb)=extcurve(i)
          extvaluebb(countb)=extvalue(i)
          extvaluefixfreebb(countb)=num(i)
          tempbb(countb)=temp(i)
          screenmixedBB(countb)=model(i)

          ii=i+1
          countabsorb=0

          while dtype(ii) eq 'absorption' do begin

              taupeakBB(countb,countabsorb)=norm(ii)
              taupeakBBfixfree(countb,countabsorb)=fixfree(ii)              
              restore,path+dname(ii)
              absbbsize(countb,countabsorb)=n_elements(memstor.aar.data.wave)
              absorbBBwave(countb,countabsorb,0:absbbsize(countb,countabsorb)-1)=memstor.aar.data.wave
              absorbBBtau(countb,countabsorb,0:absbbsize(countb,countabsorb)-1)=memstor.aar.data.flux
              countabsorb=countabsorb+1
          
              ii=ii+1

          endwhile

          numberofabsBB(countb)=countabsorb
      
	  countb=countb+1
      endif
      
      if dtype(i) eq 'powerlaw' then begin  ;same for powerlaws
          powerlaw(countp)=dname(i)
          normpl(countp)=norm(i)
          fixfreepl(countp)=fixfree(i)
          fixfreeplindex(countp)=nume(i)
          extcurvepl(countp)=extcurve(i)
          extvaluepl(countp)=extvalue(i)
          extvaluefixfreepl(countp)=num(i)
          indexpl(countp)=temp(i)
          screenmixedpl(countp)=model(i)

          ii=i+1
          countabsorb=0

          while dtype(ii) eq 'absorption' do begin

              taupeakpl(countp,countabsorb)=norm(ii)
              taupeakplfixfree(countp,countabsorb)=fixfree(ii)              
              restore,path+dname(ii)
              absplsize(countp,countabsorb)=n_elements(memstor.aar.data.wave)
              absorbplwave(countp,countabsorb,0:absplsize(countp,countabsorb)-1)=memstor.aar.data.wave
              absorbpltau(countp,countabsorb,0:absplsize(countp,countabsorb)-1)=memstor.aar.data.flux
              countabsorb=countabsorb+1
              ii=ii+1

          endwhile

          numberofabspl(countp)=countabsorb
      
	  countp=countp+1
      endif
      
      if dtype(i) eq 'extinction' then begin  ;checks for extinction files, names them appropriately and places them
          extinctionfiles(counte)=dname(i)         ;in the extinction files array appropriately
          extcurvesynonym(counte)=extcurve(i)
	  counte=counte+1
      endif
      
      if dtype(i) eq 'absorption' then begin  ;repeat for absorption files not already attributed to a template,
          absorptionfiles(counta)=dname(i)         ;blackbody or powerlaw
          absorpsynonym(counta)=extcurve(i)
	  counta=counta+1
      endif
      
  endfor


; read in the templates
;--------------------------------------------------------------------------------------------------------
  maxwavelengthticks=0
  minwavelengthticks=1e10 ; will be reseted in first loop (asumes there is  no template with ticks >=1e10)
  for i=0,numberoftemplates-1 do begin
      
      templatefiles(i)=path+templatefiles(i)            ;reads number of wavelength ticks from all templatefiles
      restore, templatefiles(i)                         ;after restoring them and calling them from the memstor variable
      wavelengthticks=n_elements(memstor.aar.data.wave)
      
      if wavelengthticks gt maxwavelengthticks then begin
          maxwavelengthticks=wavelengthticks
      endif
      
      if wavelengthticks lt minwavelengthticks then begin  ;replaces maxwavelengthticks or minwavelengthticks should a template file
          minwavelengthticks=wavelengthticks               ;with a larger or lesser amount of max or min ticks be read in
      endif
  endfor

  templatewave=fltarr(50,maxwavelengthticks)*!values.F_NaN ;asuming not more then 50 templates
  templateflux=fltarr(50,maxwavelengthticks)*!values.F_NaN ;initialize everything with NaNs
  templatestdev=fltarr(50,maxwavelengthticks)*!values.F_NaN
  
  for i=0,numberoftemplates-1 do begin                        ;reads through each template and sets its wavelength
      restore,templatefiles(i)                                ;flux and deviation equal the appropriate values
      sizetemp(i)=n_elements(memstor.aar.data.wave)           ;in appropriate positions in the arrays
      
      templatewave(i,0:sizetemp(i)-1)=memstor.aar.data.wave
      templateflux(i,0:sizetemp(i)-1)=memstor.aar.data.flux
      templatestdev(i,0:sizetemp(i)-1)=memstor.aar.data.stdev


  endfor
  

;----------------------------------------------------------------------------------------------


; read in the extinction curves
; repeats same process documented for templates
 
  maxwavelengthtickse=0
  minwavelengthtickse=1e10 ; will be reseted in first loop (asumes there is  no template with ticks >=1e10)
  for i=0,numberofext-1 do begin

      extfiles(i)=path+extinctionfiles(i)
      restore, extfiles(i)
      wavelengthticks=n_elements(memstor.aar.data.wave)

      if wavelengthticks gt maxwavelengthtickse then begin
          maxwavelengthtickse=wavelengthticks
      endif

      if wavelengthticks lt minwavelengthtickse then begin
          minwavelengthtickse=wavelengthticks
      endif
  endfor

  extwave=fltarr(50,maxwavelengthtickse)*!values.F_NaN ;asuming not more then 50 curves
  extvalue=fltarr(50,maxwavelengthtickse)*!values.F_NaN ;initialize everything with NaNs
 
  for i=0,numberofext-1 do begin
      restore,extfiles(i)
      sizeext(i)=n_elements(memstor.aar.data.wave)
 
      extwave(i,0:sizeext(i)-1)=memstor.aar.data.wave
      extvalue(i,0:sizeext(i)-1)=memstor.aar.data.flux    
  endfor
;---------------------------------------------------------------------------------------------------

;read in source
;updated so as to be able to read both .fits files and .xdr files
;cannot read multiple .xdr files yet, only .fits
;should not need multiple .xdr files
wave=[]
flux=[]
stdev=[]
for i=0, numberofsources-1 do begin
  filename=path+sourcefiles(i) ;finds sourcefile path and restores it
  if strmatch(filename, '*gto.xdr') OR $
    strmatch(filename, '*qst.xdr') OR $
    strmatch(filename,'*prep.xdr') then begin
     restore,filename

     sourcewave=memstor.aar.data.wave ;sets the wavelengths, fluxes and deviations appropriately using memstor
     sourceflux=memstor.aar.data.flux
     sourcestdev=memstor.aar.data.stdev
     ID=sourcefiles[0]
     cluster='QUEST'
  endif else if strmatch(filename, '*ideos.xdr', /fold_case) eq 1 then begin
     restore,filename
     sourcewave=sed.wave ;sets the wavelengths, fluxes and deviations appropriately using memstor
     sourceflux=sed.flux.jy
     sourcestdev=sed.eflux.jy
     z=sed.z
     title=sed.object
     ID=sed.id
     cluster=sed.clusterid
  endif else if strmatch(filename, '*.fits', /fold_case) eq 1 then begin
     data=readfits(filename) ;read in fits file data
     ;below, wavelengths, fluxes and stdevs are concatenated into one array
     wave=[[wave],[data[0,*]]]
     flux=[[flux],[data[1,*]]]
     stdev=[[stdev],[data[2,*]]]
     ;VV below, repeated wavelengths are removed, and flux and
     ;stdev's are sorted based on indecies of new sorted wavelengths
     ;without duplicates so as to remove the stdevs and fluxes 
     ;associated with the repeated duplicates that could be in fits files
     sourcewave=wave[uniq(wave, sort(wave))]
     sourceflux=flux[uniq(wave, sort(wave))]
     sourcestdev=stdev[uniq(wave, sort(wave))]
  endif
endfor
  indzero=where(sourcestdev eq 0.) ;check for 0's in stdevfile and set them to 1. 
  if indzero(0) ne -1 then begin
     sourcestdev(indzero)=1.
     countindzero=n_elements(indzero)
     if countindzero eq n_elements(sourcestdev) then $
        print,'All error values in spectrum are 0; converting to 1.' $
     else $
        print,string(countindzero,' of ',n_elements(sourcestdev),$
                    format='(I0,A0,I0)'),$
              ' error values in spectrum are 0; converting to 1.'
  endif

  return
end
