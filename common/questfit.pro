; docformat = 'rst'
;
;+
;
; Fit MIR spectra with BB + powerlaws + templates.
;
; The superposition of individually extincted template spectra, black bodies,
; and powerlaws are fitted to a spectrum. Screen or mixed extinction can be 
; chosen. In addition, each component can be absorbed by ice features.
; The fitting parameters are: normalization factors, BB temperatures, power law
; indices, visual extinction coefficients (Av), and tau peak of the ice features.
; Each fitting parameter can in addition be set free or fixed for the fitting. 
; Weighting is done via the stdev values of the source file. stdev-values which 
; are 0 will be set to 1 apriori. The stdev's can be used to exlcude certain 
; wavelength regions from the fitting, by setting the stdev values for this 
; range high. All template spectra and extinction curves are rebined to the
; source spectrum. By manipulating the binning of the input source spectrum 
; one can influence the weighting, but this must be done before fitting.
;
; MPFIT minimizes the difference between the model and the measured data.
; SAP_REBIN rebins the real data to the model data. The rebinning is flux 
;    conserving.
; IDEOS_READCF reads in control file.
; 
; In the current version it is nessesary to put in at least one template and BB
; and use absorption for at least one of them. At least one extinction file 
; is also needed. If you don t need  e.g. absorption just force it to be 
; zero via the norm-factor.
;
; example of controlfile: test.cf 
;
;
; source                 test.xdr -1.0 35. dummy        0.0  0.0  X    0.0  0.0
; template                si1.xdr  0.4  0  LUTZ99       0.1  0.0  S    0.0  0.0
; template                si2.xdr  0.4  0  DRAINE03     0.1  0.0  S    0.0  0.0
; template                si3.xdr  0.4  0  DRAINE03     0.1  0.0  S    0.0  0.0
; absorption            H2ice.xdr  0.0  1  H2ice        0.0  0.0  X    0.0  0.0
; blackbody                  cold  0.1  0. DRAINE03     0.0  1.0  M  150.0  0.0
; blackbody                   hot  0.1  0  DRAINE03     0.0  1.0  M 1500.0  0.0
; absorption            H2ice.xdr  0.0  1  H2ice        0.0  0.0  X    0.0  0.0
; absorption           anyice.xdr  0.0  1  anyice       0.0  0.0  X    0.0  0.0
; powerlaw                  steep  0.0  1  DRAINE03     0.0  1.0  M    3.0  0.0
; powerlaw                   flat  0.0  1  DRAINE03     1.0  0.0  S    1.0  0.0
; absorption            H2ice.xdr  0.0  1  H2ice        0.0  0.0  X    0.0  0.0
; extinction         draine03.xdr  1.0  1  DRAINE03     0.0  1.0  X    0.0  0.0
; extinction           lutz99.xdr  1.0  0  LUTZ99       0.0  1.0  X    0.0  0.0
; 
; col 1: datatype
; col 2: filename (if nessesary; path hardcoded in readcf.pro)
; col 3: lower wavelength limit or normalization factor
; col 4: upper wavelength limit or fix/free parameter (1 or 0) for normalization
; col 5: name of ext. curve or ice feature
; col 6: initial guess for Av
; col 7: fix/free parameter (1/0) for Av
; col 8: S,M = screen or mixed extinction
; col 9: initial guess for BB temperature or powerlaw index
; col 10: fix/free parameter (1/0) for BB temperature or powerlaw index
;                       
; :Categories:
;    IDEOS
;
; :Returns:
;    EPS file:'fit:'+sourcename+'par:'+controlfile+'.eps'
;       Note: will name after first source listed if multiple sources are used
;    ASCII file:
;       sourcename+'par:controlfile'_BBx.fit or _PLx.fit or _templatex.fit 
;          or _total.fit'
;       ([x number of contributor] contains the fitted spectra of each 
;          contributor)
;    Data file 'fitvariables.sav' where the session data (variables)
;       are stored for further processing.
;    Data file with PAH fluxes.
;
; :Params:
;
; :Keywords:
;     controlfile: in, required, type=string
;        ASCII file containing start parameters + fit infos          
;     res: in, optional, type=byte
;        Plot residual along with fit. 
;     ps: in, optional, type=byte
;        Create EPS plot and print best fit parameters.
;     data: in, optional, type=byte
;        Output file with model spectrum.
;     log: in, optional, type=byte
;        use logarithmic axis 
;     ident: in, optional, type=byte
;        show indentification index for spectra
;     z: in, optional, type=double
;        de-redshifts the source spectra by given redshift
;     pathin: in, optional, type=string, default='./'
;     pathout: in, optional, type=string, default='./'
;     verbose: in, optional, type=byte
;        Switch on verbose output from MPFIT.
;     
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
;       2006jul06, MS, created
;       2006jul07, MS, implement mixed extinction
;       2006jul14, MS, implement absorption 
;       2006aug02, MS, implement log axis, writing spectra to files, 
;                      generate eps file, output of fitparameters on screen,
;                      integration of flux for each model, contribution (%) 
;                      plot for each model; change implementation of mixed 
;                      extinction (due to different fitresults if one uses 
;                      screen or mixed extinction but with Av forced to be 0)
;       2006aug28, MS, make use of negative PL indices possible
;       2006sep04, MS, change PL from a wavelength to a frequency powerlaw; 
;                      sorting of frequencies and spectra added for proper 
;                      integration; integration also for source spectra + 
;                      change to unit W/cm^2; add right hand index to spectra 
;                      for identification
;       2006oct30, MS, implement the ice feature wavelength range for common
;                      wavelength axis 
;       2006nov06, MS, write session data to file for further use              
;       2006nov14, MS, correction regarding the binning range
;       2006nov20, MS, read path from startup-file, write fit results into a file
;                      if ps=1, logarithm for both axis if log=1, residuum
;                      =Model/Data, new colors for eps file.
;       2013jun06, VV, z sets redshift, can accept multiple source files,
;                      can accept .fits files, removes duplicates from
;                      .fits files
;       2015nov20, DSNR, revised documentation and added copyright; changed
;                        path specification to command line
;       2016aug24, DSNR, small changes to output messages for clarity;
;                        bug fix to BB temperature lower limits to prevent
;                        floating overflow; fixed floating underflow in
;                        model function
;       2016oct25, DSNR, removed necessity of including power law
;       2021apr27, DSNR, changed fltarr --> dblarr
;    
; :Copyright:
;    Copyright (C) 2015--2021 Mario Schweitzer, Vincent Viola, David S. N. Rupke
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
;
;------------------------------------------------------------------------------
; build user defined function including BB, powerlaw, templates, mixed or
; screen extinction and absorption features 
; (maximum of all contributors normalized to 1 for easier guessing the startparameters)
; parameter order in model e.g. : 
; ext(P0)*abs(P1...)*BB(T(P2))*P3+ext(P4)*abs(P5 ...)*(X^P6)*P7+ext(P8)*abs(P9)*P10*template
; asuming one of each contributor


function modelfunc,P

common datacommon,lowerlimit,upperlimit,numberofBB,numberofpl,$
       numberoftemplates,sourcewave,sourceflux,sourcestdev,indsource,$
       indtemplate,outbinflux,parametercount,cc,numbin,count,extcurvetemp,$
       extcurvesynonym,extoutbinvalue,numberofext,countext,indextinction,$
       extcurvebb,extcurvepl,screenmixedbb,screenmixedpl,screenmixedtemplate,$
       abstempoutbintau,absbboutbintau,absploutbintau,abstempoutbinwave,$
       absbboutbinwave,absploutbinwave,numberofabsbb,numberofabspl,$
       numberofabstemp
 

model=dblarr(max(numbin))
parametercount=0

;add BB

for i=0,numberofBB-1 do begin
    
;check which extinctioncurve and which extinctionmodel is used
;edit: by matching extinctioncurve synonyms to BB, PL, and Temp
;synonyms and also checking whether extinction is mixed or screen
    if P(parametercount) ne 0 then begin
        extmodelnumberbb=where(extcurvebb(i) eq extcurvesynonym)
        extmodelnumberbbb=extmodelnumberbb(0) ; convert array(1) into integer ??? why is this nessesary
        
        if screenmixedBB(i) eq 'S' then begin
            extfactorbb=$
               10d^(-0.4d*P(parametercount)*$
                    extoutbinvalue(extmodelnumberbbb,$
                                   indextinction(extmodelnumberbbb,$
                                                 0:countext(extmodelnumberbbb)-1)))
        ;calculates deviation necessary for
        ;MPFit as well as recalculates
        ;extinctionfactor for BB
        ;taulambdabb is optical depth of black body
        endif else begin
            taulambdabb=P(parametercount)*0.4d*alog(10d)*$
                        extoutbinvalue(extmodelnumberbbb,$
                                       indextinction(extmodelnumberbbb,$
                                                     0:countext(extmodelnumberbbb)-1))
            extfactorbb=(1d -exp(-taulambdabb))/taulambdabb
        endelse
    endif else begin
        extfactorbb=1d
    endelse
    
; estimate factor for absorption
; calculates change in intensity of source wave due to absorption
    
    absfactorbb=1.
    
    for tt=0,numberofabsbb(i)-1 do begin
        
        parametercount=parametercount+1
        indexbb=where(absbboutbinwave(i,tt,*) ge lowerlimit and absbboutbinwave(i,tt,*) le upperlimit)
;    The "mask" parameter eliminates floating underflow by removing very large
;    negative exponents. See http://www.idlcoyote.com/math_tips/underflow.html
;    for more details.
        exparg = -P(parametercount)*absbboutbintau(i,tt,indexbb)
        mask = (abs(exparg) lt 80)
        absfactorbb=absfactorbb*mask*exp(exparg*mask) 
        
    endfor
    
    parametercount=parametercount+1
    BB=planck(sourcewave(indsource)*10d^4d,P(parametercount)) * $
       ((sourcewave(indsource)*10d^4d)^2d/cc) 
    parametercount=parametercount+1
    model=model+P(parametercount)*absfactorbb*extfactorbb*BB/max(BB)
    parametercount=parametercount+1
    
endfor

;add powerlaws 
;repeat blackbody code for powerlaws

for i=0,numberofpl-1 do begin
    
;check which extinctioncurve and which extinctionmodel is used
    if P(parametercount) ne 0. then begin
        extmodelnumberpl=where(extcurvepl(i) eq extcurvesynonym)
        extmodelnumberplb=extmodelnumberpl(0) ; convert array(1) into integer ??? why is this nessesary
        if screenmixedpl(i) eq 'S' then begin 

        extfactorpl=10.d^((-0.4d)*P(parametercount)*extoutbinvalue(extmodelnumberplb,indextinction (extmodelnumberplb,0:countext(extmodelnumberplb)-1))) 

        endif else begin

            
            taulambdapl=P(parametercount)*0.4*alog(10)*extoutbinvalue(extmodelnumberplb,indextinction (extmodelnumberplb,0:countext(extmodelnumberplb)-1))
            extfactorpl=(1-exp(-taulambdapl))/taulambdapl

        endelse    

    endif else begin

        extfactorpl=1.

    endelse

; estimate factor for absorption

    absfactorpl=1.
    
    for tt=0,numberofabspl(i)-1 do begin
        
        parametercount=parametercount+1
        indexpl=where(absploutbinwave(i,tt,*) ge lowerlimit and absploutbinwave(i,tt,*) le upperlimit)
        absfactorpl=absfactorpl*exp(-P(parametercount)*absploutbintau(i,tt,indexpl)) 
        
    endfor
    
    parametercount=parametercount+1
    pl=double((cc/(sourcewave(indsource)*10.d^(-6)))^P(parametercount)) ; powerlaw for frequency
    parametercount=parametercount+1
    model=model+P(parametercount)*absfactorpl*extfactorpl*pl/max(pl)
    parametercount=parametercount+1
    
endfor

;add templates
;repeat code for blackbody and powerlaw for templates


for i=0,numberoftemplates-1 do begin
    
;check which extinctioncurve and which extinctionmodel is used 
    
    if P(parametercount) ne 0. then begin   
        
        extmodelnumber=where(extcurvetemp(i) eq extcurvesynonym)
        extmodelnumberb=extmodelnumber(0) ; convert array(1) into integer ??? why is this nessesary
        
        if screenmixedtemplate(i) eq 'S' then begin
            extfactor=10.d^((-0.4d)*P(parametercount)*extoutbinvalue(extmodelnumberb,indextinction (extmodelnumberb,0:countext(extmodelnumberb)-1)))
        endif else begin
            taulambdatemp=P(parametercount)*0.4*alog(10)*extoutbinvalue(extmodelnumberb,indextinction (extmodelnumberb,0:countext(extmodelnumberb)-1))
            extfactor=(1-exp(-taulambdatemp))/taulambdatemp
        endelse
    endif else begin
        extfactor=1.
    endelse
    
    
; estimate factor for absorption
    
    absfactortemp=1.
    
    for tt=0,numberofabstemp(i)-1 do begin
        
        parametercount=parametercount+1
        indextemp=where(abstempoutbinwave(i,tt,*) ge lowerlimit and abstempoutbinwave(i,tt,*) le upperlimit)
        absfactortemp=absfactortemp*exp(-P(parametercount)*abstempoutbintau(i,tt,indextemp)) 
        
    endfor
    
    
    parametercount=parametercount+1
    model=model+absfactortemp*extfactor*P(parametercount)*outbinflux(i,indtemplate(i,0:count(i)-1))/max(outbinflux(i,indtemplate(i,0:count(i)-1)),/NaN)
    parametercount=parametercount+1
    
endfor


;creates final return value for P to go to MPFit
chi=(sourceflux(indsource)-model)/sourcestdev(indsource) 

return, chi

end


;------------------------------------------------------------------------------

pro questfit,controlfile=controlfile,res=res,ps=ps,data=data,log=log,$
             ident=ident, ztweak=ztweak, contrib=contrib, pathin=pathin, $
             pathout=pathout,verbose=verbose

common datacommon,lowerlimit,upperlimit,numberofBB,numberofpl,$
       numberoftemplates,sourcewave,sourceflux,sourcestdev,indsource,$
       indtemplate,outbinflux,parametercount,cc,numbin,count,extcurvetemp,$
       extcurvesynonym,extoutbinvalue,numberofext,countext,indextinction,$
       extcurvebb,extcurvepl,screenmixedbb,screenmixedpl,$
       screenmixedtemplate,abstempoutbintau,absbboutbintau,$
       absploutbintau,abstempoutbinwave,absbboutbinwave,absploutbinwave,$
       numberofabsbb,numberofabspl,numberofabstemp

;constants

cc=299792458.d  ; [m/s](light speed)

if not keyword_set(pathin) then pathin='./'
if not keyword_set(pathout) then pathout='./'

; read in sourcedata,templatedata,BB+PL data,range ... from controlfile
; see readcf.pro for more info on variable names
 
ideos_readcf,controlfile=controlfile,numberofBB,numberofPL,numberoftemplates,$
             numberofext,sourcefiles,range,templatefiles,templatefilesb,normtemp,$
             fixfreetemp,extcurvetemp,extvaluetemp,normbb,fixfreebb,fixfreebbtemp,$
             extcurvebb,extvaluebb,tempbb,normpl,fixfreepl,fixfreeplindex,$
             extcurvepl,extvaluepl,indexpl,templatewave,templateflux,sourcewave,$
             sourceflux,sourcestdev,maxwavelengthticks,sizetemp,sizeext,extwave,$
             extvalue,extvaluefixfreetemp,extvaluefixfreebb,extvaluefixfreepl,$
             extcurvesynonym,screenmixedBB,screenmixedpl,screenmixedtemplate,$
             taupeaktemp,taupeakbb,taupeakpl,taupeaktempfixfree,taupeakbbfixfree,$
             taupeakplfixfree, numberofabstemp, numberofabsbb, numberofabspl,$
             absorbtempwave,absorbbbwave,absorbplwave,absorbtemptau,absorbbbtau,$
             absorbpltau,abstempsize,absbbsize,absplsize,z,title,ID,cluster,$
             pathin=pathin

; check the wavelength-boundaries in extcurve,absfeature,templates and source
; => wavelengthrange (upperlimit,lowerlimit)

;------------------------------------------------------------------------------
;de-redshift spectra if necessary

if keyword_set(ztweak) then sourcewave=sourcewave/(ztweak+1d)

;goodindex=where(sourcewave lt 7.6 or sourcewave gt 7.7)
;sourcewave=sourcewave[goodindex]
;sourceflux=sourceflux[goodindex]
;sourcestdev=sourcestdev[goodindex]
;goodindex=where(sourcewave lt 10.4 or sourcewave gt 10.6)
;sourcewave=sourcewave[goodindex]
;sourceflux=sourceflux[goodindex]
;sourcestdev=sourcestdev[goodindex]
;goodindex=where(sourcewave lt 12.7 or sourcewave gt 12.9)
;sourcewave=sourcewave[goodindex]
;sourceflux=sourceflux[goodindex]
;sourcestdev=sourcestdev[goodindex]
;goodindex=where(sourcewave lt 14.2 or sourcewave gt 14.4)
;sourcewave=sourcewave[goodindex]
;sourceflux=sourceflux[goodindex]
;sourcestdev=sourcestdev[goodindex]
;goodindex=where(sourcewave lt 15.5 or sourcewave gt 15.7)
;sourcewave=sourcewave[goodindex]
;sourceflux=sourceflux[goodindex]
;sourcestdev=sourcestdev[goodindex]
;goodindex=where(sourcewave lt 16.9 or sourcewave gt 17.2)
;sourcewave=sourcewave[goodindex]
;sourceflux=sourceflux[goodindex]
;sourcestdev=sourcestdev[goodindex]
;goodindex=where(sourcewave lt 25 or sourcewave gt 26)
;sourcewave=sourcewave[goodindex]
;sourceflux=sourceflux[goodindex]
;sourcestdev=sourcestdev[goodindex]
;in source
;lower boundary

minsourcewave=min(sourcewave,/NaN) 

;find smallest common wavelengthvalue in templates

mintempwave=dblarr(numberoftemplates)

for i=0,numberoftemplates-1 do begin

    mintempwave(i)=min(templatewave(i,*),/NaN) 
    
endfor

mintemplatewave=max(mintempwave,/NaN)

;find smallest common wavelengthvalue in extinctioncurves

minextwave=dblarr(numberofext)

for i=0,numberofext-1 do begin

    minextwave(i)=min(extwave(i,*),/NaN) 
    
endfor

minextinctionwave=max(minextwave,/NaN)


;find smallest common wavelengthvalue for absorptionfeature for BB
;(a little bit complicated :-( )

;copy arrays for processing
absorbbbwavecopy=absorbbbwave
absorbtempwavecopy=absorbtempwave
if numberofpl gt 0 then absorbplwavecopy=absorbplwave $
else minabsorbpltotal=0d

for i=0,numberofbb-1 do begin
    for j=0,numberofabsbb(i)-1 do begin
        index=where(absorbbbwavecopy(i,j,*) eq 0.)
        absorbbbwavecopy (i,j,index)=1.e6
    endfor
endfor
minimuma=dblarr(1000,1000)
for i=0,numberofbb-1 do begin
    for j=0,numberofabsbb(i)-1 do begin
        minimuma(i,j)=min(absorbbbwavecopy(i,j,*)) 
     endfor
endfor
minabsorbbbtotal=max(minimuma)

;find smallest common wavelengthvalue for absorptionfeature for template
for i=0,numberoftemplates-1 do begin
    for j=0,numberofabstemp(i)-1 do begin
        index=where(absorbtempwavecopy(i,j,*) eq 0.)
        absorbtempwavecopy (i,j,index)=1.e6
    endfor
endfor

minimumb=dblarr(1000,1000)

for i=0,numberoftemplates-1 do begin
    for j=0,numberofabstemp(i)-1 do begin
        minimumb(i,j)=min(absorbtempwavecopy(i,j,*)) 
     endfor
endfor
minabsorbtemptotal=max(minimumb)

;find smallest common wavelengthvalue for absorptionfeature for PL
for i=0,numberofpl-1 do begin
    for j=0,numberofabspl(i)-1 do begin
        index=where(absorbplwavecopy(i,j,*) eq 0.)
        absorbplwavecopy (i,j,index)=1.e6
    endfor
endfor

minimumc=dblarr(1000,1000)

for i=0,numberofpl-1 do begin
    for j=0,numberofabspl(i)-1 do begin
        minimumc(i,j)=min(absorbplwavecopy(i,j,*)) 
     endfor
endfor
minabsorbpltotal=max(minimumc)

;estimate smallest common wavelength for all absorption features:

commonabsmin=max([minabsorbbbtotal,minabsorbtemptotal,minabsorbpltotal])

;-------

minofall=max([minsourcewave,mintemplatewave,minextinctionwave,commonabsmin],/NaN) ;smallest common wavelength  of all data


;------

if range(0) eq -1 then begin   ;this if-else block checks to see if the lower wavelength
    lowerlimit=minofall        ;limit set in the control file is too small
endif else begin               ;and sets it to smallest values if sentinal '-1' is given
    lowerlimit=range(0)
    
    if range(0) lt minofall then begin
        print,'lowerlimit lt smallest possible value :',minofall
        stop
    endif

endelse


;upper boundary
;----------------------------------------------------------------

;repeat min code but for max
;in source

maxsourcewave=max(sourcewave,/NaN)

;find largest common wavelengthvalue in templates

maxtempwave=dblarr(numberoftemplates)

for i=0,numberoftemplates-1 do begin
    
    maxtempwave(i)=max(templatewave(i,*),/NaN)

endfor

maxtemplatewave=min(maxtempwave,/NaN)

;find largest common wavelengthvalue for absorptionfeature for BB

maxabsorbBB=max(absorbbbwave,dimension=3)
indbbabsorb=where(maxabsorbBB gt 0.)
maxabsorbBBtotal=min(maxabsorbBB(indbbabsorb))

;find largest common wavelengthvalue for absorptionfeature for PL

if numberofpl gt 0 then begin
   maxabsorbPL=max(absorbplwave,dimension=3)
   indplabsorb=where(maxabsorbpl gt 0.)
   maxabsorbpltotal=min(maxabsorbpl(indplabsorb))
endif else maxabsorbpltotal=999d

;find largest common wavelengthvalue for absorptionfeature for template

maxabsorbtemp=max(absorbtempwave,dimension=3)
indtempabsorb=where(maxabsorbtemp gt 0.)
maxabsorbtemptotal=min(maxabsorbtemp(indtempabsorb))


;estimate greatest common wavelength for all absorption features:

commonabsmax=min([maxabsorbBBtotal,maxabsorbtemptotal,maxabsorbpltotal])


;find largest common wavelengthvalue in extinctioncurves

maxextwave=dblarr(numberofext)

for i=0,numberofext-1 do begin
    
    maxextwave(i)=max(extwave(i,*),/NaN)

endfor

maxextinctionwave=min(maxextwave,/NaN)

;----
maxofall=min([maxsourcewave,maxtemplatewave,maxextinctionwave,commonabsmax],/NaN);greatest common wavelength of all data
;----

if range(1) eq -1 then begin  ;does same process as when lower limit was set
    upperlimit=maxofall
endif else begin
    upperlimit=range(1)

    if range(1) gt maxofall then begin
        print,'upperlimit gt largest possible value :',maxofall
        stop
    endif
    
endelse




;index of wavelengthboundaries in source data

indsource=where (sourcewave ge lowerlimit and sourcewave le upperlimit)

;if indsource equals 0, you will have an error later so check and make sure
;your sourcewave is less than your upper limit and greater than your
;lower limit.  you can do this by restoring your source file and
;printing the value of 'memstor.aar.data.wave'


;rebinning
;---------------------------------------------------------------------------


;absorptionfeatures

;define output variables

abstempoutbinwave=dblarr(50,50,n_elements(sourcewave)) ;not more then 50 templates and 50 abs-features per template !
abstempoutbintau=dblarr(50,50,n_elements(sourcewave))
abstempnumbin=dblarr(numberoftemplates,total(numberofabstemp))

for i =0,numberoftemplates-1 do begin
    
    for ii=0,numberofabstemp(i)-1 do begin

;define and reset variables

        tpltwave=dblarr(abstempsize(i,ii))  ; creates an array for the purpose of inputting 
        tpltvalue=dblarr(abstempsize(i,ii)) ; all data for absorption features into sap_rebin to be rebinned
        tpltstdev=dblarr(abstempsize(i,ii)) ; just 0's
        
        tpltwave(*)=absorbtempwave(i,ii,0:abstempsize(i,ii)-1) ;conversion of 3d inputvalues to 1d array
        tpltvalue(*)=absorbtemptau(i,ii,0:abstempsize(i,ii)-1) ;because MPFit and sap_rebin only take 1d arrays
        

        statusabstemp=sap_rebin(tpltwave,tpltvalue,tpltstdev,absbinwave,absbintau,absbinstdev,ref=sourcewave(indsource)) ;absorption data is rebinned at this point to fit the source spectra

        abstempnumbin(i,ii)=n_elements(absbinwave)
        
        abstempoutbinwave(i,ii,0:abstempnumbin(i,ii)-1)=absbinwave ;reconversion of outputvalues to 3d array
        abstempoutbintau(i,ii,0:abstempnumbin(i,ii)-1)=absbintau
 
    endfor    

endfor



;BB

;repeat rebinning process for absorption
;define output variables

absBBoutbinwave=dblarr(50,50,n_elements(sourcewave)) ;not more then 50 templates and 50 abs-features per template !
absBBoutbintau=dblarr(50,50,n_elements(sourcewave))
absBBnumbin=dblarr(numberofBB,total(numberofabsBB))

for i =0,numberofBB-1 do begin

    for ii=0,numberofabsBB(i)-1 do begin

        ;define and reset variables
        tpltwave=dblarr(absBBsize(i,ii))
        tpltvalue=dblarr(absBBsize(i,ii))
        tpltstdev=dblarr(absBBsize(i,ii)) ; just 0's
        
        tpltwave(*)=absorbbbwave(i,ii,0:absBBsize(i,ii)-1) ;conversion of 3d inputvalues to 1d array
        tpltvalue(*)=absorbbbtau(i,ii,0:absBBsize(i,ii)-1)
   
        statusabsbb=sap_rebin(tpltwave,tpltvalue,tpltstdev,absbinwave,absbintau,absbinstdev,ref=sourcewave(indsource))

        absBBnumbin(i,ii)=n_elements(absbinwave)
        
        absBBoutbinwave(i,ii,0:absBBnumbin(i,ii)-1)=absbinwave ;reconversion of outputvalues to 3d array
        absBBoutbintau(i,ii,0:absBBnumbin(i,ii)-1)=absbintau
        
    endfor    
    
endfor

;PL

;repeat rebinning for power law data
;define output variables

if numberofpl gt 0 then begin
  absploutbinwave=dblarr(50,50,n_elements(sourcewave)) ;not more then 50 templates and 50 abs-features per template !
  absploutbintau=dblarr(50,50,n_elements(sourcewave))
  absplnumbin=dblarr(numberofpl,total(numberofabspl))

  for i =0,numberofpl-1 do begin

    for ii=0,numberofabspl(i)-1 do begin
      ;define and reset variables
      tpltwave=dblarr(absplsize(i,ii))
      tpltvalue=dblarr(absplsize(i,ii))
      tpltstdev=dblarr(absplsize(i,ii)) ; just 0's
    
      tpltwave(*)=absorbplwave(i,ii,0:absplsize(i,ii)-1) ;conversion of 3d inputvalues to 1d array
      tpltvalue(*)=absorbpltau(i,ii,0:absplsize(i,ii)-1)
   
      statusabspl=sap_rebin(tpltwave,tpltvalue,tpltstdev,absbinwave,absbintau,absbinstdev,ref=sourcewave(indsource))

      absplnumbin(i,ii)=n_elements(absbinwave)

      absploutbinwave(i,ii,0:absplnumbin(i,ii)-1)=absbinwave ;reconversion of outputvalues to 3d array
      absploutbintau(i,ii,0:absplnumbin(i,ii)-1)=absbintau
 
    endfor    
  endfor
endif


;rebin extinctioncurves to source
;---------------------------------------------------------------------------

;same as previous rebinning
;define output variables

extoutbinwave=dblarr(50,n_elements(sourcewave))
extoutbinvalue=dblarr(50,n_elements(sourcewave))
extnumbin=dblarr(numberofext)

for i=0,numberofext-1 do begin
    ;define and reset variables
    tpltwave=dblarr(sizeext(i))
    tpltvalue=dblarr(sizeext(i))
    tpltstdev=dblarr(sizeext(i)); just 0's
    
    tpltwave(*)=extwave(i,0:sizeext(i)-1) ;conversion of 2d inputvalues to 1d array
    tpltvalue(*)=extvalue(i,0:sizeext(i)-1)
   

    statusext=sap_rebin(tpltwave,tpltvalue,tpltstdev,extbintemplatewave,extbintemplatevalue,extbintemplatestdev,ref=sourcewave(indsource))
    extnumbin(i)=n_elements(extbintemplatewave)

    extoutbinwave(i,0:extnumbin(i)-1)=extbintemplatewave ;reconversion of outputvalues to 2d array
    extoutbinvalue(i,0:extnumbin(i)-1)=extbintemplatevalue
    
endfor

;-----------------------------------------------------------------------------



;rebin templates to source
;---------------------------------------------------------------------------

;same as previous rebinning
;define output variables

outbinwave=dblarr(50,n_elements(sourcewave))
outbinflux=dblarr(50,n_elements(sourcewave))
outbinstdev=dblarr(50,n_elements(sourcewave))
numbin=dblarr(numberoftemplates)

for i=0,numberoftemplates-1 do begin
    ;define and reset variables
    tpltwave=dblarr(sizetemp(i))
    tpltflux=dblarr(sizetemp(i))
    tpltstdev=dblarr(sizetemp(i)); just 0's
    
    tpltwave(*)=templatewave(i,0:sizetemp(i)-1) ;conversion of 2d inputvalues to 1d array
    tpltflux(*)=templateflux(i,0:sizetemp(i)-1)
    
    statustemp=sap_rebin(tpltwave,tpltflux,tpltstdev,bintemplatewave,bintemplateflux,bintemplatestdev,ref=sourcewave(indsource))

    numbin(i)=n_elements(bintemplatewave)

    outbinwave(i,0:(numbin(i)-1))=bintemplatewave ;reconversion of outputvalues to 2d array
    outbinflux(i,0:(numbin(i)-1))=bintemplateflux
 

endfor


;-----------------------------------------------------------------------------


;-------------------------------------------------------------------------------------


; find indices in spectra,  and extinction curves for the above defined wavelengthrange
; and set all values in rebined spectra and ext.curves not belonging to this
; range to NaNs. (Indexrange for absorption included in userdefined function)
;-----------------------------------------------------------------------------------------

;In english: remove all values outside of wavelength range
;index of wavelengthboundaries in template data

indtemplate=dblarr(50,max(numbin))
count=dblarr(50)


for i=0,numberoftemplates-1 do begin

    index=where(outbinwave(i,*) ge lowerlimit and outbinwave(i,*) le upperlimit)
    count(i)=n_elements(index)

    indtemplate(i,0:count(i)-1)=where(outbinwave(i,*) ge lowerlimit and outbinwave(i,*) le upperlimit)
    ; ^^proper values are indexed and counted
    ; set all other values in outbinwave,flux,stdev to NaN
    ; in other words, all values outside of the lowerlimit and the upperlimit
    ; are removed
    indexelsewhere=where(outbinwave(i,*) lt lowerlimit or  outbinwave(i,*) gt upperlimit)

    if indexelsewhere(0) ne -1 then begin

        outbinwave(i,indexelsewhere)=!values.F_NaN
        outbinflux(i,indexelsewhere)=!values.F_NaN
        outbinstdev(i,indexelsewhere)=!values.F_NaN
   
    endif
    
endfor


;index of wavelengthboundaries in extinctioncurves

indextinction=dblarr(50,max(extnumbin))
countext=dblarr(50)


for i=0,numberofext-1 do begin

    index=where(extoutbinwave(i,*) ge lowerlimit and extoutbinwave(i,*) le upperlimit)
    countext(i)=n_elements(index)

    indextinction(i,0:countext(i)-1)=where(extoutbinwave(i,*) ge lowerlimit and extoutbinwave(i,*) le upperlimit)

    ; set all other values in extoutbinwave,value to NaN
    ; same as removing values in previous block
    indexelsewhereext=where(extoutbinwave(i,*) lt lowerlimit or  extoutbinwave(i,*) gt upperlimit)

    if  indexelsewhereext(0) ne -1 then begin
        extoutbinwave(i,indexelsewhereext)=!values.F_NaN
        extoutbinvalue(i,indexelsewhereext)=!values.F_NaN
    endif


endfor



;----------------------------------------------------------------------------


;set fixfreeparameters and startparameters for BB, PL and templates 
;
;---------------------------------------------------------------------------------------------------------------------------------

; 2 arrays are created below: one for the fixfree prameter for
; blackbodies, powerlaws and templates, and the other is the array for
; the start values of the black bodies, powerlaws and templates
fixfreeblackbody=dblarr(3*numberofBB+total(numberofabsbb)) 
if numberofpl gt 0 then fixfreepowerlaw=dblarr(3*numberofpl+total(numberofabspl))
fixfreetemplates=dblarr(2*numberoftemplates+total(numberofabstemp))

pstartBB=dblarr(3*numberofBB+total(numberofabsbb))
if numberofpl gt 0 then pstartpl=dblarr(3*numberofpl+total(numberofabspl))
pstarttemplates=dblarr(2*numberoftemplates+total(numberofabstemp))

if numberofpl gt 0 then indexparametercount=dblarr(numberofpl) $
else indexparametercount=0

;BB 
;the fixfree value and the start value for each black body is indexed
;and recorded in an array here
c=0
for i=0,numberofBB-1 do begin               
countbbabsorb=0    
    pstartBB(c)=extvaluebb(i)
    fixfreeblackbody(c)=extvaluefixfreeBB(i)
    c=c+1

    for ii=0,numberofabsbb(i)-1 do begin
        pstartBB(c)=taupeakBB(i,countbbabsorb)
        fixfreeblackbody(c)=taupeakBBfixfree(i,countbbabsorb)
        countbbabsorb=countbbabsorb+1   
        c=c+1 
    endfor

    pstartBB(c)=tempbb(i)
    fixfreeblackbody(c)=fixfreebbtemp(i)
    c=c+1
    pstartBB(c)=normbb(i)
    fixfreeblackbody(c)=fixfreebb(i)
    c=c+1    
    
endfor

;PL 
;same for powerlaw array as blackbodies
c=0 
for i=0,numberofpl-1 do begin 
countplabsorb=0

    pstartPL(c)=extvaluepl(i)
    fixfreepowerlaw(c)=extvaluefixfreepl(i)
    c=c+1

    for ii=0,numberofabspl(i)-1 do begin
        pstartpl(c)=taupeakpl(i,countplabsorb)
        fixfreepowerlaw(c)=taupeakplfixfree(i,countplabsorb)
        countplabsorb=countplabsorb+1   
        c=c+1 
    endfor

    pstartPL(c)=indexpl(i)
    fixfreepowerlaw(c)=fixfreeplindex(i)
    indexparametercount(i)=c+1 ; remember the number of index-parameter for later setting the parameterlimits separately
    c=c+1
    pstartPL(c)=normpl(i)
    fixfreepowerlaw(c)=fixfreepl(i)
    c=c+1
endfor



c=0
;Template 
;same for template array as black body and power law
for i=0,numberoftemplates-1 do begin
counttempabsorb=0

    pstarttemplates(c)=extvaluetemp(i)
    fixfreetemplates(c)=extvaluefixfreetemp(i)
    c=c+1

   for ii=0,numberofabstemp(i)-1 do begin
        pstarttemplates(c)=taupeaktemp(i,counttempabsorb)
        fixfreetemplates(c)=taupeaktempfixfree(i,counttempabsorb)
        counttempabsorb=counttempabsorb+1   
        c=c+1 
    endfor

    pstarttemplates(c)=normtemp(i)
    fixfreetemplates(c)=fixfreetemp(i)
    c=c+1
endfor








; create start and fixfree array
;--------------------------------------------------------------------------------------------

; all start values are combined into one array
if numberofpl gt 0 then pstart=double([pstartBB,pstartpl,pstarttemplates]) $
else pstart=double([pstartBB,pstarttemplates])

parinfo=replicate({value:0.D,fixed:0,limited:[0,0],limits:[0.D,0.D]},2*numberoftemplates+total(numberofabstemp)+3*numberofBB+total(numberofabsbb)+3*numberofPL+total(numberofabspl))

if numberofpl gt 0 then $
  fixfreetotal=[fixfreeblackbody,fixfreepowerlaw,fixfreetemplates] $
else fixfreetotal=[fixfreeblackbody,fixfreetemplates]

parinfo(*).fixed=fixfreetotal(*)

;limit parameters to be ge 0 (powerlaw index free)

parinfo(*).limited(0)=1b
parinfo(*).limits(0)=0d ; all parameters limited to be gt 0
ind_par_bbtemp = intarr(numberofBB)
ind_par_bbtemp[0] = 1 + numberofabsBB[0]
for i=1,numberofBB-1 do $
   ind_par_bbtemp[i] = ind_par_bbtemp[i-1] + 3 + numberofabsbb[i]
; lower limit for BB temperature must be > 0, or call to PLANCK chokes with
; divide by 0
parinfo(ind_par_bbtemp).limits(0)=1d
parinfo(3*numberofBB+total(numberofabsbb)+indexparametercount-1).limited(0)=0b ; no limits for the PL indexparameter

; Parinfo is a specific structure for use in MPFit.  Similarly, value,
; fixed, limited and limits are specific tags for MPFit.  See MPFit
; documentation for specific meanings.  General idea is that the
; fixfree values are set to specific replications of fixed for use in
; parinfo, then all of the first values of limited are set to zero to
; bind all parameters on the lower side.  The lower bound is then set
; to zero by using the limits tag.  Finally, this is undone for some
; values so that there are no powerlaw limits for its indexparameter

;do the fit

;--------------------------------------------------

if keyword_set(verbose) then quiet=0 else quiet=1
result=MPFIT('modelfunc',pstart,STATUS=status,ERRMSG=errmsg,PARINFO=parinfo,$
             BESTNORM=chisq,quiet=quiet)
print, errmsg

;--------------------------------------------------------------------------------------


;create eps file

if keyword_set(ps) then begin
    set_plot,'ps'
    epsfile = ID+'_'+cluster+'.eps'
    device,/portrait,/color,/encapsulated,file=pathout+epsfile,xsize=25,ysize=45

    loadct,5,/silent
    red=  [0,1,1,0,0.0,0.8]
    green=[0,1,0,0,0.7,0.4]
    blue= [0,1,0,1,0.0,0.0]
    tvlct,255*red,255*green,255*blue
    colora=2                  ;sum
    colorb=[cgcolor('green'),cgcolor('purple'),cgcolor('brown'),cgcolor('orange')]                  ;BB
    colorc=5                  ;Pl
    colord=[cgcolor('blue'),cgcolor('dark green'),cgcolor('gray'),$
      cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),$
      cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),$
      cgcolor('gray')]               ;template

endif

colora=cgcolor('red')                  ;sum
colorb=[cgcolor('green'),cgcolor('purple'),cgcolor('brown'),cgcolor('orange')]            ;BB
colorc=5                  ;Pl
    colord=[cgcolor('blue'),cgcolor('dark green'),cgcolor('gray'),$
      cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),$
      cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),cgcolor('gray'),$
      cgcolor('gray')]               ;template

!P.multi=[0,2,2]



;define flux limits for plot

;plot range

upperlimitflux=max(sourceflux(indsource)) ; takes flux values from source file and finds minimum value
lowerlimitflux=min(sourceflux(indsource)) ; same for max flux value

rangex=[lowerlimit,upperlimit]
rangey=[0.,upperlimitflux]     ;cannot have flux less than zero

;type of axes are set below based on keywords
;plots upper and lower flux limits on x-axis and just 0 to upper limit
;on y-axis

if keyword_set(res) then begin

    rangey=[-0.5*upperlimitflux,upperlimitflux] 
    ynull=0.65

    if keyword_set(log) then begin
         plot,rangex,rangey, XTickformat='(A1)', pos=[0.2, ynull, 0.95, 0.95],xstyle=1,ystyle=1,psym=1,ytitle='flux [Jy]',title=title,xrange=[lowerlimit, upperlimit],yrange=[0.0001, upperlimitflux],/ylog,/xlog

    endif else begin

         plot,rangex,rangey, XTickformat='(A1)', pos=[0.2, ynull, 0.95, 0.95],xstyle=1,ystyle=1,psym=1,ytitle='flux [Jy]',title=title,xrange=[lowerlimit, upperlimit],yrange=[0.0001, upperlimitflux]

    endelse
  
endif else begin

    ynull=0.45

    if keyword_set(log) then begin

         plot,rangex,rangey, pos=[0.2, ynull, 0.95, 0.95],xstyle=1,ystyle=1,psym=1,xtitle='rest-wavelength [�m]',ytitle='flux [Jy]',title=title,xrange=[lowerlimit, upperlimit],yrange=[0.0001, upperlimitflux],/ylog,/xlog
    endif else begin

     plot,rangex,rangey, pos=[0.2, ynull, 0.95, 0.95],xstyle=1,ystyle=1,psym=1,xtitle='rest-wavelength [�m]',ytitle='flux [Jy]',title=title,xrange=[lowerlimit, upperlimit],yrange=[0.0001, upperlimitflux]

    endelse
endelse



;sum everything up and plot components!!!
;---

;add+plot BBs 

parametercount=0
P=result ; use fitresult for parameters

 
;-----------------------------------------
; creates arrays with the dimensions of a row per BB, PL or Temp by
; the number of elements in the wave array from the source file
sum=dblarr(max(numbin))
BBmem=dblarr(numberofBB,n_elements(sourcewave(indsource)))
BBmemnorm=dblarr(numberofBB)
if numberofpl gt 0 then begin
   PLmem=dblarr(numberofpl,n_elements(sourcewave(indsource)))
   PLmemnorm=dblarr(numberofpl)
endif
Tempmem=dblarr(numberoftemplates,n_elements(sourcewave(indsource)))
Tempmemnorm=dblarr(numberoftemplates)

for i=0,numberofBB-1 do begin
    
;check which extinctioncurve and which extinctionmodel is used
;edit: by matching extinctioncurve synonyms to BB, PL, and Temp
;synonyms and also checking whether extinction is mixed or screen
   
    if P(parametercount) ne 0. then begin
        
        extmodelnumberbb=where(extcurvebb(i) eq extcurvesynonym)
        extmodelnumberbbb=extmodelnumberbb(0) ; convert array(1) into integer ??? why is this nessesary
        
        if screenmixedBB(i) eq 'S' then begin
       
        ;calculates the extinction factor in the blackbody     
        extfactorbb=10.d^((-0.4d)*P(parametercount)*extoutbinvalue(extmodelnumberbbb,indextinction (extmodelnumberbbb,0:countext(extmodelnumberbbb)-1)))

        endif else begin

        ;calculates deviation necessary for
        ;MPFit as well as recalculates
        ;extinctionfactor for BB
        ;taulambdabb is optical depth of black body
        taulambdabb=P(parametercount)*0.4*alog(10)*extoutbinvalue(extmodelnumberbbb,indextinction (extmodelnumberbbb,0:countext(extmodelnumberbbb)-1))
        extfactorbb=(1-exp(-taulambdabb))/taulambdabb

          endelse

    endif else begin 

        extfactorbb=1.

    endelse
    
;------- 

;factor for absorption

    absfactorbb=1.    
    
    for tt=0,numberofabsbb(i)-1 do begin
        
       parametercount=parametercount+1
       indexbb=where (absbboutbinwave(i,tt,*) ge lowerlimit and absbboutbinwave(i,tt,*) le upperlimit)
       exparg = -P(parametercount)*absbboutbintau(i,tt,indexbb)
;      The "mask" parameter eliminates floating underflow by removing very large
;      negative exponents. See http://www.idlcoyote.com/math_tips/underflow.html
;      for more details.
       mask = (abs(exparg) lt 80)
       absfactorbb=absfactorbb*mask*exp(mask*exparg) 
        
    endfor
    ; ^^indexing absorbtion for black bodies
;-------
    

    ; below, the blackbody effect is
    ; calculated and plotted and fitted
    ; lots of mathmathmath

    parametercount=parametercount+1
    BB=planck(sourcewave(indsource)*10.d^4.d,P(parametercount))*((sourcewave(indsource)*10d^4.d)^2.d/cc) 
    parametercount=parametercount+1
    sum=sum+P(parametercount)*absfactorbb*extfactorBB*BB/max(BB)
    cgoplot, sourcewave(indsource),P(parametercount)*absfactorbb*extfactorbb*BB/max(BB),color=colorb[i],linestyle=i
    BBmem(i,0:n_elements(BB)-1)=P(parametercount)*absfactorbb*extfactorBB*BB/max(BB) ; store the single BB's
    BBmemnorm(i)=P(parametercount)/max(BB) ; store the normalizationfactor


       if keyword_set(ident) then begin 
          if P(parametercount) ne 0. then begin
            convert=strtrim(string(i+1),2) 
            maxi=max(sourcewave(indsource))
            xyouts,maxi,BBmem(i,n_elements(BB)-1),'BB'+convert
        endif
    endif

    parametercount=parametercount+1

endfor



;add+plot powerlaws 

for i=0,numberofpl-1 do begin
    
;check which extinctioncurve is used
;same process done for blackbodies is repeated for power laws
    if P(parametercount) ne 0. then begin
        
        extmodelnumberpl=where(extcurvepl(i) eq extcurvesynonym)
        extmodelnumberplb=extmodelnumberpl(0) ; convert array(1) into integer ??? why is this nessesary
        
        if screenmixedpl(i) eq 'S' then begin            
            
            extfactorpl=10.d^((-0.4d)*P(parametercount)*extoutbinvalue(extmodelnumberplb,indextinction (extmodelnumberplb,0:countext(extmodelnumberplb)-1))) ;mathmathmath
        endif else begin
            taulambdapl=P(parametercount)*0.4*alog(10)*extoutbinvalue(extmodelnumberplb,indextinction (extmodelnumberplb,0:countext(extmodelnumberplb)-1))  ;mathmathmath
            extfactorpl=(1-exp(-taulambdapl))/taulambdapl
        endelse  
        
    endif else begin
        extfactorpl=1.
    endelse

;-------

; estimate factor for absorption

     absfactorpl=1.

     for tt=0,numberofabspl(i)-1 do begin

         parametercount=parametercount+1
         indexpl=where(absploutbinwave(i,tt,*) ge lowerlimit and absploutbinwave(i,tt,*) le upperlimit)
         absfactorpl=absfactorpl*exp(-P(parametercount)*absploutbintau(i,tt,indexpl)) ;mathmathmath

     endfor
;---------


     parametercount=parametercount+1
     pl=(cc/(sourcewave(indsource)*10.d^(-6)))^P(parametercount) ;mathmathmath
     ; ^^compute powerlaw using frequency
     ; VV Normalization of powerlaw
     parametercount=parametercount+1
     sum=sum+P(parametercount)*absfactorpl*extfactorPL*pl/max(pl) ;mathmathmath
     cgoplot,sourcewave(indsource),P(parametercount)*absfactorpl*extfactorpl*pl/max(pl),color=colorc+i,linestyle=i
     PLmem(i,0:n_elements(pl)-1)=P(parametercount)*absfactorpl*extfactorpl*pl/max(pl) ; store the single pl's
     PLmemnorm(i)=double(P(parametercount)/max(pl))                                           ; store the normalizationfactor
  
     if keyword_set(ident) then begin 
         if P(parametercount) ne 0. then begin
             convert=strtrim(string(i+1),2) 
             maxi=max(sourcewave(indsource))
             xyouts,maxi,PLmem(i,n_elements(pl)-1),'PL'+convert
         endif
     endif

     parametercount=parametercount+1

endfor

;add+plot templates
;again for templates

for i=0,numberoftemplates-1 do begin

;check which extinctioncurve is used
     
    if P(parametercount) ne 0. then begin
        
        extmodelnumber=where(extcurvetemp(i) eq extcurvesynonym)
        extmodelnumberb=extmodelnumber(0) ; convert array(1) into integer ??? why is this nessesary
        
        if screenmixedtemplate(i) eq 'S' then begin
            
            extfactor=10.d^((-0.4d)*P(parametercount)*extoutbinvalue(extmodelnumberb,indextinction (extmodelnumberb,0:countext(extmodelnumberb)-1)))
            
        endif else begin
            taulambdatemp=P(parametercount)*0.4*alog(10)*extoutbinvalue(extmodelnumberb,indextinction (extmodelnumberb,0:countext(extmodelnumberb)-1)) ;repeats
            extfactor=(1-exp(-taulambdatemp))/taulambdatemp ;repeats
        endelse

    endif else begin
        extfactor=1.
    endelse

;----------------   
  
  ; estimate factor for absorption

     absfactortemp=1.

     for tt=0,numberofabstemp(i)-1 do begin

         parametercount=parametercount+1
         indextemp=where(abstempoutbinwave(i,tt,*) ge lowerlimit and abstempoutbinwave(i,tt,*) le upperlimit)
         absfactortemp=absfactortemp*exp(-P(parametercount)*abstempoutbintau(i,tt,indextemp)) ;mathmathmath

     endfor

;-------


     parametercount=parametercount+1
     sum=sum+P(parametercount)*absfactortemp*extfactor*outbinflux(i,indtemplate(i,0:count(i)-1))/max(outbinflux(i,indtemplate(i,0:count(i)-1)),/NaN) ;mathmathmath
     cgoplot,sourcewave(indsource),P(parametercount)*absfactortemp*extfactor*outbinflux(i,indtemplate(i,0:count(i)-1))/max(outbinflux(i,indtemplate(i,0:count(i)-1)),/NaN),color=colord[i],linestyle=i
     Tempmem(i,0:n_elements(sourcewave(indsource))-1)=P(parametercount)*absfactortemp*extfactor*outbinflux(i,indtemplate(i,0:count(i)-1))/max(outbinflux(i,indtemplate(i,0:count(i)-1)),/NaN); store the single temp's
     Tempmemnorm(i)=P(parametercount)/max(outbinflux(i,indtemplate(i,0:count(i)-1)),/NaN) ; store the normalizationfactor
        
     if keyword_set(ident) then begin 
         if P(parametercount) ne 0. then begin
             convert=strtrim(string(i+1),2) 
             maxi=max(sourcewave(indsource))
             xyouts,maxi,TEmpmem(i,n_elements(sourcewave(indsource))-1),'TL'+convert
         endif
     endif

     parametercount=parametercount+1
endfor

;------------------------------------------------------------------------------------------------------


;plot sum-spectra
;plot the source without absorption affecting it

cgoplot, sourcewave(indsource),sourceflux(indsource),thick=1.5
cgoplot, sourcewave(indsource),sum,color=colora,linestyle=3,thick=2
cgoplot, sourcewave(indsource),sum,color=colora,psym=3


; write sum-fit to file

if keyword_set(data) then begin

        datafile = sourcefiles(0)+'par:_'+controlfile+'_total.fit'
        openw,lun,pathout+datafile,/get_lun
        
        for j=0,n_elements(indsource)-1 do $
            printF,lun,sourcewave(indsource(j)),sum(j)
        free_lun,lun

endif

;prepare and set the legend
header='Type   Temperature (K)    Av Ext       Tau Peak'
BB1temp=string(result(2), format="(I0)")
BB1ext=string(result(0), format="(D14.2)")
BB1abs=string(result(1), format="(D14.2)")
BB2temp=string(result(6), format="(I0)")
BB2ext=string(result(4), format="(D14.2)")
BB2abs=string(result(5), format="(D14.2)")
BB3temp=string(result(10), format="(I0)")
BB3ext=string(result(8), format="(D15.2)")
BB3abs=string(result(9), format="(D15.2)")
PAHnorm1=string(result(18), format="(I0)")
PAHnorm2=string(result(21), format="(I0)")
BB1item='BB1          ' +BB1temp  + '  ' +BB1ext  + $
         BB1abs
BB2item='BB2          ' +BB2temp  + '  ' +BB2ext  + $
         BB2abs
BB3item='BB3          ' +BB3temp  + '  ' +BB3ext  + $
         BB3abs
if numberofBB eq 4 then begin
   BB4temp=string(result(14), format="(I0)")
   BB4ext=string(result(12), format="(D15.2)")
   BB4abs=string(result(13), format="(D15.2)")
   BB4item='BB4          ' +BB4temp  + '  ' +BB4ext  + $
      BB4abs
   al_legend, [header,BB1item,BB2item,BB3item,BB4item,$
               'PAH 1'+PAHnorm1,'PAH 2'+PAHnorm2],$
               colors=[cgcolor('red'),$
                       colorb[0],colorb[1],colorb[2],colorb[3],$
                       colord[0],colord[1]],$
               textcolors=[cgcolor('red'),$
                           colorb[0],colorb[1],colorb[2],colorb[3],$
                           colord[0],colord[1]], box=0
endif else begin
   al_legend, [header,BB1item,BB2item,BB3item,'PAH 1'+PAHnorm1,'PAH 2'+PAHnorm2], colors=[cgcolor('red'),colorb[0],colorb[1],colorb[2],colord[0],colord[1]], textcolors=[cgcolor('red'),colorb[0],colorb[1],colorb[2],colord[0],colord[1]], box=0
endelse

;plot residuum
if keyword_set(res) then begin
   diff=(sourceflux(indsource)/sum)
   if keyword_set(log) then begin
      plot, sourcewave(indsource), diff, pos=[0.2, 0.45, 0.95, 0.65], xstyle=1,ystyle=1,xtitle='rest-wavelength [�m]',ytitle='Data/Model', xrange=[lowerlimit, upperlimit], yrange=[0., 1.999],/xlog;!xtickinterval=0.001 D.C.
   endif else begin
      plot, sourcewave(indsource), diff, pos=[0.2, 0.45, 0.95, 0.65], xstyle=1,ystyle=1,xtitle='rest-wavelength [�m]',ytitle='Data/Model', xrange=[lowerlimit, upperlimit], yrange=[0., 1.999];!xtickinterval=0.001 D.C.
   endelse

   cgoplot,sourcewave(indsource),1.*(sourcewave(indsource)/sourcewave(indsource)),linestyle=3 ;D.C
endif




;find the fitted parameters in fit result
;---------------------------------------

;create arrays for normalization factors, indecies, extinction fits
;and tau fits for BBs, PLs and templates
normBBfit=dblarr(numberofBB)
temperatureBBfit=dblarr(numberofBB)
extinctionBBfit=dblarr(numberofBB)
tauBBfit=dblarr(numberofBB,max(numberofabsbb))

if numberofpl gt 0 then begin
  normPLfit=dblarr(numberofpl)
  indexPLfit=dblarr(numberofpl)
  extinctionPLfit=dblarr(numberofpl)
  tauPLfit=dblarr(numberofpl,max(numberofabspl))
endif

normtempfit=dblarr(numberoftemplates)
extinctiontempfit=dblarr(numberoftemplates)
tautempfit=dblarr(numberoftemplates,max(numberofabstemp))

BBfittemps=dblarr(numberofBB)
BBtauext=dblarr(numberofBB)
BBtauabs=dblarr(numberofBB)


;BB 
;fill the array using MPfit data
c=0
for i=0,numberofBB-1 do begin               
    extinctionBBfit(i)=result(c)  ;sets extinction parameters calculated using MPFit to the appropriate section of the array
    convertc=strtrim(string(c),2)
    convertii=strtrim(string(i+1),2)
    print,'P('+convertc+'):  BB :'+convertii+' fitted extinction (Av): ',result(c)
    c=c+1
    for ii=0,numberofabsbb(i)-1 do begin
        convertc=strtrim(string(c),2)
        tauBBfit(i,ii)=result(c) ;sets tau peak for BB calculated in MPFit
        print,'P('+convertc+'):  BB :'+convertii+' fitted taupeak: ',result(c)
        c=c+1 
    endfor
    temperatureBBfit(i)=result(c) ;sets temperature fits caluclated in MPFit
    convertc=strtrim(string(c),2)
    
    print,'P('+convertc+'):  BB :'+convertii+' fitted temperature: ',result(c)
    c=c+1
    normBBfit(i)=result(c)   ;sets normalization fit for BB calculated in MPFit
    convertc=strtrim(string(c),2)

    print,'P('+convertc+'):  BB :'+convertii+' fitted normfactor: (/max) ',result(c)
    print,'[P('+convertc+')]:  BB :'+convertii+' fitted total normalization factor:',BBmemnorm(i)
    print,' '  
    c=c+1    
endfor


print,'------------------- '
print,' '
;PL
;repeat process for PLs

for i=0,numberofpl-1 do begin               
    extinctionPLfit(i)=result(c) ;sets extinction for power laws to appropriate position in array

    convertii=strtrim(string(i+1),2)
    convertc=strtrim(string(c),2)

    print,'P('+convertc+'):  PL :'+convertii+' fitted extinction (Av): ',result(c)
    c=c+1
    for ii=0,numberofabsPL(i)-1 do begin
        convertc=strtrim(string(c),2)
        tauPLfit(i,ii)=result(c) ;repeat tau peak assignment for power law
        print,'P('+convertc+'):  PL :'+convertii+' fitted taupeak: ',result(c)
        c=c+1 
    endfor
    indexPLfit(i)=result(c)
    convertc=strtrim(string(c),2)
    print,'P('+convertc+'):  PL :'+convertii+' fitted powerlaw-index: ',result(c)
    c=c+1
    normPLfit(i)=result(c)   ;repeat normalization assignment for power law
    convertc=strtrim(string(c),2)
    print,'P('+convertc+'):  PL :'+convertii+' fitted normfactor: (/max) ',result(c)
    print,'[P('+convertc+')]:  PL :'+convertii+' fitted total normalization factor:',PLmemnorm(i)
    print,' '  
    c=c+1    
endfor

print,'------------------- '
print,' '

;Template 
;repeat process for templates
for i=0,numberoftemplates-1 do begin
  
      extinctiontempfit(i)=result(c) ;repeat ext assign for templates
      convertii=strtrim(string(i+1),2)
      convertc=strtrim(string(c),2)
      print,'file: ',templatefilesb(i)
      print,'P('+convertc+'):  template :'+convertii+' fitted extinction (Av): ',result(c)
      c=c+1
      for ii=0,numberofabstemp(i)-1 do begin
          convertc=strtrim(string(c),2)
          tautempfit=result(c) ;tau peak assign for temp
          print,'P('+convertc+'):  template :'+convertii+' fitted taupeak: ',result(c)
          c=c+1 
      endfor
      normtempfit(i)=result(c) ;normalization assign for temp
      convertc=strtrim(string(c),2)
      print,'P('+convertc+'):  template :'+convertii+' fitted normfactor: (/max) ',result(c)
      print,'[P('+convertc+')]:  template :'+convertii+' fitted total normalization factor:',tempmemnorm(i)
      print,' '  
      c=c+1

endfor


print,'------------------- '
print,' '

; if ps=1 write fitresults to file

if keyword_set(ps) then begin
    
    parfile = ID+'_'+cluster+'.dat'
    openw,parlun,pathout+parfile,/get_lun
    
;BB 
    c=0
    for i=0,numberofBB-1 do begin               
        extinctionBBfit(i)=result(c)
        convertc=strtrim(string(c),2)
        convertii=strtrim(string(i+1),2)
        printf,parlun,'P('+convertc+'):  BB :'+convertii+' fitted extinction (Av): ',result(c)
        c=c+1
        for ii=0,numberofabsbb(i)-1 do begin
            convertc=strtrim(string(c),2)
            tauBBfit(i,ii)=result(c)
            printf,parlun,'P('+convertc+'):  BB :'+convertii+' fitted taupeak: ',result(c)
            c=c+1 
        endfor
        temperatureBBfit(i)=result(c)
        convertc=strtrim(string(c),2)
        
        printf,parlun,'P('+convertc+'):  BB :'+convertii+' fitted temperature: ',result(c)
        c=c+1
        normBBfit(i)=result(c)   
        convertc=strtrim(string(c),2)
        
        printf,parlun,'P('+convertc+'):  BB :'+convertii+' fitted normfactor: (/max) ',result(c)
        printf,parlun,'[P('+convertc+')]:  BB :'+convertii+' fitted total normalization factor:',BBmemnorm(i)
        printf,parlun,' '  
        c=c+1    
    endfor
    
    
    printf,parlun,'------------------- '
    printf,parlun,' '
;PL
    
    for i=0,numberofpl-1 do begin               
        extinctionPLfit(i)=result(c)
        
        convertii=strtrim(string(i+1),2)
        convertc=strtrim(string(c),2)
        
        printf,parlun,'P('+convertc+'):  PL :'+convertii+' fitted extinction (Av): ',result(c)
        c=c+1
        for ii=0,numberofabsPL(i)-1 do begin
            convertc=strtrim(string(c),2)
            tauPLfit(i,ii)=result(c)
            printf,parlun,'P('+convertc+'):  PL :'+convertii+' fitted taupeak: ',result(c)
            c=c+1 
        endfor
        indexPLfit(i)=result(c)
        convertc=strtrim(string(c),2)
        printf,parlun,'P('+convertc+'):  PL :'+convertii+' fitted powerlaw-index: ',result(c)
        c=c+1
        normPLfit(i)=result(c)   
        convertc=strtrim(string(c),2)
        printf,parlun,'P('+convertc+'):  PL :'+convertii+' fitted normfactor: (/max) ',result(c)
        printf,parlun,'[P('+convertc+')]:  PL :'+convertii+' fitted total normalization factor:',PLmemnorm(i)
        printf,parlun,' '  
        c=c+1    
    endfor
    
    printf,parlun,'------------------- '
    printf,parlun,' '
    
;Template 
    for i=0,numberoftemplates-1 do begin
        extinctiontempfit(i)=result(c)
        convertii=strtrim(string(i+1),2)
        convertc=strtrim(string(c),2)
        printf,parlun,'file :',templatefilesb(i)
        printf,parlun,'P('+convertc+'):  template :'+convertii+' fitted extinction (Av): ',result(c)
        c=c+1
        for ii=0,numberofabstemp(i)-1 do begin
            convertc=strtrim(string(c),2)
            tautempfit=result(c)
            printf,parlun,'P('+convertc+'):  template :'+convertii+' fitted taupeak: ',result(c)
            c=c+1 
        endfor
        normtempfit(i)=result(c)
        convertc=strtrim(string(c),2)
        printf,parlun,'P('+convertc+'):  template :'+convertii+' fitted normfactor: (/max) ',result(c)
        printf,parlun,'[P('+convertc+')]:  template :'+convertii+' fitted total normalization factor:',tempmemnorm(i)
        printf,parlun,' '  
        c=c+1
    endfor
    
endif




















;define frequency-axis for integration of Fnu over nu

nu=double(cc/(sourcewave(indsource)*10.d^(-6.d)))




; following section plots percentage portion of fit
; estimate the contribution of the single models in % depending on wavelength
; + integrate the flux under each contributor
; + write the data to files
;------------------------------------------------------------------------------------------------------

contribBB=dblarr(numberofBB,n_elements(sourcewave(indsource)))
intBB=dblarr(numberofBB)

nuxx=nu
sortnu=nu(sort(nuxx))
sortind=indsource(uniq(sort(nuxx)))
sortnu1=sortnu[uniq(sortnu)]
; VV if statement added below to skip 'integral contribution' of
; specific components because code does not work if there are
; duplicate points, thus it must be skipped to allow for data with
; duplicate points
if keyword_set(contrib) then begin
   if n_elements(sortnu) eq n_elements(sortnu1) then begin
      intsource=int_tabulated(sortnu,sourceflux(sortind),/double)
      intfluxsource=intsource*10.d^(-30) 
      print,'integral [W/cm^2] for source:',intfluxsource
      
;ERROR CORRECTION:
;plot,sourcewave(indsource),contribBB(i,*)/contribBB(i,*)*100,yrange=[0, 100], pos=[0.2,0.025,0.95,0.37],xstyle=1,ystyle=1,ytitle='[ % ]',title='contribution in %'
;(this was first written within the for next loop) => erasing of plot if more
;then 2 BBs !!!???
      
      
      
      if keyword_set(log) then begin
         plot,sourcewave(indsource),sourceflux(indsource)/sourceflux(indsource)*100,xrange=[lowerlimit, upperlimit],yrange=[0, 100], pos=[0.2,0.025,0.95,0.37],xstyle=1,ystyle=1,ytitle='[ % ]',title='contribution in %',/xlog  
      endif else begin
         plot,sourcewave(indsource),sourceflux(indsource)/sourceflux(indsource)*100,xrange=[lowerlimit, upperlimit],yrange=[0, 100], pos=[0.2,0.025,0.95,0.37],xstyle=1,ystyle=1,ytitle='[ % ]',title='contribution in %'
      endelse
      
      
      for i=0,numberofBB-1 do begin
         
         contribBB(i,*)=(BBmem(i,*)/sum(*))*100
         cgoplot,sourcewave(indsource),contribBB(i,*),color=colorb[i],linestyle=i
         
;sorting of frequencies and spectra for integration
         
         sortBB=BBmem(i,sort(nuxx))
         
         intBB(i)=int_tabulated(sortnu,sortBB,/double)
         
         if keyword_set(ps) then begin
            printf,parlun,'integral [W/cm^2] of BB no.:',i+1,' :',intBB(i)*10.d^(-30) 
            printf,parlun,'contribution [%] of BB no.:',i+1,' to integrated sourceflux :',(intBB(i)*10.d^(-30)/intfluxsource)*100
         endif
         
         print,'integral [W/cm^2] of BB no.:',i+1,' :',intBB(i)*10.d^(-30) 
         print,'contribution [%] of BB no.:',i+1,' to integrated sourceflux :',(intBB(i)*10.d^(-30)/intfluxsource)*100     
;file properties/open files and processing
         if keyword_set(data) and normBBfit(i) ne 0. then begin
            
            convert=strtrim(string(i+1),2)
            
            openw,1,pathout+sourcefiles(0)+'_par:'+controlfile+'_BB'+convert+'.fit'
            
            for j=0,n_elements(indsource)-1 do begin 
               printF,1,sourcewave(indsource(j)),BBmem(i,j)
            endfor
            close,1    
            
         endif
         
      endfor
      
      
      contribpl=dblarr(numberofpl,n_elements(sourcewave(indsource)))
      if numberofpl gt 0 then intPL=dblarr(numberofPL) else intPL=0d
      
      for i=0,numberofpl-1 do begin
         
         contribpl(i,*)=(plmem(i,*)/sum(*))*100
         
         cgoplot,sourcewave(indsource),contribpl(i,*),linestyle=i,color=colorc+i
         
         
         sortPl=PLmem(i,sort(nuxx))
         
         intPL(i)=int_tabulated(sortnu,sortPl,/double)
         
         if keyword_set(ps) then begin
            printf,parlun,'integral [W/cm^2] of PL no.:',i+1,' :',intPL(i)*10.d^(-30) 
            printf,parlun,'contribution [%] of PL no.:',i+1,' to integrated sourceflux :',(intPL(i)*10.d^(-30)/intfluxsource)*100
         endif
         
         print,'integral [W/cm^2] of PL no.:',i+1,' :',intPL(i)*10.d^(-30) 
         print,'contribution [%] of PL no.:',i+1,' to integrated sourceflux :',(intPL(i)*10.d^(-30)/intfluxsource)*100
         
         if keyword_set(data) and normPLfit(i) ne 0. then begin
            
            convert=strtrim(string(i+1),2)
            openw,1,pathout+sourcefiles(0)+'_par:'+controlfile+'_PL'+convert+'.fit'
            
            for j=0,n_elements(indsource)-1 do begin 
               printF,1,sourcewave(indsource(j)),PLmem(i,j)
            endfor
            close,1    
            
         endif
         
         
      endfor
      
      contribtemp=dblarr(numberoftemplates,n_elements(sourcewave(indsource)))
      inttemp=dblarr(numberoftemplates)
      
      for i=0,numberoftemplates-1 do begin
         
         contribtemp(i,*)=(tempmem(i,*)/sum(*))*100
         
         cgoplot,sourcewave(indsource),contribtemp(i,*),linestyle=i,color=colord[i]
         
         
         sorttemp=tempmem(i,sort(nuxx))
         
         inttemp(i)=int_tabulated(sortnu,sorttemp,/double)
         
         if keyword_set(ps) then begin
            printf,parlun,'integral [W/cm^2] of template no.:',i+1,' :',inttemp(i)*10.d^(-30) 
            printf,parlun,'contribution temp no.:',i+1,' to int:',(inttemp(i)*10.d^(-30)/intfluxsource)*100
         endif
         
         print,'integral [W/cm^2] of template no.:',i+1,' :',inttemp(i)*10.d^(-30) 
         print,'contribution [%] of template no.:',i+1,' to integrated sourceflux :',(inttemp(i)*10.d^(-30)/intfluxsource)*100
         if keyword_set(data) and normtempfit(i) ne 0. then begin
            
            convert=strtrim(string(i+1),2)
            openw,1,pathout+sourcefiles(0)+'_par:'+controlfile+'_template'+convert+'.fit'
            
            for j=0,n_elements(indsource)-1 do begin 
               printF,1,sourcewave(indsource(j)),tempmem(i,j)
            endfor
            close,1    
            
         endif
         
      endfor
      
      
      if keyword_set(ps) then begin
         printf,parlun,'Sum of all contributions [%] relative to int. sourceflux :',(total(inttemp)*10.d^(-30)/intfluxsource)*100+(total(intPL)*10.d^(-30)/intfluxsource)*100+(total(intBB)*10.d^(-30)/intfluxsource)*100
         printf,parlun,' '
         printf,parlun,'Chi^2 = ',total(((sourceflux(indsource)-sum)/sourcestdev(indsource))^2)
         free_lun,parlun
      endif
      
      
      print,'Sum of all contributions [%] relative to int. sourceflux :',(total(inttemp)*10.d^(-30)/intfluxsource)*100+(total(intPL)*10.d^(-30)/intfluxsource)*100+(total(intBB)*10.d^(-30)/intfluxsource)*100
      print,' '
   endif
endif



if keyword_set(ps) then begin
   device,/close
   set_plot,'x'
endif

savefile = ID+'_'+cluster+'.sav'
save,/VARIABLES,filename=pathout+savefile

print,'Output files found in dir. '+pathout
print,''
print,'Writing all program variables to '+savefile
if keyword_set(data) then print,'Writing model spectrum to '+datafile
if keyword_set(ps) then begin
   print,'Writing best fit parameters to '+parfile
   print,'Plotting best fit in '+epsfile
endif


!P.multi=[0,1,1]
close,/all
end
