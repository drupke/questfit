; docformat = 'rst'
;
;+
;
; Purpose: Compute PAH and continuum fluxes from spectrum fits for
; IDEOS project.
;
; Inputs:
;   1. Parameters of template spectra (outputs from PAHFIT)
;   2. PAH peak intensities for each fit
;   3. blackbody temperatures, peak intensities, extinctions for each
;      full-spectrum fit
;   4. extinction curves used in fits (silicate and ice)
;
; Outputs:
;   1. PAH feature fluxes, in log(10^-26 W/m^2)
;   2. continuum flux densities, in log(10^-26 W/m^-2/um), extincted
;      and unextincted
;
; Procedure:
;   1. Import necessary inputs: see above.
;   2. Compute PAH feature intensities in templates at the following
;      wavelengths: 6.2, 7.7, 8.6, 11.3, 12.7, 17 um.  Compute total
;      PAH spectrum fluxes for each template.
;   3. Interpolate extinction curves to binning of working total PAH
;      and continuum spectra. Use binning from pure PAH template
;      files.
;   4. For each galaxy, compute PAH spectra, continuum spectra, and
;      PAH feature fluxes, extincted and unextincted in each case.
;      a. Get parameters from fit.
;      b. Compute template spectra.
;      c. Compute PAH feature fluxes from PAHFIT output.
;      d. Find wavelength of peak PAH emission (in Jy).
;      e. Using template (total spectra and individual feature PAH
;         fluxes) and actual fitted fluxes, compute actual feature PAH
;         fluxes and blackbody flux densities (a simple scaling). PAH
;         template feature strengths units = W/m^2 for a peak flux of
;         1 MJy as output by PAHFIT.
;      f. Compute relative contribution of each template to the total
;         flux.
;   5. Write outputs to a single file.
;
; :Categories:
;    IDEOS
;
; :Returns:
;    Data file with PAH fluxes.
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
;      2014jan15, DSNR, copied from quest_pah.pro for IDEOS project
;      2015mar18, DSNR, modified to output PAH equivalent widths
;    
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
pro ideospah

;========
; INPUTS
;========

; PAH templates
  readcol300,'smith_nftemp3.pahpars.ext.dat',$
             pp1wl,pp1fwhm,pp1peak,pp1flux,/silent
  readcol300,'smith_nftemp4.pahpars.next.dat',$
             pp2wl,pp2fwhm,pp2peak,pp2flux,/silent
; Extinction curves
  readcol300,'extinction_chiar06_original.dat',ec0wl,ec0arat,/silent,skipline=1
  readcol300,'extinction_ice+hc.dat',icewl,icetau,/silent,skipline=1
; normalize ice tau profile to 1 at peak
  icetau=icetau/max(icetau)
; Galaxy fit results tables
  readcol300,'fit_results_table.dat',id1,$
             bbf1av,bbf1tau,bbf1temp,bbf1norm,$
             bbf2av,bbf2tau,bbf2temp,bbf2norm,$
             bbf3av,bbf3tau,bbf3temp,bbf3norm,$
             pf1av,pf1tau,pf1norm,$
             pf2av,pf2tau,pf2norm,$
             fitlo,$
             fithi,$
             format='A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,',/silent

;=======
; MAIN
;=======

; Common wavelength scale
  wl_i=1
  wl_f=40
  wl_delt=0.01
  wlcommon=dblarr((wl_f-wl_i)/wl_delt)
  wlcommon[0]=wl_i
  for i=1,n_elements(wlcommon)-1 do $
     wlcommon[i]=wlcommon[i-1]+wl_delt

; Unextincted PAH feature fluxes in two templates, in W/m^2 (with
; normalization from noise-free SINGS template fit)
  p1f6=pp1flux[2]
  p1f7=pp1flux[4]+pp1flux[5]+pp1flux[6]
  p1f8=pp1flux[8]
  p1f11=pp1flux[10]+pp1flux[11]
  p1f12=pp1flux[13]+pp1flux[14]
  p1f17=pp1flux[19]+pp1flux[20]+pp1flux[21]+pp1flux[22]
  p2f6=pp2flux[2]
  p2f7=pp2flux[4]+pp2flux[5]+pp2flux[6]
  p2f8=pp2flux[8]
  p2f11=pp2flux[10]+pp2flux[11]
  p2f12=pp2flux[13]+pp2flux[14]
  p2f17=pp2flux[19]+pp2flux[20]+pp2flux[21]+pp2flux[22]

; Total unextincted PAH fluxes in templates, in W/m^2 (with
; normalization from noise-free SINGS template fit)
  p1ftot=0
  p2ftot=0
  for i=0,n_elements(pp1flux)-1 do begin
     p1ftot=p1ftot+pp1flux[i]
     p2ftot=p2ftot+pp2flux[i]
  endfor

; Interpolate extinction curves
  ecarat_int=interpol(ec0arat,ec0wl,wlcommon)
  icetau_int=interpol(icetau,icewl,wlcommon)

; Open output file
  openw,outlun1,'ideos_pahfluxes.txt',/get_lun
  openw,outlun2,'ideos_pahweq.txt',/get_lun

;
; Cycle through all fits
;

  ngal=n_elements(id1)          ; number of galaxies

; BB temp and AV arrays
  bbtemp=dblarr(3)
  bbav=dblarr(3)
;
  for i=0,ngal-1 do begin

     print,id1[i]

;
;    Create empty scalars/arrays for synthetic PAH and blackbody fit
;    components
;
;    unextincted full PAH spectrum
     ps1flux     =dblarr(n_elements(wlcommon))
     ps2flux     =dblarr(n_elements(wlcommon))
;    blackbodies
     bbs1flux    =dblarr(n_elements(wlcommon))
     bbs2flux    =dblarr(n_elements(wlcommon))
     bbs3flux    =dblarr(n_elements(wlcommon))
     bbs1flux_ext=dblarr(n_elements(wlcommon))
     bbs2flux_ext=dblarr(n_elements(wlcommon))
     bbs3flux_ext=dblarr(n_elements(wlcommon))
     bbs1flux_ext_noice=dblarr(n_elements(wlcommon))
     bbs2flux_ext_noice=dblarr(n_elements(wlcommon))
     bbs3flux_ext_noice=dblarr(n_elements(wlcommon))
;    water + hydrocarbon feature equivalent width
     h2ohc_weq=0.

;    Cycle through wavelengths
     for j=0,n_elements(wlcommon)-1 do begin

;       This factor multiplies by c/lambda^2*10^-20*delta_lambda to
;       convert from flux densities to fluxes (the former in MJy, the
;       latter in W/m^2)
        fac = 2.99792/wlcommon[j]^2/1d6*wl_delt


;       Cycle through individual PAH features
        for k=0,n_elements(pp1wl)-1 do begin
;          Unextincted PAH spectra of individual features, normalized
;          to SINGS template fit
           pp1flux=pp1peak[k]*pp1fwhm[k]^2 / $
                   ((wlcommon[j]/pp1wl[k] - pp1wl[k]/wlcommon[j])^2 + $
                    pp1fwhm[k]^2)
           pp2flux=pp2peak[k]*pp2fwhm[k]^2 / $
                   ((wlcommon[j]/pp2wl[k] - pp2wl[k]/wlcommon[j])^2 + $
                    pp2fwhm[k]^2)
;          Unextincted total PAH spectra, normalized to SINGS template
;          fit
           ps1flux[j]=ps1flux[j]+pp1flux
           ps2flux[j]=ps2flux[j]+pp2flux
        endfor

;       Unextincted and extincted blackbody spectra, arbitrary normalization
        bbs1flux[j]=1/wlcommon[j]^3/(exp(14387.87/wlcommon[j]/bbf1temp[i])-1)
        bbs2flux[j]=1/wlcommon[j]^3/(exp(14387.87/wlcommon[j]/bbf2temp[i])-1)
        bbs3flux[j]=1/wlcommon[j]^3/(exp(14387.87/wlcommon[j]/bbf3temp[i])-1)
        bbs1flux_ext_noice[j]=bbs1flux[j]*exp(-bbf1av[i]*ecarat_int[j]/1.086)
        bbs2flux_ext_noice[j]=bbs2flux[j]*exp(-bbf2av[i]*ecarat_int[j]/1.086)
        bbs3flux_ext_noice[j]=bbs3flux[j]*exp(-bbf3av[i]*ecarat_int[j]/1.086)
        bbs1flux_ext[j]=bbs1flux_ext_noice[j]*exp(-bbf1tau[i]*icetau_int[j])
        bbs2flux_ext[j]=bbs2flux_ext_noice[j]*exp(-bbf2tau[i]*icetau_int[j])
        bbs3flux_ext[j]=bbs3flux_ext_noice[j]*exp(-bbf3tau[i]*icetau_int[j])

     endfor

;    Maxima in unextincted spectra
     ps1flux_max=max(ps1flux,ps1ind_max)
     ps2flux_max=max(ps2flux,ps2ind_max)
     bbs1flux_max=max(bbs1flux,bbs1ind_max)
     bbs2flux_max=max(bbs2flux,bbs2ind_max)
     bbs3flux_max=max(bbs3flux,bbs3ind_max)

;    If black-body maxima occur at fit boundaries, adjust indices
;    accordingly.
     for ctr=0,n_elements(fitlo) do begin
        fitloind=value_locate(wlcommon,fitlo[i]) ;for loop this section
        fithiind=value_locate(wlcommon,fithi[i])
     endfor
        if (bbs1ind_max lt fitloind) then bbs1ind_max=fitloind
        if (bbs1ind_max gt fithiind) then bbs1ind_max=fithiind
        if (bbs2ind_max lt fitloind) then bbs2ind_max=fitloind
        if (bbs2ind_max gt fithiind) then bbs2ind_max=fithiind
        if (bbs3ind_max lt fitloind) then bbs3ind_max=fitloind
        if (bbs3ind_max gt fithiind) then bbs3ind_max=fithiind
        bbs1flux_max=bbs1flux[bbs1ind_max]
        bbs2flux_max=bbs2flux[bbs2ind_max]
        bbs3flux_max=bbs3flux[bbs3ind_max]

;    Locations of desired blackbody fluxes
     w5ind=value_locate(wlcommon,5.)
     w6ind=value_locate(wlcommon,6.22)
     w7ind=value_locate(wlcommon,7.7)
     w8ind=value_locate(wlcommon,7.90)
     w11ind=value_locate(wlcommon,11.2)
     w15ind=value_locate(wlcommon,15.00)
     w25ind=value_locate(wlcommon,25.00)
     w30ind=value_locate(wlcommon,30.00)
     wsilind=value_locate(wlcommon,9.7) ; not sure about this wavelength exactly

;
; Scale PAH feature and continuum fluxes from template normalization
; to real flux scale.  Also compute fractional contribution to total
; PAH flux from each template.
;
; Scaling:
;
;  F_real(lambda)=F_template(lambda)*[F_real(lambda_max)/F_template(lambda_max)]
;                =F_template(lambda)*fluxrat_max
;
; Scaling is independent of extinction, since fluxrat_max is constant
; if F_real and F_template at lambda_max are extincted simultaneously.
; Scaling applies to fluxes as well as flux densities, since fluxes
; are just integrals of flux densities.
;
; Factor 10^6 converts template normalization from MJy to Jy.
; Factor 10^26 converts template fluxes from W/m^2 to 10^-26 W/m^2.
; Factor c/lambda^2 converts from 10^-26 W/m^2/Hz to 10^-26 W/m^2/um.
;

;
;    Properly fluxed blackbody spectra
;
;    Normalizations
     bb1fluxrat_max=bbf1norm[i]/bbs1flux_max
     bb2fluxrat_max=bbf2norm[i]/bbs2flux_max
     bb3fluxrat_max=bbf3norm[i]/bbs3flux_max
;    Spectra
     bb1_ext_noice = bbs1flux_ext_noice*bb1fluxrat_max
     bb2_ext_noice = bbs2flux_ext_noice*bb2fluxrat_max
     bb3_ext_noice = bbs3flux_ext_noice*bb3fluxrat_max
     bb1_ext = bbs1flux_ext*bb1fluxrat_max
     bb2_ext = bbs2flux_ext*bb2fluxrat_max
     bb3_ext = bbs3flux_ext*bb3fluxrat_max
     bb1 = bbs1flux*bb1fluxrat_max
     bb2 = bbs2flux*bb2fluxrat_max
     bb3 = bbs3flux*bb3fluxrat_max
     bb_ext_noice = bb1_ext_noice + bb2_ext_noice + bb3_ext_noice
     bb_ext = bb1_ext + bb2_ext + bb3_ext
     bb = bb1 + bb2 + bb3

;
;    Properly fluxed PAH spectra
;
;    Normalizations
     p1fluxrat_max = pf1norm[i]/ps1flux_max
     p2fluxrat_max = pf2norm[i]/ps2flux_max
;    Fluxes
     if (pf1norm[i] eq 0. AND pf2norm[i] eq 0.) then begin
        p6  =-1000.
        p7  =-1000.
        p8  =-1000.
        p11 =-1000.
        p12 =-1000.
        p17 =-1000.
        ptot=-1000.
        p1frac=0.
        p2frac=0.
        p6_ext=-1000.
        p7_ext=-1000.
        p8_ext=-1000.
        p11_ext=-1000.
        p12_ext=-1000.
        p17_ext=-1000.
        ptot_ext=-1000
        p1frac_ext=0.
        p2frac_ext=0.
     endif else begin
        p6=alog10(p1fluxrat_max*p1f6  +p2fluxrat_max*p2f6 )-6.+26.
        p7=alog10(p1fluxrat_max*p1f7  +p2fluxrat_max*p2f7 )-6.+26.
        p8=alog10(p1fluxrat_max*p1f8  +p2fluxrat_max*p2f8 )-6.+26.
        p11=alog10(p1fluxrat_max*p1f11 +p2fluxrat_max*p2f11)-6.+26.
        p12=alog10(p1fluxrat_max*p1f12 +p2fluxrat_max*p2f12)-6.+26.
        p17=alog10(p1fluxrat_max*p1f17 +p2fluxrat_max*p2f17)-6.+26.
        ptot=alog10(p1fluxrat_max*p1ftot+p2fluxrat_max*p2ftot )-6.+26.
        p1frac=p1fluxrat_max*p1ftot / (p1fluxrat_max*p1ftot+p2fluxrat_max*p2ftot)
        p2frac=p2fluxrat_max*p2ftot / (p1fluxrat_max*p1ftot+p2fluxrat_max*p2ftot)
;
     endelse
;    Spectra
     p1 = ps1flux*p1fluxrat_max
     p2 = ps2flux*p2fluxrat_max
        
;
;    Properly fluxed total spectra
;
     tot_ext = p1+p2+bb_ext
     tot = p1+p2+bb

;    Blackbody flux densities (10^-26 W/m^2 per micron)
     bb6 =alog10(bb[w6ind])+alog10(2.99792)+14-2*alog10(6.22)
     bb8 =alog10(bb[w8ind])+alog10(2.99792)+14-2*alog10(7.90)
     bb11 =alog10(bb[w11ind])+alog10(2.99792)+14-2*alog10(11.2)
     bb15 =alog10(bb[w15ind])+alog10(2.99792)+14-2*alog10(15.00)
     bb30 =alog10(bb[w30ind])+alog10(2.99792)+14-2*alog10(30.00)
     bb6_ext =alog10(bb_ext[w6ind])+alog10(2.99792)+14-2*alog10(6.22)
     bb8_ext =alog10(bb_ext[w8ind])+alog10(2.99792)+14-2*alog10(7.90)
     bb11_ext =alog10(bb_ext[w11ind])+alog10(2.99792)+14-2*alog10(11.2)
     bb15_ext =alog10(bb_ext[w15ind])+alog10(2.99792)+14-2*alog10(15.00)
     bb30_ext =alog10(bb_ext[w30ind])+alog10(2.99792)+14-2*alog10(30.00)
     bb6_ext_noice =alog10(bb_ext_noice[w6ind])+alog10(2.99792)+14-2*alog10(6.22)
;    Integrated optical depth at 9.7 um
     bbsiltau=-alog(bb_ext[wsilind]/bb[wsilind])
;    Maximum value of A_V
     bbavmax=max([bbf1av[i],bbf2av[i],bbf3av[i]])

;    mid-IR blackbody fluxes, in units of W/m^2
     bbtot=0
     bbtot_ext=0
;    convert to flux densities per micron
     for j=w5ind,w25ind do begin
        bbtot=bbtot+(bb[j])/wlcommon[j]/wlcommon[j]
        bbtot_ext=bbtot_ext+(bb_ext[j])/wlcommon[j]/wlcommon[j]
     endfor
     bbtot=alog10(bbtot)+alog10(2.99792)+14+alog10(wl_delt)
     bbtot_ext=alog10(bbtot_ext)+alog10(2.99792)+14+alog10(wl_delt)

;
;    Sort BB temperatures and extinctions
;
;    Set temp to 0 if component zeroed out
     if (bbf1norm[i] eq 0) then begin
        bbtemp[0] = bbf2temp[i]
        bbtemp[1] = bbf3temp[i]
        bbtemp[2] = 0
        bbav[0] = bbf2av[i]
        bbav[1] = bbf3av[i]
        bbav[2] = -1
     endif else begin
        bbtemp[0] = bbf1temp[i]
        bbtemp[1] = bbf2temp[i]
        bbtemp[2] = bbf3temp[i]
        bbav[0] = bbf1av[i]
        bbav[1] = bbf2av[i]
        bbav[2] = bbf3av[i]
     endelse
     bbtemp_i = sort(bbtemp)
     bbtemp = bbtemp[bbtemp_i]
     bbav = bbav[bbtemp_i]
;
;    H2O + hydrocarbon equivalent width
;
;    Integrate over wavelengths
     for j=0,n_elements(wlcommon)-1 do $
        h2ohc_weq = $
        h2ohc_weq + ((bb_ext_noice[j] - bb_ext[j])/bb_ext[j])*wl_delt

;    PAH equivalent widths [observed frame]
     p6weq = p6 - bb6
     p7weq = p7 - bb8
     p11weq = p11 - bb11
     
;    Print to output file
     printf,outlun1,id1[i],p6,p7,p8,p11,p12,p17,$
            ptot,p1frac,p2frac,pf1av[i],pf2av[i],$
            format='(A-30,11D-14.6)'

     printf,outlun2,id1[i],p6weq,p7weq,p11weq,format='(A-30,3D-14.6)'

  endfor
  
  free_lun,outlun1
  free_lun,outlun2
  
end
