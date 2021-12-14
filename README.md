## OVERVIEW

This software package was originally used to fit the Spitzer
mid-infrared spectra of the QUEST (**Q**uasar **U**LIRG and
**E**volution **ST**udy) sample, as described in [Schweitzer et
al. 2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...649...79S/abstract),
[Schweitzer et
al. 2008](https://ui.adsabs.harvard.edu/abs/2008ApJ...679..101S/abstract),
and [Veilleux et
al. 2009](https://ui.adsabs.harvard.edu/abs/2009ApJS..182..628V/abstract). It
uses two PAH templates from [Smith et
al. 2007](https://ui.adsabs.harvard.edu/abs/2007ApJ...656..770S/abstract)
atop an extincted and absorbed continuum model to fit the mid-IR
spectra of galaxies that are heavily-absorbed and AGN with silicate
emission.

The current version of `QUESTFIT` is optimized for processing spectra
from the CASSIS (**C**ombined **A**tlas of **S**ources with
**S**pitzer **I**RS **S**pectra)
[portal](https://cassis.sirtf.com/atlas/welcome.shtml) to produce PAH
fluxes for heavily absorbed sources. This method is described in Spoon
et al. 2021 (submitted). These PAH fluxes will appear in the IDEOS
(**I**nfrared **D**atabase of **E**xtragalactic **O**bservables from
**S**pitzer) [portal](http://ideos.astro.cornell.edu/).

## REQUIREMENTS

IDL (tested with v8.5; may work with pre-v8.0 versions)

IDL libraries:
- [IDL Astronomy User's Library](http://idlastro.gsfc.nasa.gov)
- [Coyote](http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD), for graphics AND undefine.pro
  - or from the [GitHub repository](https://github.com/davidwfanning/idl-coyote/tree/master/coyote)
- [ISAP](https://old.ipac.caltech.edu/iso/isap/isap.html)

Notes:
- The Spitzer IRS data reduction software SMART ships with a slightly
modified version of the ISAP library. At least one instance of errors
has been documented using this version of the library.
- The IDL Astronomy User's Library ships with some Coyote
routines. However, it's not clear how well these libraries keep track
of each other, so it may be preferable to download each package
separately and delete the redundant routines that ship within other
packages.

## QUICK TEST

See files in the `questfit/test` subdirectory:

IDEOS ID 4978688_0 = Mrk 231

IDL> questfit, control file='4978688_0.ideos.cf',pathin='[path-to-questfit]/questfit/test/',pathout='[output-directory]',/res,/ps,/data,/log

IRAS21219

IDL> questfit, control file='IRAS21219m1757_dlw_qst.cf',pathin='[path-to-questfit]/questfit/test/',pathout='[output-directory]',/res,/ps,/data,/log

## DETAILED USAGE

Contents of "control file" (*.cf):

A                          B           C     D      E           F    G    H     I    J
source                 spectrum.xdr    -1   -1.   dummy         0.0   0.   X    0.0  0.0
template                    si1.xdr    0.1   0.   DRAINE03      10.0  0.   S    0.0  0.0
absorption                 tau1.xdr    1.0   0.   H2Oice6       0.0   0.   S    0.0  0.0
blackbody                     BB       0.1   0.   DRAINE03      1.0   0.   S	100. 0.0
absorption                 tau1.xdr    0.0   1.   H2Oice6       1.0   0.   S    0.0  0.0
powerlaw                      PL       0.1   0.   DRAINE03      1.0   0.   S   -3.4  0.0
absorption                 tau1.xdr    1.0   0.   H2Oice3       1.0   0.   S    0.0  0.0
absorption                 tau2.xdr    1.0   0.   H2Oice6       1.0   0.   S    0.0  0.0
extinction             draine03.xdr    0.0   0.   DRAINE03      0.0   0.   X    0.0  0.0
(format is not sensitive to exact column width)

Meaning of the columns:

A : The type of data. Include at least one of each datatype
    (template,BB,PL,absorption,extinction) You have to use for at
    least one of each datatype (templates,BB,PL) an absorption (as in
    the file above). And at least one extinction file has to be
    given.

    If you want to use more then one absorption on a certain datatype
    just add them up under the regarding datatype and they will all
    work on the datatype above (like for the powerlaw in the example
    above).

B : in case of source,template,absorption,extinction put in the
    filename.
    in case of BB,PL use a string-dummy, which is just to aid
    your memory and has no effect on program execution.

C : source : this is the lower wavelength limit
            -1 will use the lowest possible common wavelength
            The program will tell you if you use a too small value

   template,blackbody,powerlaw : This is the normalization factor
   
   absorption: This is the taupeak

   extinction: floatpoint dummy (not used)

D : source : this is the upper wavelength limit
            -1 will use the largest possible common wavelength
  	    The program will tell you if you use a too large value

   template,blackbody,powerlaw : This is the fixfree-parameter for the
                                 norm-factor =1 fixed / =0 free while
                                 fitting
 
   absorption: This is the fixfree-parameter for the taupeak =1 fixed
               / =0 free while fitting

   extinction: floatpoint dummy (not used)

E: source,absorption : string-dummy (not used)

   template,blackbody,powerlaw : This is the synonym for the used
   			         extinction curve
 
   extinction : here the synonym for the extinctioncurve is defined
   	        and connected to the xdr.file

F: source,extinction,absorption : floatpoint-dummy 
   template,blackbody,powerlaw  : AV-value for extinction
 
G: source,extinction,absorption : floatingpoint-dummy

   template,blackbody,powerlaw  : fixfree-parameter for AV-value 

H: source,extinction,absorption : string-dummy (not used)

   template,blackbody,powerlaw : S=screen extinction, M=mixed extinction

I: source,template,absorption,extinction : floatpoint-dummy (not used)

   blackbody: This is the temperature in K

   powerlaw : This is the powerlaw-index 

J: source,template,absorption,extinction : floatpoint-dummy (not used)

   blackbody: This is the fixfree-paramter for the temperature

   powerlaw : This is the fixfree-paramter for the powerlaw-index 

After the set up off your control file (saved for example as
control.cf):

IDL> tempfit,control file='control.cf'

The fit will start and produce a plot with the model spectra fitted to
the source. In addition the contribution of each model depending on
wavelength is plotted. This is not in multicolor. If you want to have
a nicer view, create a eps file by typing :

IDL> tempfit,control file='control.cf',ps=1

This creates an eps-file into the 'pathout' directory At the same time
the fitoutput which appears in your shell (after the fit) will be
written into a file if you use ps=1.

The filename for the eps file will be :
fit:sourcename.xdr_par:control file.cf.eps

The filename for the data file will be :
fitresult:sourcename.xdr_par:control file.cf.dat

In the eps plot each data type (template,BB,PL) has a common color.
Within one type the spectra differ from each other by the linestyle.

At the same time on the command shell at the end of the fit some
information appears:

Iter    199   CHI-SQUARE =       165076.83          DOF = 290
    P(0) =              0.00000
    P(1) =              0.00000
    P(2) =              223.071
    P(3) =             0.128759
    P(4) =              0.00000
    P(5) =              0.00000
    P(6) =          2.00909E-13
    P(7) =              1.20905
    P(8) =            0.0300259
    P(9) =              508.749
    P(10) =              44.9747
    P(11) =              27.2432

P(0):  BB :1 fitted extinction (Av):        0.0000000
P(1):  BB :1 fitted taupeak:        0.0000000
P(2):  BB :1 fitted temperature:        222.94748
P(3):  BB :1 fitted normfactor: (/max)       0.12857978
[P(3)]:  BB :1 fitted total normalization factor:   0.00194784

-------------------

P(4):  PL :1 fitted extinction (Av):        0.0000000
P(5):  PL :1 fitted taupeak:        0.0000000
P(6):  PL :1 fitted taupeak:    2.0090873e-13
P(7):  PL :1 fitted powerlaw-index:        1.1869074
P(8):  PL :1 fitted normfactor: (/max)      0.030032732
[P(8)]:  PL :1 fitted total normalization factor:   2.1419353e-18

-------------------

P(9):  template :1 fitted extinction (Av):        508.87861
P(10):  template :1 fitted taupeak:        44.992630
P(11):  template :1 fitted normfactor: (/max)        27.253750
[P(11)]:  template :1 fitted total normalization factor:  1.99867e-09

-------------------

integral [W/cm^2] for source:   2.8561847e-18
integral [W/cm^2] of BB no.:       1 :   2.2627287e-18
contribution [%] to integrated sourceflux :       79.222072
integral [W/cm^2] of PL no.:       1 :   5.3385511e-19
contribution [%] to integrated sourceflux :       18.691197
integral [W/cm^2] of template no.:       1 :   9.6453444e-20
contribution [%] to integrated sourceflux :       3.3770030
Sum of all contributions [%] relative to int. sourceflux :       101.29027


- Above you see all the parameters after the fit.

- The first parameter P(X) is the program internal definition for the
  fitting parameter.

- The numbers of the contributors describe the order of appearance in
  the control file.

- norm factor : =total normalization factor/local maximum.  This is
                the factor you have to put fixed into the control file
                if you like the result. Be careful if you fix the
                parameters and you change the wavelength range.

- total normalization factor: this is the value which you have to
                              multiply to the original template to get
                              the fitted result.

- integral [W/cm^2] for XX: shows the flux integrated over
  	   	    	    frequency-axis


- contribution [%] to integrated source flux : shows the contribution
                                               of the integrated flux
                                               for each contributor

- Sum of all contributions [%] relative to int. sourceflux : adds up
                                                             all
                                                             contributions,
                                                             should be
                                                             100% in
                                                             the case
                                                             of a
                                                             perfect
                                                             fit.


Additional options:

IDL> tempfit,control file='control.cf',res=1 shows the residuum (Data/Model)

IDL> tempfit,control file='control.cf',log=1 uses log xy-axis

IDL> tempfit,control file='control.cf',data=1 will produce output-data
                                        files of each model with
                                        norm-factors not equal to 0.
                                        
                                        The filenames are for example:
                                        sourcename.xdr_par:controlfile.cf_template1.fit
                                        or
                                        sourcename.xdr_par:controlfile.cf_BB2.fit
                                        or
                                        sourcename.xdr_par:controlfile.cf_total.fit
                                        (this is the sum spectrum)
                                        

content of output file : wavelength in µm , flux in Jy

IDL> tempfit,controlfile='control.cf',ident=1 shows an identification index
     					      for the spectra


## QUESTIONS? BUGS? WANT TO MODIFY THE CODE?

Feel free to contact David Rupke at drupke@gmail.com with questions,
bug reports, etc.

Modifications are encouraged, but subject to the license.

## LICENSE AND COPYRIGHT

Copyright (C) 2016--2021 David S. N. Rupke, Vince Viola, Henrik Sturm,
Mario Schweitzer, Dieter Lutz, Eckhard Sturm, Dong-Chan Kim, Sylvain
Veilleux

These programs are free software: you can redistribute them and/or
modify them under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License or any later version.

These programs are distributed in the hope that they will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with these programs.  If not, see http://www.gnu.org/licenses/.
