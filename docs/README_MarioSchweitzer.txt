 Identifier   tempfit

 Purpose      fit spectra with BB + powerlaws + templates 
              (additional use of individual extinctions and individual absorptionfeatures is possible)

              BEWARE: by default the data points are weighted by 1/stdev^2. If you do not absolutely
              trust the stdev values in your data structure, you may need to set stdev by different code
              before using tempfit!

              Templates, blackbodies etc. are rebinned to the sampling of the source spectrum before
              fitting
              
 Synopsis     tempfit,controlfile=controlfile,[res=res],[ps=ps],[data=data],[log=log],[ident=ident]

 Some information about the program and its use can be drawn from the
 header of tempfit.pro. Here I will give a general guide how to use it. The
 generated fitfunction is described in fitfunction.pdf

 Version of this help file: Novemeber 21, 2006


IDL files needed for fitting:
-----------------------------

 - tempfit.pro     main fittingprogram

 - readcf.pro      reads the data from controlfile

 - startup.dat     contains the pathsettings

 - weight.pro     (not nessesary but nice for weighting of certain
                   wavelengthregions without changing stdev values outside of
                   defined wavelengthrange, Can also be used to exclude wavelengthranges)
                   DESCRIPTION SEE UNDER ADDITIONAL OPTIONS

 - weightr.pro    (not nessesary but nice for relative weighting of certain
                   wavelengthregions [will set all stdev values outside of
                   defined wavelengthrange to 1.])
                   DESCRIPTION SEE UNDER ADDITIONAL OPTIONS

 - weightres.pro  (not nessesary but nice for weighting with respect to the estimated
                   spectral resolution)
                   DESCRIPTION SEE UNDER ADDITIONAL OPTIONS

 - dlweight.pro   (not nessesary but nice for weighting taking the spectral
                   slope into account. Weighting for lower fluxes higher then
                   for higher fluxes using a powerlawfit to the spectrum)
                   DESCRIPTION SEE UNDER ADDITIONAL OPTIONS


 - chiplot.pro    (not nessesary but nice for creating a chi^2 plot with changing
                   fitparameters one- or two-dimensional optional )
                   DESCRIPTION SEE UNDER ADDITIONAL OPTIONS

               

 - extinction files 

 - absorption files 

 - template files (e.g. silicate emission models)

 for information about individual templates see the files:

 INFO.txt           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 silicatemodels.txt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 in the data_in directory  


 - the controlfile: control.cf

 - source file (of course)

  

 Which format should the files have ?
 
 - All files are given in .xdr format as produced with smart
 - The program automatically estimates the common wavelengthregion to do the fit 



 What you need (to do) before you can do the fitting:
-------------------------------------------------------------------------------------------

The fitting-package includes all nessesary pro-files from IPAC and
ASTROLIB. Put them into your IDL-path !

The "data_in" directory includes all files :templates, absorptionfeatures,
silicate models, extinctioncurves and should contain the source-xdrs you wanna fit !

For testing two QSOs ( 1ZW1 and PG1440+356) are included (both without lines and
in restframe)
Also for testing three ULIRGs are included (IRAS00091m0738nl.xdr,IRASF00456m2904SWnl.xdr
IRASF01166m0844_SEnl.xdr) in the "data_in" directory

Change the path in the startup file : startup.dat
All *.pros will read the path from here ! 

 - use e.g. emacs to put in your input path :pathin /.../data_in
                          your output path :pathout /.../data_out


example.cf is an example for a controlfile you can use and change as you want!

    
Set up the controlfile control.cf:
---------------------------------------------------------------------------------------------
Use e.g. emacs to produce a file like:


  A                          B          C     D      E           F    G    H     I    J

source       RebinNGC3998_final.xdr    -1   -1.   dummy         0.0   0.   X    0.0  0.0
template                    si1.xdr    0.1   0.   DRAINE03      10.0  0.   S    0.0  0.0
absorption                 tau1.xdr    1.0   0.   H2Oice6       0.0   0.   S    0.0  0.0
blackbody                     BB       0.1   0.   DRAINE03      1.0   0.   S	100. 0.0
absorption                 tau1.xdr    0.0   1.   H2Oice6       1.0   0.   S    0.0  0.0
powerlaw                      PL       0.1   0.   DRAINE03      1.0   0.   S   -3.4  0.0
absorption                 tau1.xdr    1.0   0.   H2Oice3       1.0   0.   S    0.0  0.0
absorption                 tau2.xdr    1.0   0.   H2Oice6       1.0   0.   S    0.0  0.0
extinction             draine03.xdr    0.0   0.   DRAINE03      0.0   0.   X    0.0  0.0

(format is not sensitive to exact column)

Meaning of the columns:
-------------------------------------------------------------------------------------------


A : contains the type of data.
    IMPORTANT: right now it is important to include at least one of each
    datatype (template,BB,PL,absorption,extinction)
    You have to use for at least one of each datatype (templates,BB,PL) an
    absorption (like done in the file above). And at least one extinctionfile has to
    be given. If not an ERROR will occur.

    If you want to use more then one absorption on a certain datatype just add
    them up under the regarding datatype and they will all work on the
    datatype above (like for the powerlaw in the example above).


B : in case of source,template,absorption,extinction put in the filename
    in case of BB,PL use a string-dummy, which is just to aid your memory and has
    no effect on program execution.


C : source : this is the lower wavelength limit
            -1 will use the lowest possible common wavelength
            The program will tell you if you use a too small value

   template,blackbody,powerlaw : This is the normalization factor
   
   absorption: This is the taupeak

   extinction: floatpoint dummy (not used)


D : source : this is the upper wavelength limit
            -1 will use the largest possible common wavelength
  	    The program will tell you if you use a too large value

   template,blackbody,powerlaw : This is the fixfree-parameter for the norm-factor
                                 =1 fixed / =0 free while fitting
 
   absorption: This is the fixfree-parameter for the taupeak
                            =1 fixed / =0 free while fitting

   extinction: floatpoint dummy (not used)


E: source,absorption : string-dummy (not used)

   template,blackbody,powerlaw : This is the synonym for the used extinction curve
 
   extinction : here the synonym for the extinctioncurve is defined and
   connected to the xdr.file


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



DOING THE FIT
------------------------------------------------------------------------------------------------

After the set up of your controlfile (saved for example as control.cf) switch to idl and type:


tempfit,controlfile='control.cf'

The fit will start and produce a plot with the modelspectra fitted to the
source. In addition the contribution of each model depending on wavelength is
plotted. This is not in multicolor. If you want to have a nicer view, create a
eps file by typing :

tempfit,controlfile='control.cf',ps=1

This creates an eps-file into the 'pathout' directory
At the same time the fitoutput which appears in your shell (after the fit) will
be written into a file if you use ps=1.
The filename for the eps file will be  : fit:sourcename.xdr_par:controlfile.cf.eps
The filename for the data file will be : fitresult:sourcename.xdr_par:controlfile.cf.dat

In the eps plot each data type (template,BB,PL) has a common color.
Within one type the spectra differ from each other by the linestyle.

(There are still some color problems, that means it could be that after
producing the eps files the plots on your screen may change their color)


At the same time on the command shell at the end of the fit
some information appear:

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
reached end of program



- Above you see all the parameters after the fit !

- The first parameter P(X) is the program internal definition for the
  fittingparameter. You have to use this parameter with the chiplot command
  to create a chi^2 plot for the regarding parameter P(X)

- The numbers of the contributors describe the order of appearance in the
  controlfile.

- norm factor : =totalnormalizationfactor/local maximum. 
                 This is the factor you have to put fixed into the controlfile
                 if you like the result. Be careful if you fix the parameters
                 and you change the wavelength range (see ATTENTION)

- total normalizationfactor means: this is the value which you have to multiply
                                   to the original template to get the fitted
                                   result.

- integral [W/cm^2] for XX: shows the flux integrated over frequency-axis 


- contribution [%] to integrated sourceflux : shows the contribution of the
                                              integrated flux for each contributor

- Sum of all contributions [%] relative to int. sourceflux : adds up all contributions, should be 
                                                             100% in the case of a perfect fit.


- Sometimes some Error messages can follow : I think they are due to too large
  or small values in the plot commands. They should not matter. If "reached
  end of program" appears the program did not interrupt and everything should
  be ok!



Additional options:
-----------------------------------------------------------------------------------------------

tempfit,controlfile='control.cf',res=1 shows the residuum (Data/Model)

tempfit,controlfile='control.cf',log=1 uses log xy-axis

tempfit,controlfile='control.cf',data=1 will produce output-data files of each
                                        model with norm-factors not equal to
                                        0.
                                        
                                        The filenames are for example:
                                        sourcename.xdr_par:controlfile.cf_template1.fit
                              or        sourcename.xdr_par:controlfile.cf_BB2.fit
                              or
                                        sourcename.xdr_par:controlfile.cf_total.fit
                                        (this is the sum spectrum)
                                        

content of output file : wavelength in µm , flux in Jy

tempfit,controlfile='control.cf',ident=1 shows an identification index for the spectra






WEIGHT.PRO

Weighting of certain wavelengthregions 
------------------------------------------------------------------------------------
Weight.pro sets the stdev values for certain wavelengthregion to a specified
value. All other values (if not 0.) will not be touched. You have to choose a
wavelengthregions and apply any stdev-value for the values within this
regions. (The program automatically sets stdev-values which are 0. to 1.)

write e.g. :weight,sourcename='sourcefile.xdr',x=10.,y=20.,weight=2.

a new file named weight:sourcefile.xdr is generated in the input directory.
The stdev values for the wavelength values between 10 and 20 um are set to
2. All other stdev values are not changed !

Now put weight:sourcefile.xdr instead of sourcefile.xdr into your controlfile
and do the fit again.

To add a 2nd, third ... region for weighting use the weightcommand again like:

weight,sourcename="weight:sourcefile.xdr",x=20.,y=22.,weight=3.

weight,sourcename="weight:sourcefile.xdr",x=22.,y=24.,weight=4.

The result will be:

stdev(10-20)=2. , stdev(20-22)=3., stdev(22-24)=4.  

If you want to choose a new region for weighting write again:

weight,sourcename="sourcefile.xdr",x=24.,y=30.,weight=5.

The result will be:

stdev(24-30)=5 (all other stedev values are not touched)

If you dont know the limits for the wavelengthrange just put too small or to
large values for x or y and the program tells you the min,max values of the
wavelengthaxis.

You can use weight.pro to exclude a certain wavelengthrange from the fit by
setting the stdev values for this range high (e.g. 1e7)






WEIGHTR.PRO

Relative weighting of certain wavelengthregions
-----------------------------------------------------------------------------------
By hand one can overwrite the original stdev-values. If you use weightr.pro 
all stdev values will be set to 1. In addition you have to choose a
wavelengthregions and apply any stdev-value for the values within this
regions. (The program automatically sets stdev-values which are 0. to 1.)

write e.g. :weightr,sourcename='sourcefile.xdr',x=10.,y=20.,weight=2.

a new file named weightr:sourcefile.xdr is generated in the input directory.
The stdev values for the wavelength values between 10 and 20 um are set to
2. All other stdev values are now 1 !

Now put weightr:sourcefile.xdr instead of sourcefile.xdr into your controlfile
and do the fit again.

To add a 2nd, third ... region for relative weighting use the weightrcommand again like:

weightr,sourcename="weightr:sourcefile.xdr",x=20.,y=22.,weight=3.

weightr,sourcename="weightr:sourcefile.xdr",x=22.,y=24.,weight=4.

The result will be:

stdev(10-20)=2. , stdev(20-22)=3., stdev(22-24)=4.  all other stdev's=1

If you want to choose a new region for weighting write again:

weightr,sourcename="sourcefile.xdr",x=24.,y=30.,weight=5.

The result will be:

stdev(24-30)=5 all other stdev's are now 1 !

If you dont know the limits for the wavelengthrange just put too small or to
large values for x or y and the program tells you the min,max values of the
wavelengthaxis.






WEIGHTRES.PRO (!!!! has to be checked !!!!)

Weighting by using the spectral resolution of the wavelengthregion estimated
from the wavelength-structure of the source xdr.
-----------------------------------------------------------------------------------------------

write e.g. :weightres,sourcename='sourcefile.xdr'

the program will set the stdev-values to normalized weights (normalized to
maximum R)  estimated from the spectral resolution (R=lambda/deltalambda)

the way how the stdev values are produced is:

sourcestdev(i)=1/((sourcewave(i+1)-sourcewave(i))/sourcewave(i))
sourcestedev=sourcestdev/max(sourcestdev)

replace the sourcename in your controlfile by weightres:sourcename.xdr
and do the fit again







CHIPLOT.PRO 

Creates a chi^2-plot with variing parameter P(X)
(one or two parameters optional)
-----------------------------------------------------------------------------------------------

After you used tempfit.pro the session information (variables) are stored in
the data_out directory in fitvariables.sav

This file is restored by chiplot. So each time you do a new fit, the session
information is updated.

1 parameter:
------------

write for example:

chiplot,para=2.,percenta=100.

The program will do a chi^2 plot of parameter P(2) which in the upper example
would be the BB temperature. The value 100 means that P(2) will be changed by +/-
100 % of its own value. The default is 50 % .
A plot with the fitresult will appear and show you how the spectral shape
changes if P(2) changes by +/- 100%. After this the chi^2 plot is produced.
The original value of P(2) is written in the title of the plot as well as in
the shell.


2 parameters :
--------------

write for example:

chiplot,para=2.,parb=3.,percenta=100.,percentb=200.

A chi^2-surface plot will be generated variing P(2) by +/- 100% and P(3)
(which is the normalization factor of the BB in the upper example) by +/- 200%
By moving the mouse on the plot and clicking the left mousebutton one can
change the line of sight. By clicking the right mousebutton a contourplot will
be generated. By moving the mouse on the plot and clicking the left mousebutton one can
change the scaling of the contourlines. Right click on the plot will quit the
program.







DLWEIGHT.PRO 

Weighting by taking the spectral slope into account 
---------------------------------------------------

Reads input spectrum. If z is specified, wavelength scale is first converted
to rest frame. Stdev values are set so that for 1/stdev^2 weighting,a 
similar relative flux deviation delta_flux/flux over a similar resolution
element delta_lambda/lambda will give the same chi2 change. The flux trend in this
computation is from a large scale power law fit and does NOT trace individual features.
Returns an output xdr with "dlweight:" added to the filename

write for example:

dlweight,'sourcename.xdr',z=0.123 (if not given in restframe)

then exchange the sourcename in your controlfile by
dlweight:sourcename.xdr and do the fit again.










IMPORTANT COMMENTS:
-----------------------------------------------------------------------------------------------



ATTENTION:

fix parameters:
---------------
If you change the wavelengthrange and you want to keep the old fitresult
(e.g. from a smaller wavelength range) you must not set the fitparameters
fixed. Since each contributor is divided by its local maximum (within the
present wavelengthrange) the old fitresult may change due to a greater maximum
which could appear within the larger wavelength range.
If you want to see how the fit result from a smaller range (without any change)
looks on a larger wavelength intervall set the fitparameters free and set the
stdev values of the sourcefile for the wavelengthrange outside the old
(smaller) range to high (e.g. 10e6) values (using exclude.pro). 
This allows the program in the presence of a new greater maximum to optimize the
fitparameters again (same result as before => same total normalizationfactor
but different norm-factor) without paying attention to the wavelengthrange outside the old one.
   

(This could be changed by using the global maximum for normalization ! But if
the flux of the spectra varies over some orders it may be better to use the
local maximum for an easier estimation of the startparameter in limited
wavelength regions)
     

the chi^2 value:
----------------

The chi^2 value will not be in its typical range since we don t know the real
propagated errors. On the other hand we are not taking a combined error into
account which would also pay attention to the errors within the templates.
Nevertheless the chi^2 value is an indicator of the quality of the fit. 
Then lower it is (in a given configuration) then better the fit is.
Pay attention if you change the weighting.Then the chi^2 value will not be
comparabel to the one with old weighting !!!



Have fun ...
