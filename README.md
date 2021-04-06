QUESTFIT relies on some subroutines from the ISAP repository, which
can be downloaded here:

https://old.ipac.caltech.edu/iso/isap/isap.html

[Note that the Spitzer IRS data reduction software SMART ships with a
slightly modified version of the ISAP library. At least one instance
of errors has been documented using this version of the library.]

Some tests:

IDEOS ID 4978688_0 = Mrk 231

IDL> questfit, controlfile='4978688_0.ideos.cf',pathin='[path-to-questfit]/questfit/test/',pathout='[output-directory]',/res,/ps,/data,/log

IRAS21219

IDL> questfit, controlfile='IRAS21219m1757_dlw_qst.cf',pathin='[path-to-questfit]/questfit/test/',pathout='[output-directory]',/res,/ps,/data,/log
