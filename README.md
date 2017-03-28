# JordanBransDicke
This contains a patch for RAMSES to run any model that has a modified H(a) and GeffG(a). We include a code to compute the background and scalar-field evolution of the simple Jordan-Brans-Dicke model.

To use the RAMSES patch simply place the folder [ramses\_patch/modified_geff_hubble] in the [ramses/patch] directory
and compile using the Makefile provided here. Alternatively modify your makefile by adding:
 - [PATCH = ../patch/modified_geff_hubble]
 - [-DMODIFIEDGEFFHUBBLE] to DEFINES 
 - [splines.o modified_geff_hubble.o] to MODOBJ

To use it add: [filename_geff_hubble_data='/path/to/file.txt'] to AMR_PARAMS in your namelist file.
This file must contain a header-line, then the number of lines and then [log(a_i),  GeffG(a_i)  H(a_i)/H0] with
a_i in increasing order, a_{i+1} > a_i. Alternatively modify [modified_geff_hubble.f90] to specify the
two functions directly.
