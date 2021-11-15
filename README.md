# FRobUP_body
This repository contains a MATLAB implementation to compute and evaluate frequency robust kT-points pTx pulses in the human heart.

##### Authors:
- Christoph S. Aigner  (<christoph.aigner@ptb.de>)
- Felix Krüger  (<felix.krueger@ptb.de>)
- Max Lutz  (<max.lutz@ptb.de>)
- Sebastian Dietrich   (<sebastian.dietrich@ptb.de>)
- Sebastian Schmitter  (<sebastian.schmitter@ptb.de>)

Usage
--------

Run script main_frequRobUPs.m


Contents
--------

##### Test scripts (run these):
    main_frequRobUPs.m          test script to compute and evaluate tailored respiration specific and respiration robust pulses at 7T

##### Routines called by the test scripts:
    cmap.mat                256x3 double matrix with the used colormap of the manuscript
    kTrandphases.mat        200x8 double matrix with randomized phase initials
    
##### External data files used by the test scripts:
    B1R.zip          

Dependencies
------------
These routines were tested under MATLAB R2019a under Windows, but should also run under older versions.

The 36 channel-wise invivo B1+ datasets of the human body at 7T are available at: TBD and were computed as described in [2].

The optimization of the kT-points is performed using code by Will Grissom and Zhipeng Cao ([3,4] and https://bitbucket.org/wgrissom/acptx/) who have given permission for inclusion within this package. 

Please cite appropriately.

License
-------

This software is published under GNU GPLv3. 
In particular, all source code is provided "as is" without warranty of any kind, either expressed or implied. 
For details, see the attached LICENSE.

Reference
---------

1) Christoph S. Aigner, Sebastian Dietrich, Felix Krüger, Max Lutz and Sebastian Schmitter, Towards frequency robust tailored and universal pulses in the human heart at 7T, submitted to ISMRM 2022

Created by Christoph S. Aigner, PTB, November 2021.
Email: christoph.aigner@ptb.de
