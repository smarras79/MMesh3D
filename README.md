MMesh3D
--------
Elliptic, structured grid generator for simply connected domains with real topography:

##Compile and run
Tested on Mac OS X and multiple versions of Ubuntu using gcc-4* and later and ifort 8.1

*Compile:
```
cd MMesh3D/src
>> make
```

*Run
```
cd MMesh3D/runs
>> ./MMesh3D-V2.a Input_meshparam.inp
```


## Used or mentioned in:
Some mention of MMesh3D:

"Only MMesh3D provides development and testing setup information. The developer indicates that MMesh3D has been installed
without problems on Mac 10.5.4 and Ubuntu 8.10 with gcc version
4.3.2. The other software products do not provide this information,
so other developers or potential users cannot be certain that their
results are reproductions of the creator’s original results. Furthermore, there is no evidence of automated tools to capture the experimental context so “works-on-my-computer” issues cannot be
easily diagnosed" @smithEtAl201

@article{meierbachtolEtAl2017,
  author = {Meierbachtol, C.S. and Svyatskiy, D. and Delzanno, G.L. and Vernon, L. and Moulton, J.D.},
  title = {An electrostatic Particle-In-Cell code on multi-block structured meshes},
  journal = {J. Comput. Phys.},
  year = {2017},
  volume = {350},
  pages = {796-823},
},

@article{luEtAl2017,
  author = {Lu, F. and Pang, Y. and Jiang, X. and Sun, J. and Huang, Y. and Wang, Z. and Ju, J.},
  title = {Automatic generation of structured multiblock boundary layer mesh for aircrafts},
  journal = {Advances in Engineering Software},
  year = {2017},
  volume = {115},
  pages = {},
  doi = {https://doi.org/10.1016/j.advengsoft.2017.10.003},
},

@article{smithEtAl201,
  author = {Smith, W.C. and Lazzarato, D.A. and Carette, J.},
  title = {State of the practice for mesh generation and mesh processing software},
  journal = {Advances in Engineering Software},
  year = {2016},
  volume = {100},
  pages = {},
  doi = {https://doi.org/10.1016/j.advengsoft.2016.06.008},
},

@incollection{ICMEhandbook,
   title = {Numerical Methods}
   author = {Agelet de Saracibar, C. and Boman, R. and Bussetta, P. and Cajas, J.C. and Cervera, M. and Chiumenti, M. and Coll, A. and Dadvan, P. and Hernand\'ez Ortega, J.A. and Houzeaux, G. and Pasenau de Riera, M.A. and Ponthot J.P.}, 
   editor = {Schmitz, G.J. and Ulrich, P.},
   booktitle = {{Handbook for software solutions for ICME}},
   publisher = {Wiley-VCH},
   year  = {2017},
   pages  = {487--532},
},

@incollection{joque2017,
   booktitle = {Creating Data Literate Students},
   author = {Joque, J},
   editor = {Fontichiaro, K. and Oehrli, J.A. and Lennex, A.},
   title = {Chapter 6: Making sense of data visualization},
   publisher = {Michigan Publishing, University of Michigan Library},
   year  = {2017},
   doi={http://dx.doi.org/10.3998/mpub.9873254},
},

@article{horgueEtAl2018,
  author = {Horgue, P. and Franc, J. and Guibert, R. and Debenest, G.},
  title = {{An extension of the open-source porousMultiphaseFoam toolbox dedicated to groundwater flows solving the Richards equation}},
  journal = {arXiv:1510.01364v1 [cs.CE]},
  year = {2018},
  volume = {},
  pages = {},
},

@misc{robertschneidersWeb
    author = {Schneiders, R.},
    title = {www.robertschneiders.de/meshgeneration/software.html},
    howpublished = {web},
},

