MMesh3D
--------
Elliptic, structured grid generator for simply connected domains with real topography:

(Created in 2011. Posted on GitHub in 2020)

##MPI version:
#Option 1 (Recommended): use homebrew to install gcc and openmpi
```
>> brew install gcc openmpi metis
```
Brew will install them in `/opt/homebrew/Cellar` or similar.

#Option 2: Open-MPI 4.1.0 compiled from source with gcc-9

```
>> ./configure --prefix=/MY/PATH/TO/openmpi-4.1.0/build_gcc9/ CC=gcc-9 CXX=g++-9 F77=gfortran-9 FC=gfortran-9
>> make -j 6 all
>> make install
```

## Compile and run
Tested on Mac OS X and multiple versions of Ubuntu using gcc-4* and later and ifort 8.1

*Compile
```
cd MMesh3D
>> make
```

*Run
```
>> cd MMesh3D/runs
>> mpirun -np 4 ./MMesh3D.a Input_meshparam.inp
```

*Run with GDB debugger:
First make sure that you have the proper debug flags activated in your compilation flags (-g and others)
```
>> cd MMesh3D/run
>> gdb --args ./MMesh3D.a Input_meshparam.inp
>> (gdb) run
```

*Visualization
```
MMesh3D generates a VTK file that can be visualized with ParaView or VisIt.
```

NOTICE: The PDF manual is for V1.0 and still works with the current V2.0 that is in this repo. However, it does not contain information on the elliptic smoothing. 
I will post instructions as soon as I have some spare time. 

## License
MMesh3D is a free and open source software according to the Free Software Foundation under a GPL license which permits any use or modification of MMesh3D as long as any derived work is open-source and preserves the same terms. For commercial, closed-source use, such as distributing a binary executable without the source, contact me directly to discuss options.

It fulfills the definition of Open Source as per the OSI definition: https://opensource.org/osd



## Citing MMesh3D
Please, cite MMesh3D by citing:
```
@phdthesis{marras2012variational,
  title={Variational multiscale stabilization of finite and spectral elements for dry and moist atmospheric problems},
  author={Marras, Simone},
  year={2012},
  school={PhD thesis, Universitat Polit{\'e}cnica de Catalunya}
},

@article{marras2016review,
  title={A review of element-based Galerkin methods for numerical weather prediction: Finite elements, spectral elements, and discontinuous Galerkin},
  author={Marras, Simone and Kelly, James F and Moragues, Margarida and M{\"u}ller, Andreas and Kopera, Michal A and V{\'a}zquez, Mariano and Giraldo, Francis X and Houzeaux, Guillaume and Jorba, Oriol},
  journal={Archives of Computational Methods in Engineering},
  volume={23},
  number={4},
  pages={673--722},
  year={2016},
},
```

## Used or mentioned in:

```
"Only MMesh3D provides development and testing setup information. 
The developer indicates that MMesh3D has been installed
without problems on Mac 10.5.4 and Ubuntu 8.10 with gcc version
4.3.2. The other software products do not provide this information,
so other developers or potential users cannot be certain that their
results are reproductions of the creator’s original results" [Smith et al. 2016]
```

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

@article{smithEtAl2016,
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

