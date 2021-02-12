---
title: 'MMesh3D: structured, elliptic mesh generator for Geophysical applications (and beyond)'
authors:
  - name: Simone Marras^[smarras@njit.edu]
    orcid: 0000-0002-7498-049X
    affiliation: "New Jersey Institute of Technology, California Institute of Technology" # (Multiple affiliations must be quoted)
affiliations:
 - name: Simone Marras, Assistant Professor, NJIT
   index: 1

date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary


# Statement of need 

`MMesh3D` is a simple, automatic, numerical grid generation tool designed
to build logically structured volumeteric grids made of both hexahedra and tringual base prisms.
The code is written in C, commented throughout, and can be easily modified
or enhanced by the user’s own functions. I wrote this code for my own purposes but
tried to write it in such a way that it could be re-usable or changed by other users. The
code is free and can be modified and redistributed under the GNU GENERAL PUBLIC
LICENSE.
The visit_writer library to write the visualization output file is included with the package.
It was written by Hank Childs at the Lawrence Livermore National Laboratory
for the visualization software VisIt [6], and writes grids and data into the Visualization
ToolKit format VTK. The license of this library is reproduced in the header of the library
itself and in the function wrt2plotfile.c that calls it.
The dynamic allocation of the arrays is done on a 1-index base through functions from
the library nrutil.c. This function is part of the NR library and is can be freely redistributed.
The NR library is a numerical library that comes with the successful volume
Numerical Recipes in C: the Art of Scientific Computing [3].


`MMesh3D` was used for different applications that go beyond atmospheric flows for which it was originally designed. We find it used to build boundary layer girds in [@luEtAl2017], for electrostatic simulations that required a free multi-block structured grid generator easy to modify [@meierbachtolEtAl2017], as well as porous medium flows using OpenFOAM [@horgueEtAl2018]. It is mentioned in the following studies and review papers [@joque2017],[@ICMEhandbook], [@smithEtAl2016], [@robertschneidersWeb].

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Fenced code blocks are rendered with syntax highlighting:
```python
for n in range(10):
    yield f(n)
``` 

# Acknowledgements


# References
