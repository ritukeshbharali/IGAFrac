**IGAFrac** is a MATLAB code for simulating fracture within the Iso-Geometric Analysis paradigm. The geometry is represented by PHT-splines (polynomial splines over hierarchical T-meshes), and is based on the [nurbs](https://octave.sourceforge.io/nurbs/index.html) package (distributed under the GNU General Public License v3.0). The code adopts a modular style, similar to that of [PyFEM](https://github.com/jjcremmers/PyFEM). 

The following models are available:

#### IGA Models 
  - Linear Elasticity
  - Phase-field Fracture

#### Materials
  - Isotropic Elastic

#### Nonlinear solver
  - Newton-Raphson
  - Staggered
  - Quasi Newton-Raphson with extrapolation
  - Arc-length

#### Other options
  - Realtime output
  - VTK writer

### Getting started
  - Run one of the examples in 'inputData' folder!  

### Contribute
  - Pull request/Issues

Please consider citing the following paper if you use our code:
  - Bharali et. al. (2022) A robust monolithic solver for phase-field fracture integrated with fracture energy based arc-length method and under-relaxation. [link](https://www.sciencedirect.com/science/article/pii/S0045782522001992)
