# Compressible multi-component flow solver 
(latest edit 2023/02/03 fujiwara)

<!--(Comment out）-->
<!--(以下マークダウン記法の基本）
**Boldstyle**
## Titles
`inline`
*  : list
1. : decimal
<br/> : line break
___   : split line
-->

## How to use
1. Compiling source files

    Excute the following command in `solver_3D_omp/Src/`.
    ```html:sample
       make
    ```
    The executor file `job.out` is generated in `solver_3D_omp/Src/`.

2. Setting simulation parameters

    Edit the input file `stdin` in `solver_3D_omp/`.
    The meanings of each params are written [below](#stdin).

3. Start calculation

    Excute the following command in `solver_3D_omp/`.
    ```html:sample
       ./job.out
    ```
    
4. Result

    The calculation results are output to files `output.*******.dat` in `solver_3D_omp/`.
    You can visilize the results with `plot_*****.py` in `solver_3D_omp/plots/`.
    
## Code preview

### stdin
```html:sample
    100000           : iter  
    1000             : nout  
    0                : restart( 0: init calc, 1: restart flow data) -> memo : "from ~~ step" 
    101,3,3          : jmax, kmax, lmax 
    1.0d-4           : dt 
    1.0d0            : Mach number (Mach)   
    0.7d0            : Prandtl number (Pr)   
    100              : Reynolds number (Re) 
    2                : num species 
    5                : rhs (1:divergence, 2:upwind, 3:slau, 4:KEEP, 5:proposed ...) 
    0                : muscl (1:use, 0:not)
    2                : lhs (1:1steuler, 2:RK4 ) 
    0                : vflag (0:Euler , 1:Euler + LAD, 2:physical NS (in progress)) 
    1                : dim_outf (1:1D -> j, 2:2D -> j,k, 3:3D -> j,k,l)
    0                : mconst (1:Molculer weight const, 0:not)  
    3,3,3,3,3,3      : boundary condition ( 1:outflow, 2:inflow, 3,periodic )
```

- num species

    Define the number of gas species.
    The molar weight & specific heat ratio of species is defined at `read_data.f90` for 5 gas (N2, He, CO2, CH4, SF6).
    You can edit these property.
    
- rhs 

    Select the numerical flux used at the advection term in the right hand side.
    You can edit & add schemes in `scheme.f90`.
    
- muscl 

    Switch MUSCL scheme flag for high order accuration. 
    MUSCL is only used when use upwind schemes.
    
- vflag 

    Select the evaluation of the viscous term. 

- dim_outf 

Select output data type defiend at `outf.f90`. 

### Source
- main.f90
- read_data.f90
- init.f90
- step.f90
- bc.f90
- flux.f90
- visflux.f90
- muscl.f90
- outf.f90
- sumdf.f90

### modules
- scheme.f90
- param.f90
