# Compressible multi-component flow solver 
(latest edit 2022/12/07 fujiwara)

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

<!--
[![MemoTube](https://user-images.githubusercontent.com/28818747/98420702-13594700-20cb-11eb-9268-8c304fdb7cb2.png)](https://tohoku.memotube.xyz)<br/><br/>
-->

## How to use
1. Compiling source files

    Excute the following command in `solver_omp/Src/`.
    ```html:sample
       make
    ```
    The executor file `job.out` is generated in `solver_omp/Src/`.

2. Setting simulation parameters

    Edit the input file, `stdin`, in `solver_omp/`.
    The meanings of each params are written [below](#stdin).

3. Start calculation

    Excute the following command in `solver_omp/`.
    ```html:sample
       ./job.out
    ```
    
4. Result

    The calculation results are output to files `output.*******.dat` in `solver_omp/`.
    You can visilize the results with `plot_*****.py` in `solver_omp/plots/`.
    
## Code preview
  ### stdin
  ### Makefile
  ### Source
    #### main.f90
    #### read_data.f90
    #### init.f90
    #### step.f90
    #### bc.f90
    #### flux.f90
    #### visflux.f90
    #### muscl.f90
    #### outf.f90
    #### sumdf.f90
  ### modules
    #### scheme.f90
    #### module.f90
    
