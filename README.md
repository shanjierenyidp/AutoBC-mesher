# AutoBC-mesher
## Introduction

This mesher generate 2d unstructured mesh with boundary condition (BC) given geometry information and BC information.

![airfoil](/Users/chengruisun/Documents/GitHub/AutoBC-mesher/airfoil.png)

## Get Started

The installation of AutoBC-mesher is as easy.

First, download the .zip file that contains mesh.py 

Then Extract the file to the directory you like

Nevigate into the directory, type ` pwd` 

Do

` export PATH=$PATH:/Path/to/the/output/of/pwd `	

Then type 

` mesh.py -h`

to see if it works

## Tutorial

To see all the arguments, go to terminal window and type:

` python mesh.py -h`

and following is the detailed explaination

### Input file `-ifile`

To specify input file, type

` -ifile path/to/the/input/file` 

For example, use

` -ifile /aorta/aorta.npy  `

**Note that the shape of the input geometry must be a n*2 matrix**

### Output file `-ofile`

To specify output file, type

` -ofile path/to/the/output/filename.format` 

For example,use

` -ofile ./aorta.stl  `

The format of the ouput can be chosen among:

- .stl: ` -ofile ./aorta.stl  `
- .su2: ` -ofile ./aorta.su2  `
- openfoam polymesh: ` -ofile ./aorta  `
- .vtk: ` -ofile ./aorta.vtk`
- .vtp: ` -ofile ./aorta.vtp  `

### Segments `-s`

Argument segment represents the point index at which the boundary ends

For example, use 

`-s 1 10 109 118 `

### Resolution ` -r`

By specifying resolution, the resolution of each boundary is defined

For example, use

`-r 30 50 30 50`

**Note the resolution must be an integer**

### Labels `-l`

By specifying labels, the boundary label is determined

For example, use

`-l  inlet upperWall outlet lowerWall`

### Max area `-ma`

The maximum area of the units in the mesh, the max area can be a float number

For example, use

`-ma 4e-4`

### Number of inner loops `-nl`

This arguments specifies the number of inner loops in the geometry. The default value of this argument is 0

After specifying the number of loops, the user will be required to give the parameters neccessary to specify each inner loop, for example:

` Please enter the points of loop, seperated by space./airfoil/airfoil.npz` 

`
Please enter the segments of loop, seperated by space1`

`Please enter the labels of loop, seperated by spaceairfoil`

`
Please enter the resolutions of loop, seperated by space200`

#### 

