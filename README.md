# OVI 
![example workflow](https://github.com/rugtres/OVI/actions/workflows/build.yml/badge.svg)

Oncolytic Virus Immune simulation

Repository for code and executable simulating the interplay between a tumour and an Oncolytic Virus, including anti-viral resistance and immune response.

## GUI

The GUI version is available as Windows installer in the folder 'release'. Furthermore, the folder GUI contains a QT project file and the necessary ui files to be able to compile the program

## Console

For demanding simulations, it can be beneficial to use the program without the GUI components, in order to improve computational speed. For this purpose, the files in the folder 'console' are available. These include some code to do automated analyses and code to read parameters from a config file (in contrast to reading them from the GUI). One can use the QT project file, but alternatively the folder also includes a Makefile to compile the code. 

## Testing

Testing is performed using the Catch2 framework, based on the console version. The 'test' folder contains all relevant files to be able to perform testing. 

Branch|Code Coverage
---|---
master|[![codecov](https://codecov.io/gh/rugtres/OVI/branch/main/graph/badge.svg)](https://codecov.io/gh/rugtres/OVI)

## Screenshots

Default coloring of cells:
![](https://github.com/rugtres/OVI/blob/main/Screenshots/GUI.png)

Inflammation factor visualisation:
![](https://github.com/rugtres/OVI/blob/main/Screenshots/Inflammation_factor.png)


