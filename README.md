# OVR 
![example workflow](https://github.com/rugtres/ovr/actions/workflows/build.yml/badge.svg)

Oncolytic Virus Resistance simulation

Repository for code and executable simulating the interplay between a tumour and an Oncolytic Virus, including anti-viral resistance.

## GUI

The GUI version is available as Windows installer in the folder 'release'. Furthermore, the folder GUI contains a QT project file and the necessary ui files to be able to compile the program

## Console

For demanding simulations, it can be beneficial to use the program without the GUI components, in order to improve computational speed. For this purpose, the files in the folder 'console' are available. These include some code to do automated analyses and code to read parameters from a config file (in contrast to reading them from the GUI). One can use the QT project file, but alternatively the folder also includes a Makefile to compile the code. 

## Testing

Testing is performed using the Catch2 framework, based on the console version. The 'test' folder contains all relevant files to be able to perform testing. 

Branch|Code Coverage
---|---
master|[![codecov.io](https://codecov.io/github/rugtres/ovr/coverage.svg?branch=master)](https://codecov.io/github/rugtres/ovr/coverage.svg?branch=master)


![](https://github.com/rugtres/ovr/blob/main/Screenshots/GUI_tabbed.png)
