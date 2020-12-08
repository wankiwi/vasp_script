# vasp_script  
Some scripts to perform post-processing or generate input files on VASP calculation.  
***
### neb_generate.py  
```
usage: neb_generate.py [-h] [-v] [-i INITIAL_STATE_CARFILE] [-f FINAL_STATE_CARFILE] [-m {line,idpp}] [-n NUMBER_OF_IMAGES]
```
Takes initial and final CARfiles, and linear or idpp method interpolationthe specified number of images between them.
The interpolated files are written to thedirectories 00 to NI+1, where NI is the number of specified images.

optional arguments:  
  **-h, --help**            show this help message and exit  
  **-v, --version**         Display version  
  **-i INITIAL_STATE_CARFILE, --initial_state_carfile INITIAL_STATE_CARFILE**  
                        The storage location of initial state CARfile. [Optional] [default="is\CONTCAR"]  
  **-f FINAL_STATE_CARFILE, --final_state_carfile FINAL_STATE_CARFILE**  
                        The storage location of final state CARfile. [Optional] [default="fs\CONTCAR"]  
  **-m {line,idpp}, --interpolation_method {line,idpp}**  
                        The method of interpolation. [Optional] [default="idpp"]  
  **-n NUMBER_OF_IMAGES, --number_of_images NUMBER_OF_IMAGES**  
                        The number of interpolation. [Optional] [default=5]             
***
