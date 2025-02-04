# MagTile
This repo contains the hardware and software implementation of the actuation system. Additionally the repo contains a computational implementation of the system on MATLAB.

A total of 4 subsirectories consolidate the necessary information and code to reconstruct the actuation system. 
- data: The subdirectory contains the raw data used in calibration.
- hardware: The subdirectory contains two folders. The first, titled 'cad', contains all the CAD files in STL format describing the 3D printed parts used in the assembly of the actuation system. The second folder, titled 'pcb', contains a file with the detailed schematics of describing the design of the PCB used to actuate and control the coils.
- simulations: the subdirectory contains all the MATLAB files required to run the simulations used in this study and demonstrated in the paper.
- software_implementation: The subdirectory contains all the scripts necessary to operate the platform, including the python implementation responsible of tracking and streaming of data, and the arduino script operating the coils.
