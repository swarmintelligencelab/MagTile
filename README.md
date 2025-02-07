# MagTile

This repository contains the hardware and software implementation of the MagTile actuation system. Additionally, it includes a computational model of the system implemented in MATLAB.

## Repository Structure
The repository is organized into four main subdirectories, each containing the necessary information and code to reconstruct and operate the actuation system.

### **1. Data**
This subdirectory contains the raw data used for system calibration, ensuring accurate and repeatable actuation behavior.

### **2. Hardware**
This subdirectory includes the physical design files and electronic schematics required for building the MagTile system:
- **`cad/`**: Contains all CAD files in STL format, representing the 3D-printed components used in the assembly of the actuation system.
- **`pcb/`**: Contains detailed schematics outlining the design of the printed circuit board (PCB) used for actuating and controlling the solenoids.

### **3. Computational Implementation**
This subdirectory contains all the MATLAB codes necessary for the computational implementation of MagTile. These scripts allow for the simulation and control of the system's behavior

### **4. Software Implementation**
This subdirectory includes all scripts required for real-time operation and control of the actuation platform:
- **Python Scripts**: Responsible for tracking and streaming data.
- **Arduino Scripts**: Handle the control logic for operating the solenoids.

## Usage
1. Navigate to the relevant subdirectory for the required functionality (e.g., computational implementation, hardware design, or software operation).
2. Run the MATLAB scripts to simulate the system's behavior.
3. Utilize the Python and Arduino scripts to control the actuation system in real-time.

For detailed implementation and additional information, refer to the accompanying research paper.

