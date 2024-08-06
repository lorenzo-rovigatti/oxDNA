# Linear dsDNA with sticky ends

Author: Jinhuang Zhou  
Last update: Aug 2024

The config and topology files are in the "classic" format. 

## Purpose

This example seeks to demonstrate the generation and circularization of a linear dsDNA with sticky ends. 

We define a linear dsDNA with sticky ends as a linear dsDNA with two overhanging sequences on either ends (for a visual explanation, take a look at initial_struct.png). 

This example uses the external force called mutual traps to facilitate faster circularization (the mutual traps file is called ext.dat). To read more on external forces and mutual traps, head to the official documentation page: (https://lorenzo-rovigatti.github.io/oxDNA/forces.html)

Other than the seed and print_conf_interval options found in the input file, the rest of the simulation input options are the exact same as the simulation input options found in the CIR_ASSEMBLY example (note the lack of refresh_vel and log_file; they were deleted in this input file since it was not used in the CIR_ASSEMBLY input file).

*Note: Circularization did not occur even with mutual traps likely because the length of the dsDNA with sticky ends (35 bp) is too short. 

## How to generate a linear dsDNA with sticky ends

First create a .txt file. In this example, we chose to name the file sticky_sequence.txt.

The format of sticky sequences is:

    STICKY 4 TTTTGGGGGGGGGGGGGGGGGGGGAGCTGCAGTCC AAAA

Format: STICKY length_of_sticky_overhang complete_first_strand_seq sec_overhang_seq

The first strand in the above example is TTTTGGGGGGGGGGGGGGGGGGGGAGCTGCAGTCC.

The sticky end sequence of the first strand is TTTT.

The sticky end sequence of the second strand is AAAA.

The keyword STICKY denotes that the following sequence in the same line will be a sticky sequence.

The number after STICKY denotes the length of the sticky end.
    
In the above example, the sticky ends are TTTT and AAAA (they must match). In this example the length is 4.

It is important to note that the sticky ends do not have to be complementary.
    
The code will automatically generate the complementary sequences making up the second strand after the "sticky end" of the first strand. In the above example, the code will automatically generate the complementary sequence to GGGGGGGGGGGGGGGGGGGGAGCTGCAGTCC. The complementary sequence makes up the second strand. The user should not type the complementary sequence that make up the second strand.

## Generate and Run Instructions

To generate your own intial configuration, the initial positions, orientation, and etc of each nucleotide in the model, write the following in the terminal:

```../../utils/generate-sa.py 40. sequences```

*Note: Format: path/to/generate-sa.py side_length_of_model_box path/to/sequences_file
*Note: Make sure you are in the STICKY_ENDS_dsDNA_CIRCULARIZATION directory.
*Note: ../../utils/generate-sa.py simply means the path to the generate-sa.py file
*Note: the 40 specifies the side length of the box in the model.

Running the above command will generate generated.dat, generated.dat.pyidx, and generated.top files. 

Read more on those files here: https://lorenzo-rovigatti.github.io/oxDNA/configurations.html#configuration-file. 

To run the model using the provided input file, type the following into the command line:

```../../build/bin/oxDNA ../../examples/STICKY_ENDS_dsDNA_CIRCULARIZATION/input```

*Note: Format: path/to/oxDNA_executable path/to/input_file
*Note: Make sure the build path is correct
*Note: Make sure all dependencies are installed properly.

## Outputs

The simulation will output the following files:

trajectory.dat - contains the information of the positions, orientations, and etc of all nucleotides in the simulation. It prints every 1e3 time steps as specified in the input file by the option print_conf_interval = 1e3. 

last_conf.dat - This file contains the information of the positions, orientations, and etc of all nucleotides at the last printed time step of the simulation. 

energy.dat - contains the information of the energy of the overall system. It prints every 1e3 time steps as specified in the input file by the option print_energy_every = 1e3.

Since the original simulation was not run with a set seed, the output files ran by the user will be different than the output files provided in the example directory. 

## Visualization

visualization.ipynb is a file that aids in visualizing the simulation. To see the initial configuration of the system. Simply run the first block of code (press the triangle button on the left). It will generate an interactive visualization in the white block below. Double click on a nucleotide to center it to the screen. 

*Note: Make sure generated.dat and generated.top exists before running the visualization code.

The second block of code visualizes the final configuration of the system. Run it in the same way as the first block of code. 

*Note: Make sure generated.top and last_conf.dat exists in this example directory.

initial_struct.png and final_png are in the examples directory for viewing pleasure. 

