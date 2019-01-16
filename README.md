# epistasis_landscape_map
Python script to create an epistasis fitness landscape diagram for one or more evolved metabolic pathways.

Usage: $ python3 epistasisLandscapeGenerator.py filePath fitCol --pathSigs 'sig1' 'sig2'... --pathNames 'name1' 'name2'...

Ex: $ python3 epistasisLandscapeGenerator.py example_fitness_data.csv 6 --pathSigs '\*\*\*\*2' '\*\*\*\*1' '2\*\*\*1' --pathNames 'Foreign Pathway Alone' 'Both Pathways in Combination' 'Native Pathway Alone'

epistasisLandscapeGenerator.py takes a .csv file, a column number, a list of pathway genotype signatures, and a list of names for the pathways. It generates a window containing a tkinter canvas with a diagram of of the fitness landscapes for those pathways. 
A button at the bottom of the window saves a timestamped .png screenshot of the canvas to the current working directory.

Requires the modules 'pandas' and 'pyscreenshot' ('ImageGrab' from PIL if on Windows) to be installed prior to running.
The other modules should already be included with python3.

The argument 'filePath' is the relative or absolute path to the .csv file containing the data to plot.
The first .csv column should be an ID, the next n columns indicate the alleles, followed by the fitness data.
An example .csv is included in this repository; it will generate the example_landscape.png image likewise included if the script is run with the above example usage in a directory containing both the script and the example.csv file.

The argument 'fitCol' is the column number (0-indexed) of the .csv file containing the fitness data.

The argument 'pathSigs' is a list (space-delimited) of strings each of length n indicating a pathway genotype signature.

Use 0, 1, and 2 to represent ancestral, evolved, and absent versions of alleles, respectively.

Use these same numeric indicators in the .csv file. Use * to represent either 0 or 1 in the path signatures.

The argument 'pathNames' is a list of names for each pathway; must have the same number of elements as 'pathSigs'.

All four arguments are required.

The script can handle up to 10 genes and 6 pathways by default, though these can be increased arbitrarily by simply expanding the color palletes with standard tkinter-recognized color names. I would question the sanity of doing so.

Note: epistasisLandscapeGenerator.py should ideally not be run through an Anaconda python environment on Linux, as modern fonts referenced by tkinter will not be recognized by Anaconda's font libraries. It will work, it will just be ugly. 

An Anaconda rep. has said they have no intention of fixing this. Use the system's python3 instead. 

Check what python you are currently using with "$ python (or python3) --version". Anaconda python environments will identify themselves as such. Calling the full path to the system python (i.e. "$ /usr/bin/python(3)") will guarantee the use of system's python. Use "$ source deactivate" (perhaps more than once, depending on your environment setup) to close the current python environment to potentially move from an Anaconda environment to a system one. To access the system python without this explicit full path call, you may potentially need to alter your ~/.bashrc file, specifically the 'export PATH="/path/to/anaconda/bin:$PATH"' line, automatically appended by Anaconda during installation to adjust which version of python BASH sets a path to by default upon opening a terminal.
