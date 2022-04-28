## Philosophy

The new oxDNA suite test is being designed with the following ideas in mind

1. Implementing standard tests should be reasonably easy and the procedure
should be automated as much as possible
2. It should be possible to implement custom tests, although the library will
provide a set of default tests that should be able to fulfill all basic needs
3. It should be possible to design and run difference "classes" (here called
levels) of tests, ranging from very short simulations to scientific ones

You can have as many test levels as you want. As of now, you can run all the
"quick" tests from your build folder issuing the "make test_quick" command.

## How to implement a new test

Here is a brief tutorial to implement your own test:

1. Create a new folder (or a subfolder)
2. Add that folder to the test_folder_list.txt file (relative to the TEST_LR folder)
3. Prepare one input file for each level of testing you want to provide for that specific folder. These should be named as LEVEL_input (e.g. quick_input)
4. Test each input file: you should be able to run oxDNA LEVEL_input from within the folder for each input file
5. Prepare one compare file for each level of testing you want to provide for that specific folder. These should be named as LEVEL_compare (e.g. quick_compare). Compare files tell the test suite which tests should be performed (look at [the list of available tests](#available-tests). The layout is simple: one line per test. You can find an example at the end of this file.
6. Run the PrepareCompareFiles.py script like this:  ./PrepareCompareFiles.py YOUR_FOLDER PATH_TO_OXDNA LEVEL. For example, you created a MY_SYSTEM folder within the TEST_LR folder and you compile oxDNA out-of-source in the build folder (because you like good-practice policies!). If you wanted to setup a test of level 'quick' you would run the command "../PrepareCompareFiles.py . ../../build/bin/oxDNA quick" from within the MY_SYSTEM folder or "./PrepareCompareFiles.py MY_SYSTEM/ ../build/bin/oxDNA quick" from within the TEST_LR folder. This command will overwrite your LEVEL_compare file with the right reference values.
7. Run the test suite. Go to the TEST_LR folder and issue the command "./TestSuite.py test_folder_list.txt ../build/bin/oxDNA LEVEL". Note that the test suite will only run simulations and tests in folders having the LEVEL_compare and LEVEL_input files.

Here is a compare file example. All fields should be separated by a double colon. The
first field is the test class name. All the other fields are parameters that
will be passed to the class instance.

ColumnAverage::energy.dat::2::-2.09155866667::0.0237939725575
FileExists::energy.dat::True

## Available tests

* FileExists checks that a file exists. The test expects one mandatory and one optional field: the name of the file and whether it should be checked that the file should be non-empty (defaults to True).
* ColumnAverage computes the average over a column stored in a file and compares it with a reference value. The test fails if the computed value is outside a tolerance range. The test expects 4 fields: a file name, a column index, a reference value and the associated tolerance.
* DiffFiles checks whether two files are identical. The test expects 2 fields: the name of the reference file and the name of the data file that should be checked against it.