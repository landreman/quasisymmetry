#!/usr/bin/env python3

# This script is designed to be called by "make test" in the parent directory,
# not to be called directly.  The reason is that there are several system-dependent
# variables used which are set in the makefiles.

import os
import subprocess
from sys import stdout

def verifyVariableExists(str):
    try:
        temp = os.environ[str]
    except:
        print("Error!  Variable "+str+" is not set.  This error may be caused by calling runExamples.py directly rather than by calling 'make test'.")
        raise

    return temp

retestStr = verifyVariableExists("QUASISYMMETRY_RETEST")

if retestStr=="yes":
    retest=True
elif retestStr=="no":
    retest=False
else:
    print("Error! QUASISYMMETRY_RETEST must be either 'yes' or 'no'. There is likely an error in the main makefile.")
    exit(1)

wereThereAnyErrors = False
examplesWithErrors = []

try:
    submitCommand = os.environ["QUASISYMMETRY_COMMAND_TO_SUBMIT_JOB"]
except:
    submitCommand = ""

print("QUASISYMMETRY_COMMAND_TO_SUBMIT_JOB=",submitCommand)

# Get a list of the subdirectories:
subdirectories = filter(os.path.isdir, os.listdir("."))

# From this list, keep only subdirectories that have both a "tests.py" file and a job.SFINCS_SYSTEM file:
examplesToRun = []
directoriesThatArentExamples = []
for subdirectory in subdirectories:
    print("Examining subdirectory "+subdirectory)
    if os.path.isfile(subdirectory+"/tests.py"):
        if os.path.isfile(subdirectory+"/quasisymmetry_in."+subdirectory):
            examplesToRun.append(subdirectory)
        else:
            print("WARNING: directory "+subdirectory+" contains a tests.py file but no quasisymmetry_in.XXX file of the same name as the directory.")
    else:
        directoriesThatArentExamples.append(subdirectory)

print()
if len(examplesToRun) == 0:
    print("Error: No subdirectories of examples/ found containing a tests.py and one quasisymmetry_in.XXX file.")
    exit(1)

if len(directoriesThatArentExamples)>0:
    print("The following subdirectories of /examples do not contain a tests.py file and so will be ignored:")
    for example in directoriesThatArentExamples:
        print("   " + example)
    print()

print("The following examples will be used as tests:")
for example in examplesToRun:
    print("   " + example)
print()


#if isABatchSystemUsed == "no":
if True:
    for subdirectory in examplesToRun:

        print()
        print("Preparing to check example: "+subdirectory)
        try:
            os.chdir(subdirectory)
        except:
            print("Error occurred when trying to change directory to "+subdirectory)
            raise

        print("Moved to working directory "+os.getcwd())

        if not retest:
            try:
                os.remove("quasisymmetry_out."+subdirectory+".nc")
            except:
                # If the .nc output file does not exist, there will be an exception, but we can safely ignore it.
                pass
            
            print("Launching QUASISYMMETRY...")

            inputFile = "quasisymmetry_in."+subdirectory
            submitCommand2 = submitCommand+" ../../quasisymmetry "+inputFile
            submitCommandSplit = submitCommand2.split()
            print("About to submit the following command: ",submitCommandSplit)
            # Flush everything printed to stdout so far:
            stdout.flush()

            try:
                # Next we launch QUASISYMMETRY.
                #subprocess.call(["srun","-n","1","-c","32","../../quasisymmetry",inputFile])
                #subprocess.call(["../../quasisymmetry",inputFile])
                subprocess.call(submitCommandSplit)
            except:
                print("An error occurred when attempting to launch QUASISYMMETRY.")
                raise

        print()
        print("QUASISYMMETRY execution complete. About to run tests on output.")
        stdout.flush()

        try:
            testResults = subprocess.call("./tests.py")
        except:
            print("An error occurred when attempting to run tests.py in the following directory:")
            print(os.getcwd)
            raise

        if testResults > 0:
            wereThereAnyErrors = True
            examplesWithErrors.append(subdirectory)

        # Step back one directory
        os.chdir("..")

    print("-----------------------------------------------")
    print("Done with tests.")
    print("Examples attempted:")
    for subdirectory in examplesToRun:
        print("  " + subdirectory)


print()
# Report whether any tests failed.
if wereThereAnyErrors:
    print("AT LEAST ONE TEST WAS FAILED.")
    print("Examples which failed:")
    for x in examplesWithErrors:
        print("   "+x)
else:
    print("ALL TESTS THAT WERE RUN WERE PASSED SUCCESSFULLY.")

print()
