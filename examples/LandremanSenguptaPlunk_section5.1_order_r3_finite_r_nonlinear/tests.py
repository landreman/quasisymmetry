#!/usr/bin/env python3

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main QUASISYMMETRY directory.

exec(open('../testsCommon.py').read())

numFailures = 0

outputFile = readOutputFile()
referenceFile = readReferenceFile()

absoluteTolerance = 1e-12
relativeTolerance = 1e-100
numFailures += compareToReference(referenceFile,outputFile,'RBC',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'RBS',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'ZBC',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'ZBS',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'R0',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'z0',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'curvature',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'torsion',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
#numFailures += compareToReference(referenceFile,outputFile,'Boozer_toroidal_angle',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'sigma',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'X1c',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Y1c',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Y1s',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'X1c_untwisted',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'X1s_untwisted',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Y1c_untwisted',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Y1s_untwisted',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'elongation',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'effective_nfp',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'rms_curvature',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'max_curvature',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'iota',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'max_elongation',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'axis_helicity',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'standard_deviation_of_R',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'standard_deviation_of_Z',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)

numFailures += compareToReference(referenceFile,outputFile,'B20',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'X20',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'X2s',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'X2c',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Y20',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Y2s',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Y2c',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Z20',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Z2s',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Z2c',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)

outputFile.close()
referenceFile.close()

print("numFailures:",numFailures)
exit(numFailures > 0)
