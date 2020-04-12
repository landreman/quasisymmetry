#!/usr/bin/env python3

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main QUASISYMMETRY directory.

exec(open('../testsCommon.py').read())

numFailures = 0

outputFile = readOutputFile()
referenceFile = readReferenceFile()

absoluteTolerance = 1e-13
relativeTolerance = 1e-100
numFailures += compareToReference(referenceFile,outputFile,'iotas',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'max_elongations',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'rms_curvatures',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'max_curvatures',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'max_modBinv_sqrt_half_grad_B_colon_grad_Bs',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'axis_lengths',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'standard_deviations_of_R',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'standard_deviations_of_Z',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'axis_helicities',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'scan_eta_bar',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'eta_bar_values',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)

outputFile.close()
referenceFile.close()

print("numFailures:",numFailures)
exit(numFailures > 0)
