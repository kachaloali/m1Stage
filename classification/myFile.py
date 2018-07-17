# myFile.py
# Jeanne Cambefort
# 17/03/14
# IRISA

import os, re
from myError import *


"""
This module contains functions on files.
"""

def print2file(outputFile, string):
    """
    It just prints a string in a output file. 
    """
    output = open(outputFile, "w")
    if len(string) == 0:
        printError("This string is empty!")
    # It removes the last newline. 
    output.write(string.rstrip("\n"))
    output.close()

def concat(outputFile, string):
    """
    It just appends this string to the end of the outputFile.
    """
    if len(string) == 0:
        printError("This string is empty!")
    if not os.path.exists(outputFile):
        print2file(outputFile, string)
    else:
        output = open(outputFile, "a")
        # It removes the last newline. 
        output.write("\n" + string.rstrip("\n"))
        output.close()

def checkDirectory(directory):
    """
    It checks that the specified directory exist and it prints
    an error message if it does not exist.
    """
    if not os.path.exists(directory):
        printError("This " + directory + " directory does not exist!")


def createDirectory(directory):
    """
    It just creates a directory, if it doesn't exist.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def checkTblFile(tblFile):
    """
    It checks the specified tabulated file.
    """
    with open(tblFile, 'r') as tbl:
        for line in tbl:
            lineTab = line.split("\t")
            # In this case, there is no tabulation in the input file.
            if len(lineTab) <= 1:
                printError("There is no tabulation in this "+ tblFile  
                           + " inputfile!")
