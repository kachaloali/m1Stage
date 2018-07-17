# myError.py
# Jeanne Cambefort
# 18/03/14
# IRISA

import sys
from inspect import currentframe, getframeinfo

def printError(message):
    """
    It prints an error message, then it exits.
    """
    frameinfo = getframeinfo(currentframe())
    print(frameinfo.filename, frameinfo.function, frameinfo.lineno)
    sys.exit(message)
