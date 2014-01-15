# This file does a relative import of the nhqm library to make the scripts 
# simpler to run. However, we recommend that you add the nhqm folder to your
# $PYTHONPATH environment variable.

import scipy as sp
import matplotlib.pyplot as plt
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")
