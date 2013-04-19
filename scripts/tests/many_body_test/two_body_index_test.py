
from __future__ import division
from imports import *
from nhqm.bases import many_body as mb
from nhqm.bases import many_body_spill as mbs
import scipy as sp
from scipy import linalg
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
import itertools
from nhqm.bases import mom_space as mom
from nhqm.problems import He5
from nhqm.calculations import QM as calc
from time import time

from nhqm.bases.fermion_state import FermionState


bra = FermionState([1,2])
ket = FermionState([2,3])
print mb.two_body_indexes(bra,ket)