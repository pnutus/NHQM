from imports import *
import HO_basis

@sp.vectorize
def V(r, l, s):
    return - 1 / r

print HO_basis.energies(30, V)

# Save to file
# sp.savetxt("matris3.txt", H)

