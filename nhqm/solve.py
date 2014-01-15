from nhqm.helpers.matrix import eigensolve
from nhqm.helpers.quantum import normalize
from nhqm.states import SingleParticleStationaryState

def solve(problem, quantum_numbers, basis):
    """
    Given a spherically symmetric problem, relevant quantum numbers and a basis,
    expands the problem in the basis and returns a list of stationary states, 
    in order of ascending binding energy (the real part of the energy).
    """
    H = basis.hamiltonian(problem, quantum_numbers)
    energies, eigenvectors = eigensolve(H, hermitian=basis.hermitian)
    
    # because the Berggren is weird, we need to normalize the eigenvectors
    try:
        eigenvectors = map(basis.normalize, eigenvectors.T)
    except AttributeError:
        pass # the basis has no special normalization method
    
    def create_solution(energy, eigenvector):
        return SingleParticleStationaryState(energy, eigenvector, problem, 
                               quantum_numbers, basis)
    solutions = map(create_solution, energies, eigenvectors)
    return solutions
    
