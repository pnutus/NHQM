from nhqm.QM_helpers import eigensolve

def solve(problem, quantum_numbers, basis):
    """
    Given a spherically symmetric problem, relevant quantum numbers and a basis,
    expands the problem in the basis and returns a list of stationary states, 
    in order of ascending binding energy (the real part of the energy).
    """
    H = basis.hamiltonian(problem, quantum_numbers)
    energies, eigenvectors = eigensolve(H, hermitian=basis.hermitian)
    
    
    def create_solution(energy, eigenvector):
        return StationaryState(energy, eigenvector, problem, 
                               quantum_numbers, basis)
    solutions = map(create_solution, energies, eigenvectors.T)
    return solutions
    

class StationaryState:
    """
    Stationary (or quasi-stationary) state for a single-particle problem.  
    """
    def __init__(self, energy, eigenvector, problem, quantum_numbers, basis):
        self.energy          = energy
        self.eigenvector     = eigenvector
        self.problem         = problem
        self.quantum_numbers = quantum_numbers
        self.basis           = basis
        
    def wavefunction(self, r):
        """
        The radial wavefunction for plotting.
        """
        # The wavefunction is cached so that it only has to be generated once
        try:
            return self.cached_wavefunction(r)
        except AttributeError:
            self.cached_wavefunction = self.basis.gen_wavefunction(
                                            self.eigenvector, 
                                            self.quantum_numbers)
            return self.cached_wavefunction(r)
        