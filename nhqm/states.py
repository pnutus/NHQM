
class SingleParticleStationaryState:
    """
    Stationary (or quasi-stationary) state for a single-particle problem. 
    These can be used as a basis for a many-particle problem. 
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
        