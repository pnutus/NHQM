from scipy.fftpack import fft

def bfft(samples, interval):
    """Best-Friends-Forever Transform!!!"""
    N = len(samples)
    sampfreq = N / float(interval)
    transform = fft(samples)[:(N // 2)] / N * 2
    frequency = sp.linspace(0, sampfreq / 2., N // 2) * (2 * sp.pi)
    return transform, frequency

if __name__ == '__main__':
    x = sp.linspace(0, 5, 1024)
    T, F = bfft(sp.sin(10*x), 5)
    plt.plot(F, abs(T)**2)
    plt.axes()
    plt.show()