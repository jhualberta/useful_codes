from scipy import signal
original = [0, 1, 0, 0, 1, 1, 0, 0]
impulse_response = [2, 1]
recorded = signal.convolve(impulse_response, original)
print recorded
recovered, remainder = signal.deconvolve(recorded, impulse_response)
print recovered
