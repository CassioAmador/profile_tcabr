import sys
import matplotlib.pylab as plt
sys.path.insert(0, './../src/')

from proc_sweep import ProcSweep

shot = ProcSweep(28268)
shot.read_single_sweep(channel='Ka', sweep_cur=1000)
sig = shot.signal_filter(channel='Ka', freqs=(2e3, 16e3))

plt.figure(1)
matrix_ka = shot.spectrogram(channel='Ka', figure=1000, normal=1, freqs=(2e3, 16e3))
plt.pcolormesh(shot.X['Ka'], shot.Y['Ka'], matrix_ka)

plt.figure(2)
matrix_ka = shot.spectrogram2(channel='Ka', figure=1000, normal=1, freqs=(2e3, 16e3))
plt.pcolormesh(shot.X['Ka'], shot.Y['Ka'], matrix_ka)
plt.show()

