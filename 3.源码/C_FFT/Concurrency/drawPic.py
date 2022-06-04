import numpy as np
import matplotlib.pyplot as plt

# 原始信号

ORI = []
with open("ORI.txt") as f:
    for line in f:
        ORI.extend([float(i) for i in line.split()])
# End 
ORI_X = [x for x in range(0, len(ORI))]

# 傅立叶变换

FFT = []
with open("FFT.txt") as f:
    for line in f:
        FFT.extend([float(i) for i in line.split()])
# End

plt.subplot(2, 1, 1)
plt.plot(ORI_X, ORI)
plt.subplot(2, 1, 2)
# plt.plot([i for i in range(0, 100)], [FFT[i] for i in range(0, 100)])
plt.plot(ORI_X, FFT)

plt.show()