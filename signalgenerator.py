import math
import numpy as np

def complexSignal(f1, f2, a1, a2, data_points = 3000, dT = 0.01,noisy = True,
                  mean = 0, std = 10,separate_signals = False):
    
    if noisy:
        noise = np.random.normal(mean, std, size=data_points)
    else:
        noise = np.zeros(shape=(data_points))
    
    tremor1 = []
    tremor2 = []
    
    frequencies = np.zeros(shape=(2, data_points))

    t = 0

    for i in range(data_points):
        t += dT
        tremor1.append(a1*math.sin(2 * math.pi * f1 * t))
        tremor2.append(a2* math.cos(2 * math.pi * f2 * t))
        if a1 > a2:
            frequencies[0][i] = f1
            frequencies[1][i] = f2
        else:
            frequencies[0][i] = f2
            frequencies[1][i] = f1
    if separate_signals:
        return tremor1, tremor2, noise
    else:
        return np.array(tremor1) + np.array(tremor2) + noise, frequencies