#!/usr/bin/env python3
import sys

file = sys.argv[1]
nwin = int(sys.argv[2])

coverage_dict = {i: 0.0 for i in range(nwin)}

with open(file, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        bin_idx = int(line[0].split('_')[2])  # Extract bin number
        cov = float(line[4])                  # Extract coverage value
        coverage_dict[bin_idx] += cov

noise_bins = list(range(10)) + list(range(nwin - 10, nwin))
noise = sum(coverage_dict[i] for i in noise_bins)
avg_noise = noise / 20 if total_noise != 0 else 1

for i in range(nwin):
    signal = coverage_dict[i]
    norm_signal = signal / avg_noise  # Normalize signal
    print(f'bin_{i}\t{round(norm_signal, 5)}')

