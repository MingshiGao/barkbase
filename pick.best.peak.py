#!/usr/bin/env python3
import sys

with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip().split("\t")
        signal = [float(i) for i in line[4].split(",")]
        peaks = line[3].split(",")
        print(peaks[signal.index(max(signal))])
