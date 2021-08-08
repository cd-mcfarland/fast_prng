#!/usr/bin/env python3
#
# Plots histogram data in `counts.csv` created by histogram.c
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("counts.csv")
x = df.pop('Low')
w = (df.pop("High") - x).iloc[0]
S = df[df.columns[0]]
plt.bar(x, S.values, w, align='edge')
ax = plt.gca()
ax.set(xlabel='x', title=S.name, ylabel='Counts')
plt.savefig('histogram.png', dpi=300)
