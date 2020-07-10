# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sqrt, pi
from numpy import exp, pi, sin, cosh, abs
import os
from matplotlib.ticker import ScalarFormatter
from scipy import integrate
##########CONFIG###########
# 動画の保存形式を選択
GIF = 0
MP4 = 0
PLT = 0

TEST = 0
FILE_PATH = '../figs/analytic'

LOG = 0
YMAX = 1e-1
YMIN = 1e-15
XMAX = 100.0
a = 1.5e-2
# a = sqrt(pi)
# Rmax = 1e4

def analytic(x, Rmax):
  tmp1 = sqrt(pi/24.0)/(a*Rmax)*x*x
  tmp2 = sqrt(2.0*pi*pi*pi/27.0)*abs(x*x*x)/(a*Rmax)
  return tmp1/(1.0 + cosh(tmp2))

# 初期画像を設定
def set_ax(ax):
  ax.set_ylabel("J", fontsize=21)
  ax.set_xlabel("x", fontsize=21)
  ax.set_title("Analytic", fontsize=24)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  ax.tick_params(labelsize=21)
  ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
  ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
  if(LOG):
    ax.set_yscale("log")
    ax.set_ylim(YMIN,YMAX)

x = np.linspace(-XMAX,XMAX,1000+1)

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.15)
ax = fig.add_subplot(111)
set_ax(ax)
Rmax = 1e4
for _ in range(4):
  ax.plot(x,analytic(x,Rmax), label = "Rmax=%.0e"%Rmax)
  Rmax *= 10.0
ax.legend()
if(PLT):
  plt.show()
else:
  plt.savefig(FILE_PATH)


Rmax = 1.0
for i in range(10):
  print("i:",i)
  I, _ = integrate.quad(lambda x: analytic(x,Rmax), -np.inf,np.inf)
  print(I)
  Rmax *= 10.0

  