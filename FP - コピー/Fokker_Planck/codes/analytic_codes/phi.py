# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sqrt, pi
from numpy import exp, pi, sin, cosh, abs
import os

##########CONFIG###########
# 保存形式を選択
PLT = 0
FILE_PATH = '../figs/phi_log'

a = 1.5e-2
# a = sqrt(pi)

def phi_ana(x):
  return a/(x*x)/sqrt(pi)
def phi_eps(x,eps):
  return a/(eps+x*x)/sqrt(pi)

# 初期画像を設定
def set_ax(ax):
  ax.set_ylabel(r"$\phi(x)$", fontsize=21)
  ax.set_xlabel("x", fontsize=21)
  # ax.set_title("phi", fontsize=24)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  ax.tick_params(labelsize=21)
  ax.set_yscale("log")
  ax.set_xscale("log")


x = np.linspace(0,5.0,int(1e3)+1)

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.15)
ax = fig.add_subplot(111)
set_ax(ax)

ax.plot(x,phi_ana(x),label=r"$a/\sqrt{\pi}x^2$")
eps = 1e-1
for _ in range(5):
  ax.plot(x,phi_eps(x,eps), label = r"$\varepsilon=$%.0e"%eps)
  eps *= 0.1
ax.legend()
if(PLT):
  plt.show()
else:
  plt.savefig(FILE_PATH)

