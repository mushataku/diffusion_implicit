# -*- coding: utf-8 -*-

# フォークとパラメータ a_v の温度依存性

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import sqrt, pi
from numpy import exp, pi, sin, cosh, abs
import os
from matplotlib.ticker import ScalarFormatter

##########CONFIG###########
# 保存形式を選択
PLT = 0
FILE_PATH = '../figs/av'

def av(T):
  return 4.7e-4*(T/1e4)**(-0.5)

# 初期画像を設定
def set_ax(ax):
  ax.set_ylabel(r"$a_v$", fontsize=21)
  ax.set_xlabel("T [K]", fontsize=21)
  # ax.set_title("phi", fontsize=24)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  ax.tick_params(labelsize=21)
  # ax.set_yscale("log")
  ax.set_xscale("log")
  ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
  ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))


T = np.linspace(1,1e5,int(1e5)+1)

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.15)
ax = fig.add_subplot(111)
set_ax(ax)

ax.plot(T,av(T))
if(PLT):
  plt.show()
else:
  plt.savefig(FILE_PATH)

xp = (av(10.0)*1e3/sqrt(pi))**(1.0/3.0)
print(xp)