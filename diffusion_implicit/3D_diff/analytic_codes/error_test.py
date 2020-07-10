# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin
from matplotlib.ticker import ScalarFormatter

import os

##########CONFIG###########
# 動画の保存形式を選択
GIF = 0
MP4 = 1
PLT = 0

TEST = 1
###########################

################### PARAMETER ##################
df = pd.read_csv("../data/error.csv")
T = df["time"]
error = df["error"]
FILE_PATH = '../figs/error_test.png'
kappa = 1.0
################### PARAMETER ##################

################################################################
########################### main ###############################
################################################################

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
# fig.subplots_adjust(left=0.15)
ax = fig.add_subplot(111)
# kappa 出力 後でやる
fig.text(0, 0.01, r"$\kappa$="+str(kappa),
          backgroundcolor="black",color="white", size=20)

ax.plot(T, error, "o-")
ax.set_ylabel("error", fontsize=21)
ax.set_xlabel("t", fontsize=21)
ax.set_title("error of 3D Diffusion", fontsize=24)
ax.grid(linestyle="dotted")
ax.xaxis.set_tick_params(direction='in')
ax.yaxis.set_tick_params(direction='in')
ax.tick_params(labelsize=21)
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
ax.yaxis.offsetText.set_fontsize(20)

# plt.show()
plt.savefig(FILE_PATH)