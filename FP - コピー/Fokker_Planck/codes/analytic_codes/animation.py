# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from math import sqrt, pi
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin, cosh, abs
from matplotlib.ticker import ScalarFormatter

import os

##########CONFIG###########
# 動画の保存形式を選択
GIF = 0
MP4 = 0
PLT = 1

TEST = 1
# 0:rJ 1:xJ
TYPE = 1

# 0:Full 1:途中までのデータでアニメ作成
TMP_PLOT = 0

a = 1.5e-2
Rmax = 1e5

X = "3.000"
R = str(Rmax) + "00"
# R = Rmax
# YMAX = 0
YMAX = 1e0
YMIN = 1e-11
LOG = 0
###########################

#--------------def-------------
def readcsv(datafilename):
  df=pd.read_csv(datafilename,header=None)
  dft = np.array(df).T
  dftr = dft[:,1:]
  dftru = np.array(dftr)
  return dftru

################### PARAMETER ##################
FILE_PATH = '../figs/animation'
if(TMP_PLOT):
  df_time = pd.read_csv("../data/output_time_tmp.csv")
else:
  df_time = pd.read_csv("../data/output_time.csv")
T = df_time["time"]
# df_JV = pd.read_csv("../data/conservation.csv")
# JV = df_JV["JV"]
# co = readcsv("../data/condition.csv")
################### PARAMETER ##################

# 解析解
def analytic(x):
  tmp1 = sqrt(pi/24.0)/a/Rmax*x*x
  tmp2 = sqrt(2.0*pi*pi*pi/27.0)*abs(x*x*x)/a/Rmax
  return tmp1/(1.0 + cosh(tmp2))

#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_xJ(frame):
  fp = "../data/J/xJ/r_" + R + "/%.3f.csv" % T[frame]
  df_xJ = pd.read_csv(fp)
  x = df_xJ["x"]
  J = df_xJ["J"]
  return x,J

def get_rJ(frame):
  fp = "../data/J/rJ/x_" + X + "/%.3f.csv" % T[frame]
  df_rJ = pd.read_csv(fp)
  r = df_rJ["r"]
  J = df_rJ["J"]
  return r,J

# 初期画像を設定
def init_J(ax, TYPE):
  ax.set_ylabel("J", fontsize=21)
  ax.set_title("FP eq", fontsize=24)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  ax.tick_params(labelsize=21)
  # ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
  # ax.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
  ax.yaxis.offsetText.set_fontsize(20)
  if(LOG):
    ax.set_yscale("log")
    ax.set_ylim(YMIN,YMAX)
  elif(YMAX != 0):
    ax.set_ylim(0,YMAX)
  if(TYPE == 0):
    ax.set_xlabel("r", fontsize=21)
    x,J = get_rJ(0)
  else:
    ax.set_xlabel("x", fontsize=21)
    x,J = get_xJ(0)
  if TEST == 0:
    im_J, = ax.plot(x,J, "ro", label="numeric")
    ax.legend(fontsize=20)
    return im_J
  else:
    im_J_analytic, = ax.plot(x,analytic(x), "b", label="analytic")
    im_J, = ax.plot(x,J, "r", label="numeric")
    ax.legend(fontsize=20)
    return im_J, im_J_analytic

# データ更新
def reset_J(im,frame,TYPE):
  if(TYPE == 1):
    x,J = get_xJ(frame)
    im.set_data(x,J)
  else:
    r,J = get_rJ(frame)
    im.set_data(r,J)
    

def reset_J_analytic(im,frame):
  x,_ = get_xJ(frame)
  im.set_data(x,analytic(x))


#########################描画のための関数#########################



################################################################
########################### main ###############################
################################################################

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.2)
ax_J = fig.add_subplot(111)
time_text = fig.text(0.01, 0.99, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')
# 位置の出力
if(TYPE):
  fig.text(0, 0.01, "at r = %e"%Rmax,
          backgroundcolor="black",color="white", size=20)
else:
  fig.text(0, 0.01, "photon of x = " + X,
          backgroundcolor="black",color="white", size=20)


#### アニメの初期画像生成
if TEST == 0:
  im_J = init_J(ax_J, TYPE)
else:
  im_J, im_J_analytic = init_J(ax_J, TYPE)
# time_text.set_text("time = 0.000\n"+r"$\iiint$JdV = %.3e"%JV[0])
time_text.set_text("time = 0.000")

#### 画像更新用関数
def animate(frame):
  if(frame == 0):
    return
  print("frame:",frame)
  reset_J(im_J, frame, TYPE)
  if(TEST):
    reset_J_analytic(im_J_analytic, frame)
  
  # time_text.set_text("time = %.3f\n"%T[frame]+r"$\iiint$JdV = %.3e"%JV[frame])
  time_text.set_text("time = %.3e"%T[frame])

ani = FuncAnimation(fig, animate, frames=int(len(T))
              , interval=200, repeat=True, blit=False)


if(GIF == 1):
    ani.save(FILE_PATH+".gif", writer='pillow')
if(MP4 == 1):
    ani.save(FILE_PATH+".mp4", writer="ffmpeg", fps=5)
if(PLT == 1):
    plt.show()

