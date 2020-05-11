# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin

import os

##########CONFIG###########
# 動画の保存形式を選択
GIF = 0
MP4 = 1
PLT = 0

TEST = 1
###########################

################### PARAMETER ##################
df_time = pd.read_csv("../data/output_time.csv")
T = df_time["time"]
FILE_PATH = '../figs/animation'
kappa = 1.0
################### PARAMETER ##################

# 解析解
def analytic(t, x):
  return exp(-9.0*kappa*pi*pi*t)*sin(3.0*pi*x) - exp(-1.0*kappa*pi*pi*t)*sin(1.0*pi*x)

#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_u(frame):
  fp = "../data/u/%.3f.csv" % T[frame]
  df_u = pd.read_csv(fp)
  x = df_u["x"]
  u = df_u["u"]
  return x,u

# 初期画像を設定
def init_u(ax):
  ax.set_ylabel("u", fontsize=21)
  ax.set_xlabel("x", fontsize=21)
  ax.set_title("1D Diffusion", fontsize=24)
  ax.grid(linestyle="dotted")
  ax.xaxis.set_tick_params(direction='in')
  ax.yaxis.set_tick_params(direction='in')
  ax.tick_params(labelsize=21)
  # ax.set_yscale("log")
  x,u = get_u(0)
  if TEST == 0:
    im_u, = ax.plot(x,u, "ro-", label="numeric")
    ax.legend(fontsize=20)
    return im_u
  else:
    im_u_analytic, = ax.plot(x,analytic(T[0],x), "bo-", label="analytic")
    im_u, = ax.plot(x,u, "ro-", label="numeric")
    ax.legend(fontsize=20)
    return im_u, im_u_analytic

# データ更新
def reset_u(im,frame):
  x,u = get_u(frame)
  im.set_data(x,u)

def reset_u_analytic(im,frame):
  x,_ = get_u(frame)
  im.set_data(x,analytic(T[frame],x))

#########################描画のための関数#########################



################################################################
########################### main ###############################
################################################################

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
fig.subplots_adjust(left=0.15)
ax_u = fig.add_subplot(111)
time_text = fig.text(0.01, 0.99, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')
# kappa 出力 後でやる
fig.text(0, 0.01, r"$\kappa$="+str(kappa),
          backgroundcolor="black",color="white", size=20)

#### アニメの初期画像生成
if TEST == 0:
  im_u = init_u(ax_u)
else:
  im_u, im_u_analytic = init_u(ax_u)
time_text.set_text("time = 0.000")

#### 画像更新用関数
def animate(frame):
  if(frame == 0):
    return
  print("frame:",frame)
  reset_u(im_u, frame)
  if(TEST):
    reset_u_analytic(im_u_analytic, frame)
  
  time_text.set_text("time = %.3f"%T[frame])

ani = FuncAnimation(fig, animate, frames=int(len(T))
              , interval=200, repeat=True, blit=False)


if(GIF == 1):
    ani.save(FILE_PATH+".gif", writer='pillow')
if(MP4 == 1):
    ani.save(FILE_PATH+".mp4", writer="ffmpeg", fps=5)
if(PLT == 1):
    plt.show()

