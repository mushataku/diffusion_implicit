import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin
import mpl_toolkits.axes_grid1

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
cfp = open('../data/condition.txt')
NX, NY = map(int, cfp.readline().split())
kappa, Lx, Ly = map(float, cfp.readline().split())
cfp.close
X, Y = np.meshgrid(np.linspace(0,Lx,NX), np.linspace(0,Ly,NY))
_, U = np.meshgrid(np.linspace(0,Lx,NX), np.linspace(0,Ly,NY))
FILE_PATH = '../figs/animation'
TITLE = "2D diffusion"
################### PARAMETER ##################

#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_u(frame):
  fp = "../data/u/%.3f.csv" % T[frame]
  df_u = pd.read_csv(fp)
  u = df_u["u"]
  for i in range(NX*NY):
    U[i//NY][i%NY] = u[i]

# 初期画像を設定
def init_u(ax):
  ax.set_xlabel("x", fontsize=21)
  ax.set_ylabel("y", fontsize=21)
  ax.set_title(TITLE, fontsize=24)
  ax.tick_params(labelsize=21)

  # 参照： https://qiita.com/nishimuraatsushi/items/2df8542ebf97affa26fb
  divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
  cax = divider.append_axes('right', '5%', pad='3%')

  get_u(0)
  im = ax.pcolormesh(X, Y, U, cmap="jet")
  cbar = fig.colorbar(im, cax=cax, orientation="vertical") # カラーバーの表示
  cbar.ax.tick_params(labelsize=18)
  return im

# データ更新
def reset_u(im,frame):
  get_u(frame)
  # set_array の引数は(2Dを無理やり1Dにした)1D の配列
  # じゃあ u を渡せばよいか？というわけには行かず、なぜかこうするとうまくいく
  # 参照：https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data/31490420#31490420
  # 訳が分からん…
  im.set_array(U[:-1,:-1].ravel())


#########################描画のための関数#########################


################################################################
########################### main ###############################
################################################################

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
# fig.subplots_adjust(left=0.80)
ax_u = fig.add_subplot(111)
time_text = fig.text(0.01, 0.99, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')
# kappa 出力 後でやる
# fig.text(0, 0.01, "Re="+str(Re),
#           backgroundcolor="black",color="white", size=20)

#### アニメの初期画像生成
im_u = init_u(ax_u)


#### 画像更新用関数
def animate(frame):
  if(frame == 0):
    return
  print("frame:",frame)
  reset_u(im_u, frame)
  
  time_text.set_text("time = %.3f"%T[frame])

ani = FuncAnimation(fig, animate, frames=int(len(T))
              , interval=200, repeat=True, blit=False)


if(GIF == 1):
    ani.save(FILE_PATH+".gif", writer='pillow', fps=25)
if(MP4 == 1):
    ani.save(FILE_PATH+".mp4", writer="ffmpeg", fps=25)
if(PLT == 1):
    plt.show()

