import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numpy import exp, pi, sin
from mpl_toolkits.mplot3d import Axes3D
import os

##########CONFIG###########
# 動画の保存形式を選択
GIF = 0
MP4 = 0
PLT = 1
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
FILE_PATH = '../figs/3D_animation'
TITLE = "2D diffusion"
################### PARAMETER ##################

# 解析解
def analytic(x, y, t):
  return exp(-2.0*kappa*pi*pi*t)*sin(pi*x)*sin(pi*y)

#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_u(frame, U):
  fp = "../data/u/%.3f.csv" % T[frame]
  df_u = pd.read_csv(fp)
  u = df_u["u"]
  for i in range(NX*NY):
    U[i//NY][i%NY] = u[i]

def SET_AX(ax):
  ax.set_xlabel("x", fontsize=21)
  ax.set_ylabel("y", fontsize=21)
  ax.set_zlabel("u", fontsize=21)
  ax.set_title(TITLE, fontsize=24)
  ax.tick_params(labelsize=21)
  ax.set_zlim(0,1.0)

# 初期画像を設定
def init_u(ax):
  SET_AX(ax)
  get_u(0, U)
  im_u = [ax.plot_surface(X,Y,U,cmap="bwr")]
  fig.colorbar(im_u[0],shrink=0.9,ticks=[0.0,0.2,0.4,0.6,0.8,1.1])
  return im_u

#########################描画のための関数#########################




#### 画像更新用関数
def animate(frame, U, im_u):
  if(frame == 0):
    return
  print("frame:",frame)
  # reset_u(im_u, frame)
  get_u(frame, U)
  im_u[0].remove()
  im_u[0] = ax_u.plot_surface(X, Y, U, cmap="bwr")
  
  time_text.set_text("time = %.3f"%T[frame])

################################################################
########################### main ###############################
################################################################

# 基本的な部品を宣言  
fig = plt.figure(figsize=(8, 8))
# fig.subplots_adjust(left=0.80)
ax_u = fig.add_subplot(111, projection="3d")
# ax_u = Axes3D(fig)
# im_u = init_u(ax_u)

time_text = fig.text(0.01, 0.99, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')
# kappa 出力。後でやる
# fig.text(0, 0.01, "Re="+str(Re),
#           backgroundcolor="black",color="white", size=20)

SET_AX(ax_u)
get_u(0, U)
im_u = [ax_u.plot_surface(X,Y,U,cmap="bwr")]
# fig.colorbar(im_u[0],shrink=0.9,ticks=[0.0,0.2,0.4,0.6,0.8,1.1])

#### アニメの初期画像生成
ani = FuncAnimation(fig, animate, frames=int(len(T)), fargs=(U,im_u)
              , interval=20, repeat=True, blit=False)



if(GIF == 1):
    ani.save(FILE_PATH+".gif", writer='pillow', fps=25)
if(MP4 == 1):
    ani.save(FILE_PATH+".mp4", writer="ffmpeg", fps=25)
if(PLT == 1):
    plt.show()

