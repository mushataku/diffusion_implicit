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

TEST = 0
###########################

################### PARAMETER ##################
df_time = pd.read_csv("../data/output_time.csv")
T = df_time["time"]
cfp = open('../data/condition.txt')
NX, NY = map(int, cfp.readline().split())
kappa, Lx, Ly = map(float, cfp.readline().split())
cfp.close
FILE_PATH = '../figs/3D_test'
TITLE = "2D diffusion"
################### PARAMETER ##################

# 解析解
def analytic(x, y, t):
  return exp(-2.0*kappa*pi*pi*t)*sin(pi*x)*sin(pi*y)

#########################描画のための関数#########################

# frame 枚目の写真を作るデータを取得
def get_u(frame):
  fp = "../data/u/%.3f.csv" % T[frame]
  df_u = pd.read_csv(fp)
  u = df_u["u"]
  for i in range(NX*NY):
    U[i//NY][i%NY] = u[i]


error=[]
for i in range(int(len(T))):
  get_u(i)
  tmp = np.linspace(0.0,0.0,1)
  for jx in range(NX):
    for jy in range(NY):
      x = Lx*jx/NX
      y = Ly*jy/NY
      du = abs(U[jx][jy]-analytic(x, y, T[i]))
      tmp[0] = max(du-tmp[0])
  error.append(tmp[0])


fig = plt.figure(8,8)
ax = fig.subplots(111)
ax.plot(T,error)

plt.show()