/* 最終更新：2020/02/18
[機能]
・2D Navier-Stokes 方程式を非圧縮条件のもとで Fractional Step 法で計算する。
・50 行目あたりのパラメータをいじって計算条件を設定できる
・67 行目あたりで出力する時間の間隔、長さを調整できる
・INITIAL の値を変えることで初期条件を選べる
・BOUNDARY で境界条件を選択できるが fixed は試作中
・IMCOMPRESSIBLE を 0 にすると 連続の式を切って比較することができる

[出力形式]
・相対パス "./data/condition.txt" 中に以下の形式で計算条件を出力
NX NY FRAMES
Lx Ly Re mu
ただし、NX, NY は出力する格子点の数で FRAMES は出力した回数

・"./data/u.txt" に以下の形式で u の計算結果を出力
time1
（時刻 time1 での配列 u ）
time2
（時刻 time2 での配列 u ）
（以下略）

・v, p, rot, div についても同様
*/

#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <time.h>
#include <vector>

using vd = std::vector<double>;
using vvd = std::vector<vd>;


/******************************CONFIG******************************/
// 0:diagonal 1:random 2:left 3:sin 4:static
const int INITIAL = 0;

// 0:fixed, 1:periodic (fixedは試作中)
const int BOUNDARY = 1;

// 0:2D burgers, 1: imcompressible 
const int IMCOMPRESSIBLE = 1;
/*******************************************************************/


/******************************計算条件******************************/
// nx, ny は100くらいなら1分程度で終わる 200 くらいが精度良さそうだけど7分くらいかかる
const int nx = 100+5;
const int ny = 100+5;
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx/double(nx-5);
const double dy = Ly/double(ny-5);

const double Re = 10000.0;
//計算の安定性を決めるファクターμ, mu > 0.25 だと計算が爆発する
const double mu = 0.20;

/*
以下出力枚数の調整用
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に配列 (DT, 2DT, ..., ENDTIME) をsetする
*/
const double T_EPS = 1.0e-10;
const double DT = 0.02;
const double ENDTIME = 10.0;
vd TIME;
//出力時刻を set する関数
void TIME_set();
/*************************************************************************/

/*************************poisson eq 関連の関数・定数*************************/
const double OMEGA = 1.8; // sor法の加速係数
const double P_EPS = 1.0e-4;

// 仮の流速 u*, v* からpoisson 方程式の右辺 s を求める
void get_s(vvd &u, vvd &v, vvd &s, double dt);
// SOR法で1ステップだけf更新。誤差|f - fn|の最大値を返す
double sor(vvd &f, vvd &s);
// sor をたくさん呼び出して poisson 方程式を解く。返り値は反復回数
int poisson2d(vvd &f, vvd &s);
// 残差の最大値を計算。sor反復中に呼び出すことで精度の上昇を見れる。
double residual(vvd &f, vvd &s);

// 圧力から速度を補正
void correction(vvd &u, vvd &v, vvd &p, double dt);

/****************************************************************************/

/******************************burgersを解く関数******************************/
// 初期状態を決定
void initial(vvd &u, vvd &v);
// fn に境界条件を課す
void boundary(vvd &un, vvd &vn);

void diffusion(vvd &f, vvd &fn, double dt);
void x_advection(vvd &f, vvd &fn, vvd &u, double dt);
void y_advection(vvd &f, vvd &fn, vvd &v, double dt);
void rotation(vvd &u, vvd &v, vvd &rot);
void divergence(vvd &u, vvd &v, vvd &div);
/****************************************************************************/


/*********************************ファイル関連*********************************/
FILE *condition_fp = fopen("data/condition.txt","w");
FILE *u_fp = fopen("data/u.txt", "w");
FILE *v_fp = fopen("data/v.txt", "w");
FILE *div_fp = fopen("data/div.txt", "w");
FILE *rot_fp = fopen("data/rot.txt", "w");
FILE *p_fp = fopen("data/p.txt", "w");

// その時刻における f[jy][jx] の値をファイルにアウトプット
void output(vvd &f, double t, FILE *data_fp);
// u, v, p, rot, div をまとめてoutput
void output_all(vvd &u, vvd &v, vvd &p, vvd &rot, vvd &div, double t);
/*****************************************************************************/

int main(){
  double dt, t = 0.0;
  int ti = 0, icnt = 0; //TIMEのindex
  int iteration;

  vvd u(ny,vd(nx, 0.0)), un(ny, vd(nx, 0.0));
  vvd v(ny,vd(nx, 0.0)), vn(ny, vd(nx, 0.0));
  vvd rot(ny,vd(nx, 0.0)), div(ny, vd(nx, 0.0));
  vvd p(ny,vd(nx, 0.0)), s(ny, vd(nx, 0.0));

  clock_t start_t, end_t;
  start_t = time(NULL);

  initial(u, v);
  boundary(u, v);
  TIME_set();

  dt = fmin(0.2*fmin(dx,dy)/1.0, mu*fmin(dx*dx, dy*dy)*Re);
  
  printf("NX:%d NY:%d\nRe:%f mu:%f\n", nx, ny, Re, mu);
  printf("dt:%.10f\n", dt);
  output_all(u, v, p, rot, div, 0.0); icnt++;

  do{
    
    //printf("time:%f\n", t);

    {// x方向に移流
      x_advection(u, un, u, dt);
      x_advection(v, vn, u, dt);
      boundary(un, vn);
      u = un; v = vn;
    }

    {// y方向に移流
      y_advection(u, un, v, dt);
      y_advection(v, vn, v, dt);
      boundary(un, vn);
      u = un; v = vn;
    }

    {// 拡散
      diffusion(u, un, dt);
      diffusion(v, vn, dt);
      boundary(un, vn);
      u = un; v = vn;
    }

    // poisson 方程式と結合して速度の修正（ここで非圧縮の条件が入る）
    if(IMCOMPRESSIBLE == 1){
      get_s(u, v, s, dt);
      iteration = poisson2d(p, s);
      if(iteration == -1){
        printf("poisson does not converged when time = %f\n", t);
        return 0;
      }
      correction(u, v, p, dt);
      boundary(u, v);
    }
    
    t += dt;

    // 全ての更新が終わったら出力 
    if(ti < TIME.size() && t > TIME[ti] - T_EPS){
      output_all(u, v, p, rot, div, t);
      ti++; icnt++;
    }

  } while (t < ENDTIME + dt);

  printf("number of pictures:%d\n", icnt);
  fprintf(condition_fp, "%d %d %d\n%f %f %f %f\n", nx-4, ny-4, icnt, Lx, Ly, Re, mu);

  fclose(condition_fp);
  fclose(u_fp);
  fclose(v_fp);
  fclose(div_fp);
  fclose(rot_fp);
  fclose(p_fp);

  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);

  return 0;
}

void initial(vvd &u, vvd &v){
  if(INITIAL == 0){
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        if(0.3*nx < jx && jx < 0.7*nx && 0.3*ny < jy && jy < 0.7*ny){
          u[jy][jx] = 1.0;
          v[jy][jx] = 1.0;
        }
      }
    }
  }
  if(INITIAL == 1){
    //全部の点の50%くらいをランダムに選ぶ
    for(int i = 0; i < nx*ny*0.5; i++) {
      u[rand()%ny][rand()%nx] = 1.0;
      v[rand()%ny][rand()%nx] = 1.0;
    }
  }

  if(INITIAL == 2){
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        if(0.3*ny < jy && jy < 0.7*ny && 0.3*nx < jx && jx < 0.7*nx){
          u[jy][jx] = -1.0;
          v[jy][jx] = 0.0;
        }
      }
    }
  }
  if(INITIAL == 3){
    double x, y, kx = 2.0*M_PI, ky = 2.0*M_PI;
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        x = dx*(double)(jx-2);
        y = dy*(double)(jy-2);

        u[jy][jx] = -cos(kx*x)*sin(ky*y)/kx;
        v[jy][jx] = sin(kx*x)*cos(ky*y)/ky;
        
        x += 0.3; y+= 0.7;
        u[jy][jx] = -0.6*cos(2.0*kx*x)*sin(2.0*ky*y)/kx;
        v[jy][jx] = 0.6*sin(2.0*kx*x)*cos(2.0*ky*y)/ky;
        
      }
    }
  }

  u[ny/2][nx/2] = 1.0;
  
  //check用
  // for(int jy = 0; jy < ny; jy++) {
  //   for(int jx = 0; jx < nx; jx++) {
  //     if(jx < nx-1){
  //       printf("%.1f ", f[jy][jx]);
  //     }
  //     else{
  //       printf("%.1f\n", f[jy][jx]);
  //     }
  //   }
  // }

  return;
}

void output(vvd &f, double t, FILE *data_fp){

  fprintf(data_fp,"%f\n", t);

  //端っこの境界条件のためのダミーの格子点は出力しない
  for(int jy = 2; jy <= ny-3; jy++) {
    for(int jx = 2; jx <= nx-3; jx++) {
      if(jx < nx-3){
        fprintf(data_fp, "%f ", f[jy][jx]);
      }
      else{
        fprintf(data_fp, "%f\n", f[jy][jx]);
      }
    }
  }
}

void output_all(vvd &u, vvd &v, vvd &p, vvd &rot, vvd &div, double t){
  divergence(u, v, div);
  rotation(u, v, rot);
  output(u, t, u_fp);
  output(v, t, v_fp);
  output(div, t, div_fp);
  output(rot, t, rot_fp);
  output(p, t, p_fp);
}

void diffusion(vvd &f, vvd &fn, double dt){
  //境界は更新しない
  for(int jy = 2; jy < ny-2; jy++) {
    for(int jx = 2; jx < nx-2; jx++) {
      fn[jy][jx] = f[jy][jx] + dt * 
      ( (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1])/dx/dx
         + (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx])/dy/dy
      )/Re;
    }
  }
  return;
}

void boundary(vvd &un, vvd &vn){
  if(BOUNDARY == 0){
    for(int jy=0 ; jy < ny; jy++){
      un[jy][nx-1] = un[jy][nx-3];// 流入境界を導入中
      un[jy][nx-2] = un[jy][nx-3];
      un[jy][0]    = 1.0;
      un[jy][1]    = 1.0;
      vn[jy][nx-1] = vn[jy][nx-3];
      vn[jy][nx-2] = vn[jy][nx-3];
      vn[jy][1]    = 0.0;
      vn[jy][0]    = 0.0;
    }

    for(int jx=0 ; jx < nx; jx++){
      un[ny-1][jx] = 1.0;
      un[ny-2][jx] = 1.0;
      un[1][jx]    = 1.0;
      un[0][jx]    = 1.0;
      vn[ny-1][jx] = 0.0;
      vn[ny-2][jx] = 0.0;
      vn[1][jx]    = 0.0;
      vn[0][jx]    = 0.0;
    }
  }
  if(BOUNDARY == 1){
    for(int jy=0 ; jy < ny; jy++){
      un[jy][nx-1] = un[jy][3];
      un[jy][nx-2] = un[jy][2];
      un[jy][1] = un[jy][nx-3];
      un[jy][0] = un[jy][nx-4];
      vn[jy][nx-1] = vn[jy][3];
      vn[jy][nx-2] = vn[jy][2];
      vn[jy][1] = vn[jy][nx-3];
      vn[jy][0] = vn[jy][nx-4];
    }

    for(int jx=0 ; jx < nx; jx++){
      un[ny-1][jx] = un[3][jx];
      un[ny-2][jx] = un[2][jx];
      un[1][jx] = un[ny-3][jx];
      un[0][jx] = un[ny-4][jx];
      vn[ny-1][jx] = vn[3][jx];
      vn[ny-2][jx] = vn[2][jx];
      vn[1][jx] = vn[ny-3][jx];
      vn[0][jx] = vn[ny-4][jx];
    }
  }
  //check 用
  // printf("boundary check\n");
  // for(int jy = 0; jy < ny; jy++) {
  //   for(int jx = 0; jx < nx; jx++) {
  //     if(jx < nx-1){
  //       printf("%.1f ", vn[jy][jx]);
  //     }
  //     else{
  //       printf("%.1f\n", vn[jy][jx]);
  //     }
  //   }
  // }
  return;
}

void x_advection(vvd &f, vvd &fn, vvd &u, double dt){

  double a,b,c,z;

  for (int jy = 2; jy < ny - 2; jy++){
    for (int jx = 2; jx < nx - 2; jx++){
      if (u[jy][jx] > 0.0){
        a = (f[jy][jx+1] - 3.0*f[jy][jx] + 3.0*f[jy][jx-1] - f[jy][jx-2]) / (6.0*dx*dx*dx);
        b = (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1]) / (2.0*dx*dx);
        c = (2.0*f[jy][jx+1] + 3.0*f[jy][jx] - 6.0*f[jy][jx-1] + f[jy][jx-2]) / (6.0*dx);
      }
      else{
        a = (f[jy][jx+2] - 3.0*f[jy][jx+1] + 3.0*f[jy][jx] - f[jy][jx-1]) / (6.0*dx*dx*dx);
        b = (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1]) / (2.0*dx*dx);
        c = (-f[jy][jx+2] + 6.0*f[jy][jx+1] - 3.0*f[jy][jx] - 2.0*f[jy][jx-1]) / (6.0*dx);
      }
      z = -u[jy][jx]*dt;
      fn[jy][jx] = a*z*z*z + b*z*z + c*z + f[jy][jx];
    }
  }
}

void y_advection(vvd &f, vvd &fn, vvd &v, double dt){

  double a,b,c,z;

  for (int jy = 2; jy < ny - 2; jy++){
    for (int jx = 2; jx < nx - 2; jx++){
      if (v[jy][jx] > 0.0){
        a = (f[jy+1][jx] - 3.0*f[jy][jx] + 3.0*f[jy-1][jx] - f[jy-2][jx]) / (6.0*dy*dy*dy);
        b = (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx]) / (2.0*dy*dy);
        c = (2.0*f[jy+1][jx]+ 3.0*f[jy][jx] - 6.0*f[jy-1][jx] + f[jy-2][jx]) / (6.0*dy);
      }
      else{
        a = (f[jy+2][jx] - 3.0*f[jy+1][jx] + 3.0*f[jy][jx] - f[jy-1][jx]) / (6.0*dy*dy*dy);
        b = (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx]) / (2.0*dy*dy);
        c = (-f[jy+2][jx] + 6.0*f[jy+1][jx] - 3.0*f[jy][jx] - 2.0*f[jy-1][jx]) / (6.0*dy);
      }
      z = -v[jy][jx]*dt;
      fn[jy][jx] = a*z*z*z + b*z*z + c*z + f[jy][jx];
    }
  }
}

void TIME_set(){
  double tmp = DT;
  while(tmp < ENDTIME + T_EPS){
    TIME.push_back(tmp);
    tmp += DT;
  }
}

void rotation(vvd &u, vvd &v, vvd &rot){
  for(int jy = 1; jy < ny-1; jy++) {
    for(int jx = 1; jx < nx-1; jx++) {
      rot[jy][jx] = 0.5*(v[jy+1][jx+1] - v[jy+1][jx]+ v[jy][jx+1] - v[jy][jx])/dx
                  - 0.5*(u[jy+1][jx+1] - u[jy][jx+1] + u[jy+1][jx] - u[jy][jx])/dy;
    }
  }
}

void divergence(vvd &u, vvd &v, vvd &div){
  for (int jy = 1; jy < ny-1; jy++){
    for (int jx = 1; jx < nx-1; jx++){
      div[jy][jx] = 0.5*(u[jy+1][jx+1] - u[jy+1][jx] + u[jy][jx+1] - u[jy][jx])/dx
                  + 0.5*(v[jy+1][jx+1] - v[jy][jx+1] + v[jy+1][jx] - v[jy][jx])/dy;
    }
  }
}

void get_s(vvd &u, vvd &v, vvd &s, double dt){
  for (int jy = 2; jy < ny-2; jy++){
    for (int jx = 2; jx < nx-2; jx++){
      s[jy][jx] = 0.5*(u[jy+1][jx+1] - u[jy+1][jx] + u[jy][jx+1] - u[jy][jx])/dx
                  + 0.5*(v[jy+1][jx+1] - v[jy][jx+1] + v[jy+1][jx] - v[jy][jx])/dy;
      s[jy][jx] /= dt;
    }
  }
}

double residual(vvd &f, vvd &s){
  double res, rmax = 0.0;
  //各格子点の(d^2f/dx^2 + d^2f/dy^2 - s)を計算
  for (int jy = 1; jy < ny - 1; jy++){
    for (int jx = 1; jx < nx - 1; jx++){
      res = (f[jy][jx + 1] - 2.0 * f[jy][jx] + f[jy][jx - 1]) / (dx * dx) 
          + (f[jy + 1][jx] - 2.0 * f[jy][jx] + f[jy - 1][jx]) / (dy * dy)
          - s[jy][jx];
      rmax = fmax(res, rmax);
    }
  }

  return rmax;
}

double sor(vvd &f, vvd &s){
  double fn, err = 0.0;
  for (int jy = 1; jy < ny - 1; jy++){
    for (int jx = 1; jx < nx - 1; jx++){
      fn = ((f[jy][jx + 1] + f[jy][jx - 1]) / (dx * dx) 
         + (f[jy + 1][jx] + f[jy - 1][jx]) / (dy * dy) - s[jy][jx]) 
         * 0.5 * dx * dx * dy * dy / (dx * dx + dy * dy);
      err = fmax(fabs(fn - f[jy][jx]), err);
      f[jy][jx] = (1.0 - OMEGA) * f[jy][jx] + OMEGA * fn;
    }
  }

  return err;
}

int poisson2d(vvd &f, vvd &s){
  int icnt = 0, imax = 99999;
  double err;

  // 反復
  while(icnt++ < imax){
    //fの更新と更新前後の差の最大値確認
    err = sor(f, s);
    //更新してもほとんど変わらないようなら終了
    if(P_EPS > err) return icnt;
  }

  // imax 回反復しても収束しないなら-1を返して終了
  return -1;
}

void correction(vvd &u, vvd &v, vvd &p, double dt){
  for (int j = 2; j < ny-2; j++){
    for (int i = 2; i < nx-2; i++){
      u[j][i] += -0.5*(p[j][i] - p[j][i-1] + p[j-1][i] - p[j-1][i-1])/dx*dt;
      v[j][i] += -0.5*(p[j][i] - p[j-1][i] + p[j][i-1] - p[j-1][i-1])/dy*dt;
    }
  }
}