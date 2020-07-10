#include <cmath>
#include <cstdio>
#include <vector>
#include <time.h>
#include <sys/stat.h>
using vd = std::vector<double>;
using vvd = std::vector<vd>;
using vvvd = std::vector<vvd>;
using vvvvd = std::vector<vvvd>;


/**********************config*************************/
// terminal に出力して確認 
const int OUTPUT_TERMINAL = 0;
// 初期条件 -> 0:長方形, 1:三角形 2:sinx 3:無
const int INITIAL = 3;
// 0:Euler法, 1:Cranck-Nicolson法
const int METHOD = 0;
// 0:無, 1:delta, 2:一様
const int SOURCE = 2;
/*****************************************************/

/******************計算条件********************/

const double Lx = 1.0;
const int NX = 100 + 1;
const double dx = Lx/(NX-1);
const double kappa = 1.0;
const double OMEGA = 1.8;
const double dt = 0.1;
const double TMAX = 4.0 + 1e-9;
const double lx = kappa*dt/dx/dx;
const double EPS = 1e-9;

/*
以下出力枚数の調整用
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に配列 (DT, 2DT, ..., ENDTIME) をsetする
*/
const double T_EPS = 1.0e-5;
const double DT = 0.01;
const double ENDTIME = 2.0;
vd TIME;
//出力時刻を set する関数
void TIME_set();
/*********************************************/

vd Source(NX);
void Source_set();
void Source_reset();
/*************************** diffusion を解く関数 ***************************/
// void init(vvd &u);
void init(vd &u);
void boundary(vd &nu);

int diffusion(vd &u);
double SOR(vd &u, vd &a);
/*************************** diffusion を解く関数 ***************************/

/*********************************ファイル関連*********************************/
FILE *condition_fp = fopen("../data/condition.csv", "w");
FILE *time_fp = fopen("../data/output_time.csv", "w");
FILE *u_fp;
char u_filename[50];
// その時刻における u の値をファイルにアウトプット
void output(double t, vd &u);
/*****************************************************************************/

int main(){
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);

  vd u(NX);
  int ti = 0, ite;
  double du, nu;
  double t = 0.0;
  TIME_set();
  Source_set();
  init(u);
  fprintf(time_fp, "time\n");
  output(t, u);
  printf("lx:%f\n", lx);
  printf("****************CALICULATION START****************\n");
  bool reset = 1;
  while(t < TMAX) {
    // u から nu を求める
    ite = diffusion(u);
    if(ite == -1){
      printf("Diffusion eqeation does not converge!!!\n");
      return 0;
    }
    t += dt;
    printf("iteration:%d when time=%f\n", ite, t);
    if(t > 2.0 && reset){
      Source_reset();
      printf("reset!!\n");
      reset = 0;
    }
    // nu が求め終わったのでこれを出力
    if(ti < TIME.size() && t > TIME[ti] - T_EPS){
      output(t, u);
      ti++;
    }
  }

  fclose(time_fp);
  fclose(condition_fp);
  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);
  return 0;
}



void init(vd &u){
  if(INITIAL == 0){
    for(int i = 0; i < NX; i++) {
      if(0.3*NX < i && i < 0.7*NX){
        u[i] = 1.0;
      }
      else{
        u[i] = 0.0;
      }
    }
  }
  
  if(INITIAL == 1){
    for(int i = 0; i < NX; i++) {
      u[i] = 0.5*Lx-std::abs(i*dx-0.5*Lx);
    }
  }

  if(INITIAL == 2){
    for(int i = 0; i < NX; i++) {
      u[i] = std::sin(3.0*M_PI*dx*i) - std::sin(1.0*M_PI*dx*i);
    }
  }

  if(INITIAL == 3){
    for(int i = 0; i < NX; i++) {
      u[i] = 0;
    }
  }
  boundary(u);
}

void TIME_set(){
  double t = DT;
  while(t < ENDTIME + T_EPS){
    TIME.push_back(t);
    t += DT;
  }
}

void output(double t, vd &u){
  sprintf(u_filename, "../data/u/%.3f.csv", t);
  FILE *fp = fopen(u_filename, "w");
  fprintf(time_fp, "%f\n", t);

  fprintf(fp, "x,u\n");
  if(OUTPUT_TERMINAL) printf("time:%f\n", t);
  
  for(int i = 0; i < NX; i++) {
    if(OUTPUT_TERMINAL){
      if(i < NX-1) printf("%f ", u[i]);
      else printf("%f\n", u[i]);
    }
    fprintf(fp, "%e,%e\n", i*dx, u[i]);
  }
  fclose(fp);
}

void boundary(vd &nu){
  nu[0] = nu[NX-1] = 0.0;
}

int diffusion(vd &u){
  // u から nu を求める
  vd a(NX);
  int imax = 99999;
  double du;
  for(int i = 0; i < NX; i++) {
    if(METHOD == 0){
      // u[i] += dt*Source[i];
      // a[i] = u[i];
      a[i] = u[i] + dt*Source[i];
    }
    if(METHOD == 1){
      a[i] = u[i] + lx/2.0*(u[i+1]+u[i-1]-2.0*u[i]);
    }
  }
  // nu を求める
  for(int icnt = 0; icnt < imax; icnt++) {
    // u から 仮のnuを求める
    // |u-nd| < EPS ならこの u をアクセプト
    du = SOR(u, a);
    if(du < EPS) return icnt;
  }

  // imax 回反復しても収束しないなら -1 を返して終了
  return -1;
}

double SOR(vd &u, vd &a){
  double nu, du = 0.0;
  // 仮の u から nu (仮のnu)を求める 一番外側 i==NX-1 は更新しない
  for(int i = 1; i < NX-1; i++) {
    if(METHOD == 0){
      nu = OMEGA*(lx*(u[i+1]+u[i-1])+a[i])/(1.0+2.0*lx);
    }
    if(METHOD == 1){
      nu = OMEGA*( lx/2.0/(1.0+lx)*(u[i+1]+u[i-1]) + a[i]/(1.0+lx) );
    }
    nu += (1.0-OMEGA)*u[i];
    du = std::max(std::abs(nu-u[i]), du);
    u[i] = nu;
  }
  {// i = NX-1
    int i = NX-1;
    if(METHOD == 0 && SOURCE == 2){
      double phi = 0.1;
      // nu = OMEGA*(lx*u[i-1]+a[i])/(1.0+1.0*lx);
      nu = OMEGA*(lx*u[i-1]+a[i])/(1.0+(1.0+phi*dx)*lx);
    }
    nu += (1.0-OMEGA)*u[i];
    du = std::max(std::abs(nu-u[i]), du);
    u[i] = nu;
  }
  return du;
}

void Source_set(){
  if(SOURCE == 1){
    for(int i = 0; i < NX; i++) {
      if(i == NX/2) Source[i] = 1.0/dx;
      else Source[i] = 0.0;
    }
  }
  if(SOURCE == 2){
    for(int i = 0; i < NX; i++) {
      Source[i] = 1.0;
    }
  }
}

void Source_reset(){
  for(int i = 0; i < NX; i++) {
    Source[i] = 0.0;
  }
}
