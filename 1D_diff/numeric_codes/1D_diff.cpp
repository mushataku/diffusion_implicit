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
// 初期条件 -> 0:長方形, 1:三角形 2:sinx 
const int INITIAL = 2;
/*****************************************************/

/******************計算条件********************/

const double Lx = 1.0;
const int NX = 100 + 1;
const double dx = Lx/(NX-1);
const double kappa = 1.0;
const double OMEGA = 1.8;
const double dt = 0.01;
const double TMAX = 2.0 + 1e-9;
const double lx = kappa*dt/dx/dx;
const double EPS = 1e-9;

/*
以下出力枚数の調整用
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に配列 (DT, 2DT, ..., ENDTIME) をsetする
*/
const double T_EPS = 1.0e-10;
const double DT = 0.01;
const double ENDTIME = 1.0;
vd TIME;
//出力時刻を set する関数
void TIME_set();
/*********************************************/


/*************************** diffusion を解く関数 ***************************/
// void init(vvd &u);
void init(vd &u);

void boundary(vd &nu);

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
  vd a(NX);
  int ti = 0;
  double du, nu;
  double t = 0.0;
  TIME_set();
  printf("%d\n", (int)TIME.size());
  init(u);
  fprintf(time_fp, "time\n");
  output(t, u);
  printf("lx:%f\n", lx);
  printf("****************CALICULATION START****************\n");
  while(t < TMAX) {
    // u から nu を求める
    for(int i = 1; i < NX-1; i++) {
      a[i] = u[i] + lx/2.0*(u[i+1]+u[i-1]-2.0*u[i]);
    }
    // nu を求める
    du = 1.0;
    while(du > EPS){
      // 仮の u[ti+1] から nu = (仮のu[ti+1]) を求める
      du = 0.0;
      for(int i = 1; i < NX-1; i++) {
        nu = OMEGA*( lx/2.0/(1.0+lx)*(u[i+1]+u[i-1]) + a[i]/(1.0+lx) );
        nu += (1.0-OMEGA)*u[i];

        du = std::max(std::abs(nu-u[i]), du);
        u[i] = nu;
      }
      // du < EPS ならこの u をアクセプト
    }
    t += dt;
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

void TIME_set(){
  double t = DT;
    printf("%f\n", t);
  while(t < ENDTIME + T_EPS){
    TIME.push_back(t);
    t += DT;
  }
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
      u[i] = std::sin(M_PI*dx*i);
    }
  }
  boundary(u);
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