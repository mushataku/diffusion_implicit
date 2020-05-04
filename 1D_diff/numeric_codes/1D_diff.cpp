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
const int NX = 50 + 1;
const double dx = Lx/(NX-1);
const double kappa = 1.0;
const double OMEGA = 1.8;
const double dt = 0.001;
const double TMAX = 1.0 + 1e-9;
const int NT = 100+1;
const double lx = kappa*dt/dx/dx;
const double EPS = 1e-5;
/*********************************************/


/*************************** diffusion を解く関数 ***************************/
void init(vvd &u);
void boundary(vd &nu);

/*************************** diffusion を解く関数 ***************************/
/*********************************ファイル関連*********************************/
FILE *condition_fp = fopen("../data/condition.csv", "w");
FILE *time_fp = fopen("../data/output_time.csv", "w");
FILE *u_fp;
char u_filename[50];
// その時刻における u の値をファイルにアウトプット
void output(int ti, double t, vvd &u);
/*****************************************************************************/

int main(){
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);

  vvd u(NT,vd(NX));
  vd a(NX);
  int ti = 0;
  double du, nu;
  double t = 0.0;
  init(u);
  printf("****************CALICULATION START****************\n");
  fprintf(time_fp, "time\n");
  output(ti, t, u);
  printf("lx:%f\n", lx);
  for(ti = 0, t = 0.0; ti < NT-1 && t < TMAX; ti++) {
    // u[ti] から u[ti+1] を求める
    // まず仮に u[ti+1] を u[ti] とする
    u[ti+1] = u[ti];
    for(int i = 1; i < NX-1; i++) {
      a[i] = u[ti][i] + lx/2.0*(u[ti][i+1]+u[ti][i-1]-2.0*u[ti][i]);
    }
    // u[ti+1] を求める
    du = 1.0;
    while(du > EPS){
      // 仮の u[ti+1] から nu = (仮のu[ti+1]) を求める
      du = 0.0;
      for(int i = 1; i < NX-1; i++) {
        nu = OMEGA*( lx/2.0/(1.0+lx)*(u[ti+1][i+1]+u[ti+1][i-1]) + a[i]/(1.0+lx) );
        nu += (1.0-OMEGA)*u[ti+1][i];

        du = std::max(std::abs(nu-u[ti+1][i]), du);

        u[ti+1][i] = nu;
      }
      // du < EPS ならこの u[ti+1] をアクセプト
    }
    // u[ti+1] が求め終わったのでこれを出力
    t += dt;
    output(ti+1, t, u);
  }

  fclose(time_fp);
  fclose(condition_fp);
  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);
  return 0;
}

void init(vvd &u){
  if(INITIAL == 0){
    for(int i = 0; i < NX; i++) {
      if(0.3*NX < i && i < 0.7*NX){
        u[0][i] = 1.0;
      }
      else{
        u[0][i] = 0.0;
      }
    }
  }
  
  if(INITIAL == 1){
    for(int i = 0; i < NX; i++) {
      u[0][i] = 0.5*Lx-std::abs(i*dx-0.5*Lx);
    }
  }

  if(INITIAL == 2){
    for(int i = 0; i < NX; i++) {
      u[0][i] = std::sin(M_PI*dx*i);
    }
  }
  boundary(u[0]);
  printf("%e\n", u[0][25]);
}

void output(int ti, double t, vvd &u){
  sprintf(u_filename, "../data/u/%.3f.csv", t);
  FILE *fp = fopen(u_filename, "w");
  fprintf(time_fp, "%f\n", t);

  fprintf(fp, "x,u\n");
  if(OUTPUT_TERMINAL) printf("time:%f\n", t);
  
  for(int i = 0; i < NX; i++) {
    if(OUTPUT_TERMINAL){
      if(i < NX-1) printf("%f ", u[ti][i]);
      else printf("%f\n", u[ti][i]);
    }
    fprintf(fp, "%e,%e\n", i*dx, u[ti][i]);
  }
  fclose(fp);
}

void boundary(vd &nu){
  nu[0] = nu[NX-1] = 0.0;
}