#include <cmath>
#include <cstdio>
#include <vector>
#include <time.h>
#include <sys/stat.h>
#include <limits>
using vd = std::vector<double>;
using vvd = std::vector<vd>;
using vvvd = std::vector<vvd>;
using vvvvd = std::vector<vvvd>;

/**********************config*************************/
// terminal に出力して確認 
const int OUTPUT_TERMINAL = 0;
// 初期条件 -> 0:長方形, 1:sinπx*sinπy 
const int INITIAL = 1;
/*****************************************************/

/******************計算条件********************/
const double Lx = 1.0;
const double Ly = 1.0;
const int NY = 100 + 1;
const int NX = 100 + 1;
const double dx = Lx/(NX-1);
const double dy = Ly/(NY-1);
const double kappa = 0.1;
const double OMEGA = 1.0;
const double dt = 0.01;
const double TMAX = 2.0 + 1e-9;
const double lx = kappa*dt/dx/dx;
const double ly = kappa*dt/dy/dy;
const double EPS = 1e-9;

/*
以下出力枚数の調整用
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に配列 (DT, 2DT, ..., ENDTIME) をsetする
*/
const double T_EPS = 1.0e-10;
const double DT = 0.01;
const double ENDTIME = 2.0;
vd TIME;
//出力時刻を set する関数
void TIME_set();
/*********************************************/


/*************************** diffusion を解く関数 ***************************/
// void init(vvd &u);
void init(vvd &u);
void boundary(vvd &u);

int diffusion(vvd &u);
double SOR(vvd &u, vvd &a);
double analytic(int jx, int jy, double t);
/*************************** diffusion を解く関数 ***************************/

/*********************************ファイル関連*********************************/
// FILE *condition_fp = fopen("../data/condition.txt", "w");
// FILE *time_fp = fopen("../data/output_time.csv", "w");
FILE *u_fp;
char u_filename[50];
// その時刻における u の値をファイルにアウトプット
void output(double t, vvd &u);
/*****************************************************************************/

/****************************/
double error = 0.0;
/****************************/

int main(){
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);


  vvd u(NX, vd(NY));
  int ti = 0;
  double du, nu;
  double t = 0.0;
  TIME_set();
  init(u);
  // fprintf(time_fp, "time\n");
  output(t, u);
  printf("lx:%f\n", lx);
  // fprintf(condition_fp, "%d %d\n%e %e %e", NX, NY, kappa, Lx, Ly);
  printf("****************CALICULATION START****************\n");
  while(t < TMAX) {
    // u から nu を求める
    if(diffusion(u) == -1){
      printf("Diffusion eqeation does not converge!!!\n");
      printf("when time = %f\n", t);
      return 0;
    }
    t += dt;
    // nu が求め終わったのでこれを出力
    if(ti < TIME.size() && t > TIME[ti] - T_EPS){
      output(t, u);
      ti++;
    }
  }

  // fclose(time_fp);
  // fclose(condition_fp);
  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);
  return 0;
}



void init(vvd &u){
  if(INITIAL == 0){
    for(int jx = 0; jx < NX; jx++) {
      for(int jy = 0; jy < NY; jy++) {
        if(0.3*NX < jx && jx < 0.7*NX && 0.3*NY < jy && jy < 0.7*NY){
          u[jx][jy] = 1.0;
        }
        else{
          u[jx][jy] = 0.0;
        }
      }
    }
  }

  if(INITIAL == 1){
    for(int jx = 0; jx < NX; jx++) {
      for(int jy = 0; jy < NY; jy++) {
        u[jx][jy] = std::sin(M_PI*dx*jx)*std::sin(M_PI*dy*jy);
      }
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

void output(double t, vvd &u){
  // sprintf(u_filename, "../data/u/%.3f.csv", t);
  // FILE *fp = fopen(u_filename, "w");
  // fprintf(time_fp, "%f\n", t);

  // fprintf(fp, "x,y,u\n");
  if(OUTPUT_TERMINAL) printf("time:%f\n", t);
  
  error = 0.0;
  for(int jx = 1; jx < NX-1; jx++) {
    for(int jy = 1; jy < NY-1; jy++) {
      if(OUTPUT_TERMINAL){
        if(jy < NY-1) printf("%f ", u[jx][jy]);
        else printf("%f\n", u[jx][jy]);
      }
      // fprintf(fp, "%e,%e,%e\n", jx*dx, jy*dy, u[jy][jx]);
      double tmp = std::abs(u[jx][jy]/analytic(jx,jy,t)-1.0);
      error = std::max(tmp,error);
      if(std::isnan(tmp)){
        printf("u[%d][%d] = %e %e\n", jx, jy, u[jx][jy], analytic(jx,jy,t));
        exit(-1);
      }
    }

  }
  printf("t = %e error:%e\n", t, error);
  // fclose(fp);
}

void boundary(vvd &u){
  for(int jx = 0; jx < NX; jx++) {
    u[jx][0] = u[jx][NY-1] = 0;
  }
  for(int jy = 0; jy < NY; jy++) {
    u[0][jy] = u[NX-1][jy] = 0;
  }
}

int diffusion(vvd &u){
  // u から nu を求める
  vvd a(NX, vd(NY));
  int imax = 99999;
  double du;
  for(int jx = 1; jx < NX-1; jx++) {
    for(int jy = 1; jy < NY-1; jy++) {
      a[jx][jy] = u[jx][jy] + 0.5*lx*(u[jx+1][jy]+u[jx-1][jy]-2.0*u[jx][jy]);
      a[jx][jy] += 0.5*lx*(u[jx][jy+1]+u[jx][jy-1]-2.0*u[jx][jy]);
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

double SOR(vvd &u, vvd &a){
  double nu, du = 0.0, tmp;
  // 仮の u から nu (仮のnu)を求める
  for(int jx = 1; jx < NX-1; jx++) {
    for(int jy = 1; jy < NY-1; jy++) {
      tmp = u[jx+1][jy] + u[jx-1][jy] + u[jx][jy+1] + u[jx][jy-1];
      nu = OMEGA*(lx*0.5*tmp + a[jx][jy])/(1.0+2.0*lx);
      nu += (1.0-OMEGA)*u[jx][jy];
      du = std::max(std::abs(nu-u[jx][jy]), du);

      u[jx][jy] = nu;
    }
  }
  return du;
}


double analytic(int jx, int jy, double t){
  double x = dx*jx;
  double y = dy*jy;
  return exp(-2.0*kappa*M_PI*M_PI*t)*sin(M_PI*x)*sin(M_PI*y);
}
