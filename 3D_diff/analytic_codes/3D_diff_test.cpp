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
const double Lz = 1.0;
const int NY = 100 + 1;
const int NX = 100 + 1;
const int NZ = 100 + 1;
const double dx = Lx/(NX-1);
const double dy = Ly/(NY-1);
const double dz = Lz/(NZ-1);
const double kappa = 1.0;
const double OMEGA = 1.8;
const double dt = 0.01;
const double TMAX = 0.5 + 1e-9;
const double lx = kappa*dt/dx/dx;
const double ly = kappa*dt/dy/dy;
const double lz = kappa*dt/dz/dz;
const double EPS = 1e-5;

/*
以下出力枚数の調整用
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に配列 (DT, 2DT, ..., ENDTIME) をsetする
*/
const double T_EPS = 1.0e-10;
const double DT = 0.01;
const double ENDTIME = 0.5;
vd TIME;
//出力時刻を set する関数
void TIME_set();
/*********************************************/


/*************************** diffusion を解く関数 ***************************/
void init(vvvd &u);
void boundary(vvvd &u);

int diffusion(vvvd &u);
double SOR(vvvd &u, vvvd &a);
double analytic(int jx, int jy, int jz, double t);
/*************************** diffusion を解く関数 ***************************/

/*********************************ファイル関連*********************************/
// FILE *condition_fp = fopen("../data/condition.txt", "w");
// FILE *time_fp = fopen("../data/output_time.csv", "w");
FILE *u_fp;
FILE *error_fp = fopen("../data/error.csv","w");
char u_filename[50];
// その時刻における u の値をファイルにアウトプット
void output(double t, vvvd &u);
/*****************************************************************************/

int main(){
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);


  vvvd u(NX, vvd(NY, vd(NZ)));
  int ti = 0;
  double du, nu;
  double t = 0.0;
  TIME_set();
  init(u);
  // fprintf(time_fp, "time\n");
  fprintf(error_fp, "time,error\n");
  output(t, u);
  printf("lx:%f\n", lx);
  // fprintf(condition_fp, "%d %d\n%e %e %e", NX, NY, kappa, Lx, Ly);
  printf("****************CALICULATION START****************\n");
  while(t < TMAX) {
    // u から nu を求める

    // if(diffusion(u) == -1){
    //   printf("Diffusion eqeation does not converge!!!\n");
    //   printf("when time = %f\n", t);
    //   return 0;
    // }
    printf("number of iteration :%d \n", diffusion(u));

    t += dt;
    // nu が求め終わったのでこれを出力
    if(ti < TIME.size() && t > TIME[ti] - T_EPS){
      output(t, u);
      ti++;
    }
  }

  // fclose(time_fp);
  // fclose(condition_fp);
  fclose(error_fp);
  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);
  return 0;
}



void init(vvvd &u){
  if(INITIAL == 0){
    for(int jx = 0; jx < NX; jx++) {
      for(int jy = 0; jy < NY; jy++) {
        for(int jz = 0; jz < NZ; jz++){
          if(0.3*NX < jx && jx < 0.7*NX && 0.3*NY < jy && 
                  jy < 0.7*NY && 0.3*NZ < jz && jz < 0.7*NZ){
            u[jx][jy][jx] = 1.0;
          }
          else{
            u[jx][jy][jz] = 0.0;
          }
        }
      }
    }
  }

  if(INITIAL == 1){
    for(int jx = 0; jx < NX; jx++) {
      for(int jy = 0; jy < NY; jy++) {
        for(int jz = 0; jz < NZ; jz++){
          u[jx][jy][jz] = std::sin(M_PI*dx*jx)*std::sin(M_PI*dy*jy)*std::sin(M_PI*dz*jz);
        }
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

void output(double t, vvvd &u){
  // sprintf(u_filename, "../data/u/%.3f.csv", t);
  // FILE *fp = fopen(u_filename, "w");
  // fprintf(time_fp, "%f\n", t);
  // fprintf(fp, "x,y,u\n");
  if(OUTPUT_TERMINAL) printf("time:%f\n", t);
  
  double error = 0.0;
  for(int jx = 1; jx < NX-1; jx++) {
    for(int jy = 1; jy < NY-1; jy++) {
      for(int jz = 1; jz < NZ-1; jz++) {
        if(OUTPUT_TERMINAL){
          if(jy < NY-1) printf("%f ", u[jx][jy][jz]);
          else printf("%f\n", u[jx][jy][jz]);
        }
        // fprintf(fp, "%e,%e,%e\n", jx*dx, jy*dy, u[jy][jx]);
        // error += std::abs(u[jx][jy][jz] - analytic(jx,jy,jz,t));
        error += std::abs(1.0 - u[jx][jy][jz]/analytic(jx,jy,jz,t));
      }
    }
  }
  error /= double(NX)*NY*NZ;
  printf("t = %e error:%e\n", t, error);
  fprintf(error_fp, "%e,%e\n", t, error);
  // fclose(fp);
}

void boundary(vvvd &u){
  for(int jx = 0; jx < NX; jx++) {
    u[jx][0][0] = u[jx][NY-1][0] = u[jx][0][NZ-1] = u[jx][NY-1][NZ-1] = 0;
  }
  for(int jy = 0; jy < NY; jy++) {
    u[0][jy][0] = u[NX-1][jy][0] = u[0][jy][NZ-1] = u[NX-1][jy][NZ-1] = 0;
  }
  for(int jz = 0; jz < NZ; jz++) {
    u[0][0][jz] = u[NX-1][0][jz] = u[0][NY-1][jz] = u[NX-1][NY-1][jz] = 0;
  }
}

int diffusion(vvvd &u){
  // u から nu を求める
  vvvd a(NX, vvd(NY, vd(NZ)));
  int imax = 99999;
  double du;
  for(int jx = 1; jx < NX-1; jx++) {
    for(int jy = 1; jy < NY-1; jy++) {
      for(int jz = 1; jz < NZ-1; jz++) {
        a[jx][jy][jz] = u[jx][jy][jz]
                      + 0.5*lx*(u[jx+1][jy][jz]+u[jx-1][jy][jz]-2.0*u[jx][jy][jz])
                      + 0.5*lx*(u[jx][jy+1][jz]+u[jx][jy-1][jz]-2.0*u[jx][jy][jz])
                      + 0.5*lx*(u[jx][jy][jz+1]+u[jx][jy][jz-1]-2.0*u[jx][jy][jz]);
      }
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

double SOR(vvvd &u, vvvd &a){
  double nu, du = 0.0, tmp;
  // 仮の u から nu (仮のnu)を求める
  for(int jx = 1; jx < NX-1; jx++) {
    for(int jy = 1; jy < NY-1; jy++) {
      for(int jz = 1; jz < NZ-1; jz++) {
        tmp = u[jx+1][jy][jz] + u[jx-1][jy][jz] + u[jx][jy+1][jz] + u[jx][jy-1][jz]
            + u[jx][jy][jz+1] + u[jx][jy][jz-1] ;
        nu = OMEGA*(lx*0.5*tmp + a[jx][jy][jz])/(1.0+3.0*lx);
        nu += (1.0-OMEGA)*u[jx][jy][jz];
        du = std::max(std::abs(nu-u[jx][jy][jz]), du);

        u[jx][jy][jz] = nu;
      }
    }
  }
  return du;
}


double analytic(int jx, int jy, int jz, double t){
  double x = dx*jx;
  double y = dy*jy;
  double z = dz*jz;
  return exp(-3.0*kappa*M_PI*M_PI*t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
}
