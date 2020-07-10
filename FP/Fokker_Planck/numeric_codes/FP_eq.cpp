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
// terminal に出力して時間確認 
const int OUTPUT_TERMINAL = 1;
const int DEBUG_f = 0;
// 初期条件 -> 0:0, 1:一様分布;
const int INITIAL = 0;

// mesh -> 0:linear, 1:log
const int R_LOGMESH = 0;
const int X_LOGMESH = 0;

const int CHECK_DIFFERENCE = 0;
const int ZURU = 0;
const int MAKE_INPUT_FILE = 0;
/*****************************************************/

/******************計算条件********************/
const double Rmax = 1e5;
const double Rmin = 0.0;
const int NR = 100 + 1;
const double Xmax = 50.0;
const double Xmin = 1e-5;
const int HALF_NX = 5000;
const int NX = 2 * HALF_NX + 1;
const double a = 1.5e-2;
const double Xp = std::pow(a*Rmax/std::sqrt(M_PI),1.0/3.0);
// const double phi_eps = 1e-3;
const double phi0 = 1e6;

const int ITERATION_MAX = 1000000;
const double OMEGA = 1.0;
const double dt = Rmax/100.0;
// const double dt = 1e1;
const double DT = dt;
const double TMAX = 1e7;
const int STEADY = 10;
const double SOR_EPS = 1e-4;
      // double SOR_EPS = 1e-6;

vd r(NR+1),rp(NR+1),rm(NR+1);
vd x(NX),xp(NX),xm(NX);
vd co_xp(NX), co_xm(NX);
vvd co_rp(NR,vd(NX)), co_rm(NR,vd(NX)), co_over(NR,vd(NX));
void grid_set();

double phi(double x);
vvd Source(NR,vd(NX));
void Source_set();
void Source_reset();
/*********************************************/
/*
以下出力枚数の調整用
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に配列 (DT, 2DT, ..., ENDTIME) をsetする
*/
const double T_EPS = 1.0e-3;
const double ENDTIME = TMAX;
vd TIME;
//出力時刻を set する関数
void TIME_set();
/*********************************************/


/*************************** diffusion を解く関数 ***************************/
void init(vvd &J);
void boundary(vvd &nJ);
int diffusion(vvd &J);
double SOR(vvd &J, vvd &A);
/*************************** diffusion を解く関数 ***************************/

/*********************************ファイル関連*********************************/
FILE *condition_fp = fopen("../data/condition.csv", "w");
FILE *time_fp = fopen("../data/output_time.csv", "w");
// その時刻における J の値をファイルにアウトプット

void output_rJ(double t, vvd &J);
void output_xJ(double t, vvd &J);
void output(double t, vvd &J);
void make_input_file(double t, vvd &J);
void conservation(double t, int OUT_X, vvd &J);
/*****************************************************************************/

int main(){
  // 時間計測
  clock_t start_t, end_t;
  start_t = time(NULL);

  vvd J(NR+1,vd(NX+1));
  int ti = 0, ite, cnt_steady = 0;
  double t = 0.0;
  TIME_set();
  grid_set();
  Source_set();

  init(J);
  fprintf(time_fp, "time\n");
  fprintf(condition_fp, "Rmax,dt,NR\n%e,%e,%d\n",Rmax,dt,NR);
  output(t, J);
  printf("SOR_EPS:%e\n", SOR_EPS);
  // printf("phi_eps:%e\n", phi_eps);
  printf("phi0:%e\n", phi0);
  printf("Rmax:%e\n", Rmax);
  printf("Xmax:%e\n", Xmax);
  printf("Xp:%e\n", Xp);
  printf("dt:%e\n", dt);
  printf("****************CALICULATION START****************\n");
  bool reset = 1;
  while(t < TMAX + T_EPS) {
    // J を更新する
    ite = diffusion(J);
    if(ite == -1){
      printf("FP equation does not converge when time = %e\n", t);
      return 0;
    }
    
    t += dt;
    // printf("J[%d][%d] = %e\n", NR-1,NX/2,J[NR-1][NX/2]);

    // printf("time:%f iteration:%d\n", t, ite);

    // if(t > 100.0 && reset){
      // printf("before reset, Source:%e\n", Source[0][NX/2]);
      // for(int ir = 0; ir < NR; ir++) {
      //   for(int ix = 0; ix < NX; ix++) {
      //     if(Source[ir][ix] > 0.0) printf("Source[%d][%d]=%e\n", ir,ix,Source[ir][ix]);
      //     Source[ir][ix] = 0.0;
      //   }
      // }
      // Source_reset();
      // printf("reset!! Source:%e\n", Source[0][NX/2]);
      // reset = 0;
    // }

    // nJ が求め終わったのでこれを出力
    if(ti < TIME.size() && t > TIME[ti] - T_EPS){
      if(OUTPUT_TERMINAL) printf("output when time:%f iteration:%d\n", t, ite);
      output(t, J);
      ti++;
      // 定常判定：反復回数が一回がある程度続いたら定常に達したと見なす
      if(ite == 1){
        cnt_steady++;
        if(cnt_steady == STEADY) break;
      }
    }
  }

  fclose(time_fp);
  fclose(condition_fp);
  printf("******************NORMAL END******************\n");
  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);
  return 0;
}

void init(vvd &J){
  if(INITIAL == 0){
    for(int ir = 0; ir < NR; ir++) {
      for(int ix = 0; ix < NX; ix++) {
        J[ir][ix] = 0.0;
      }
    }
  }
  if(INITIAL == 1){
    for(int ir = 0; ir < NR; ir++) {
      for(int ix = 1; ix < NX-1; ix++) {
        // if(std::abs(x[ix]) < 10.0) J[ir][ix] = 0.025;
        // else J[ir][ix] = 0.0;
        if(ir != 0) J[ir][ix] = 1.0/(4.0*M_PI*r[ir]*r[ir]);
      }
    }
  }
}

void TIME_set(){
  double t = DT;
  while(t < ENDTIME + T_EPS){
    TIME.push_back(t);
    t += DT;
  }
}

void grid_set(){
  if(R_LOGMESH){
    double logdr = (std::log10(Rmax)-std::log10(Rmin))/(NR-1);
    for(int ir = 0; ir < NR+1; ir++) {
      r[ir] = Rmin*std::pow(10.0,ir*logdr);
      
    }
    for(int ir = 0; ir < NR; ir++) {
      rp[ir] = std::sqrt(r[ir]*r[ir+1]);
    }
    for(int ir = 1; ir < NR; ir++) {
      rm[ir] = rp[ir-1];
    }
    rm[0] = 0.0;
  }
  else{
    double dr = (Rmax-Rmin)/(NR-1);
    for(int ir = 0; ir < NR+1; ir++) {
      r[ir] = dr*ir;
    }
    for(int ir = 0; ir < NR; ir++) {
      rp[ir] = 0.5*(r[ir]+r[ir+1]);
    }
    for(int ir = 1; ir < NR; ir++) {
      rm[ir] = rp[ir-1];
    }
    rm[0] = 0.0;

    // for(int ir = 0; ir < NR; ir++) {
    //   printf("r[%d]=%f rp[%d]=%f rm[%d]=%f\n", ir, r[ir], ir, rp[ir], ir, rm[ir]);
    // }

    // 初期条件設定時の 0 除算を避けるための処方箋
    // r[0] = 1e-20;
  }

  // x 方向のメッシュ生成
  if(X_LOGMESH){
    // NX--;
  //   double logdx = (std::log10(Xmax)-std::log10(Xmin))/(NR-1);
  //   for(int ir = 0; ir < NR+1; ir++) {
  //     r[ir] = Rmin*std::pow(10.0,ir*logdr);
      
  //   }
  //   for(int ir = 0; ir < NR; ir++) {
  //     rp[ir] = std::sqrt(r[ir]*r[ir+1]);
  //   }
  //   for(int ir = 1; ir < NR; ir++) {
  //     rm[ir] = rp[ir-1];
  //   }
  }
  else{
    double dx = Xmax/HALF_NX;
    for(int ix = 0; ix < NX; ix++) {
      x[ix] = -Xmax + dx*ix;
    }
    for(int ix = 0; ix < NX-1; ix++) {
      xp[ix] = 0.5*(x[ix]+x[ix+1]);
    }
    for(int ix = 1; ix < NX; ix++) {
      xm[ix] = xp[ix-1];
    }
  }

  // for(int ix = 1; ix < NX-1; ix++) {
  //   printf("x[%d]=%f xp[%d]=%f xm[%d]=%f\n", ix, x[ix], ix, xp[ix], ix, xm[ix]);
  // }

  double lx,lr;

  // // 各種更新時に使う各種係数の前計算
  for(int ir = 0; ir < NR; ir++) {
    for(int ix = 1; ix < NX-1; ix++) {
      lx = 0.5*dt/(xp[ix]-xm[ix]);
      lr = dt/phi(x[ix])/(rp[ir]*rp[ir]*rp[ir]-rm[ir]*rm[ir]*rm[ir]);
      co_xp[ix] = lx*phi(xp[ix])/(x[ix+1]-x[ix]);
      co_xm[ix] = lx*phi(xm[ix])/(x[ix]-x[ix-1]);
      // co_rp の決定
      if(ir != NR-1){
        co_rp[ir][ix] = lr*rp[ir]*rp[ir]/(r[ir+1]-r[ir]);
      }
      else{
        // 境界条件
        co_rp[ir][ix] = 0.5*3.0*dt*rp[ir]*rp[ir]/(rp[ir]*rp[ir]*rp[ir]-rm[ir]*rm[ir]*rm[ir]);
      }
      // co_rmの決定
      if(ir != 0){
        co_rm[ir][ix] = lr*rm[ir]*rm[ir]/(r[ir]-r[ir-1]);
      }
      else{
        co_rm[ir][ix] = 0.0;
      }
      co_over[ir][ix] = 1.0 + co_xp[ix] + co_xm[ix] + co_rp[ir][ix] + co_rm[ir][ix];
    }
  }

  // exit(-1);

}

void conservation(double t, int OUT_X, vvd &J){
  // 工事中
  FILE *conservation_fp;
  double sum = 0.0;
  for(int ir = 0; ir < NR; ir++) {
    double vol = 4.0*M_PI/3.0*(rp[ir]*rp[ir]*rp[ir]-rm[ir]*rm[ir]*rm[ir]);
    sum += J[ir][OUT_X]*vol;
  }
  if(OUTPUT_TERMINAL) printf("conservation:%e\n", sum);
  fprintf(conservation_fp, "%e,%e\n",t,sum);
}

void output_rJ(double t, vvd &J){
  double tmp;

  char filename[50];
  char dirname[50];
  for(int ix = NX/2; ix < NX; ix += (NX-1)/10) {
    sprintf(dirname, "../data/J/rJ/x_%.3f", x[ix]);
    sprintf(filename, "../data/J/rJ/x_%.3f/%.e.csv", x[ix], t);
    mkdir(dirname, S_IREAD|S_IWRITE);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "r,J_point,J\n");
    for(int ir = 0; ir < NR; ir++) {
      tmp = 4.0*M_PI*r[ir]*r[ir];
      fprintf(fp, "%e,%e,%e\n", r[ir], J[ir][ix], tmp*J[ir][ix]);
    }
    fclose(fp);
  }
}
void output_xJ(double t, vvd &J){
  double tmp;
  char filename[50];
  char dirname[50];
  for(int ir = 0; ir < NR; ir += (NR-1)/10*2) {
    tmp = 4.0*M_PI*r[ir]*r[ir];
    sprintf(dirname, "../data/J/xJ/r_%.e", r[ir]);
    sprintf(filename, "../data/J/xJ/r_%.e/%.e.csv", r[ir], t);
    mkdir(dirname, S_IREAD|S_IWRITE);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "x,J_point,J\n");
    for(int ix = 0; ix < NX; ix++) {
      fprintf(fp, "%e,%e,%e\n", x[ix], J[ir][ix], tmp*J[ir][ix]);
    }
    fclose(fp);
  }
}

void output(double t, vvd &J){

  output_rJ(t, J);
  output_xJ(t, J);
  fprintf(time_fp, "%f\n", t);
  
  // for(int ir = 0; ir < NR; ir++) {
  //   if(DEBUG_f){
  //     if(ir < NR-1) printf("ir:%d r:%e J:%f ", ir, r[ir], J[ir][ix]);
  //     else printf("ir:%d r:%e J:%f\n", ir, r[ir], J[ir][ix]);
  //   }
  //   // fprintf(fp, "%e,%e\n", r[ir], J[ir][ix]);
  // }
}

void make_input_file(double t, vvd &J){
  // 実装中
  char filename[50];
  sprintf(filename, "../input/Rmax%.0f.csv", t);
  FILE *fp = fopen(filename, "w");

  fclose(fp);
}


void boundary(vvd &J){
  // 使っていない（更新式に自然に取り込まれているため）
  for(int ix = 1; ix < NX-1; ix++) {
    // J[NR-1][ix] = 2.0*J[NR-2][ix]/(2.0-3.0*phi(x[ix])*(r[NR-1]-r[NR-2]));
  }
}

int diffusion(vvd &J){
  // J から nJ を求める
  vvd A(NR, vd(NX,0));
  double dJ;

  for(int ir = 0; ir < NR; ir++) {
    for(int ix = 1; ix < NX-1; ix++) {
      // J[ir][ix] += dt*Source[ir][ix];
      // A[ir][ix] = J[ir][ix];

      A[ir][ix] = J[ir][ix] + dt*Source[ir][ix];
    }
  }

  // SOR で陰的に J を更新
  for(int icnt = 1; icnt <= ITERATION_MAX; icnt++) {
    // J から 仮の nJ を求める
    dJ = SOR(J, A); // dJ = |nJ-J|

    if(icnt%(ITERATION_MAX/10) == 0) printf("dJ:%e when icnt:%d\n", dJ, icnt);

    // dJ < EPS ならこの J をアクセプト
    if(dJ < SOR_EPS){
      if(CHECK_DIFFERENCE) printf("dJ:%e iteration:%d\n", dJ, icnt);
      return icnt;
    }
  }

  // ITERATION_MAX 回反復しても収束しないなら -1 を返して終了
  return -1;
}

double SOR(vvd &J, vvd &A){
  double nJ, dJ = 0.0, tmp;
  int ir = 0, ix = 0;
  // 仮の f から nf (仮のnu)を求める
  for(ir = 0; ir < NR; ir++) {
    // xに関しては端点は更新しない
    for(ix = 1; ix < NX-1; ix++) {

      if(ir == 0){
        // J[ir-1][ix] の項はなし
        nJ = OMEGA*(co_xp[ix]*J[ir][ix+1] + co_xm[ix]*J[ir][ix-1] + co_rp[ir][ix]*J[ir+1][ix] + A[ir][ix])/co_over[ir][ix];
      }
      else if(ir == NR-1){
        // J[ir+1][ix] の項はなし
        nJ = OMEGA*(co_xp[ix]*J[ir][ix+1] + co_xm[ix]*J[ir][ix-1] + co_rm[ir][ix]*J[ir-1][ix] + A[ir][ix])/co_over[ir][ix];
      }
      else{
        nJ = OMEGA*(co_xp[ix]*J[ir][ix+1] + co_xm[ix]*J[ir][ix-1] + co_rp[ir][ix]*J[ir+1][ix] + co_rm[ir][ix]*J[ir-1][ix] + A[ir][ix])/co_over[ir][ix];
      }
      // tmp = nJ;// デバッグ用
      nJ += (1.0-OMEGA)*J[ir][ix];
      if(nJ < 0){
        printf("before minus nJ:%e at r:%e x:%e\n", tmp, r[ir], x[ix]);
        printf("after minus nJ:%e at r:%e x:%e\n", nJ, r[ir], x[ix]);
        printf("J[%d][%d] = %e at r:%e x:%e\n", ir, ix, J[ir][ix], r[ir], x[ix]);
        exit(-1);
      }
      // if(ZURU && ir == NR-1 && ix == (NX-1)/2) nJ = 0.0;
      dJ = std::max(std::abs(nJ-J[ir][ix])/std::min(nJ,J[ir][ix]), dJ);
      J[ir][ix] = nJ;
    }

  }
  return dJ;
  
}


double phi(double X){
  if(std::abs(X) < 1e-7) return phi0;
  else return a/(X*X*std::sqrt(M_PI));
  // return a/((phi_eps+x*x)*std::sqrt(M_PI));
}

void Source_set(){
  for(int ir = 0; ir < NR; ir++) {
    double Sj;
    if(ir == 0){
      double dV = 4.0*M_PI/3.0*(rp[0]*rp[0]*rp[0]-rm[0]*rm[0]*rm[0]);
      Sj = 1.0/dV;
    }
    else if(ir > 0) Sj = 0.0;
    for(int ix = 0; ix < NX; ix++) {
      Source[ir][ix] = phi(x[ix])*Sj/(4.0*M_PI);
    }
  }

  FILE *phi_fp = fopen("../data/phi.csv", "w");
  fprintf(phi_fp, "x,phi\n");
  for(int ix = 0; ix < NX; ix++) {
    fprintf(phi_fp, "%e,%e\n", x[ix], phi(x[ix]));
  }
  fclose(phi_fp);
}

void Source_reset(){
  for(int ir = 0; ir < NR; ir++) {
    for(int ix = 0; ix < NX; ix++) {
      Source[ir][ix] = 0.0;
    }
  }
}