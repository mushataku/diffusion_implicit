#define _USE_MATH_DEFINES
#include <cstdio>
#include <vector>
#include <time.h>
#include <random>
#include <direct.h>
using namespace std;

/**********************config*************************/
//熱平衡状態になっているかを確認
const int DEBUG_EQUILIBRIUM = 0;

/*****************************************************/

/******************計算条件********************/
//スピンの数
const int n = 100;
const double N = double(n);
// 温度T:k_bT/J は [0.1,5]
const double T_ini = 0.1;
const double T_fin = 8.0;
const double dT1 = 0.02;
const double dT2 = 0.02;
const double dT3 = 0.02;
//エネルギーのサンプルの数
const int energy_num = 1e6;
// 配位の独立性を高めるために繰り返す更新回数
const int between = 20;
/*********************************************/

/**********************一様乱数の定義*************************/
//非決定的な乱数生成機、これを使ってシードをランダムに選ぶ
random_device rnddev;
//シード値の指定に使う乱数
mt19937 mt(rnddev());
//一様乱数
uniform_int_distribution<int> rndn(0,n);
uniform_int_distribution<int> rnd2(0,1);
uniform_real_distribution<double> rndd(0.0,1.0);
/************************************************************/

/*******************************関数の定義*******************************/
//配列qに初期のスピン配位をランダムに入れる
void init_configuration(vector<double> &q);
void set_Tvector(vector<double> &vT);
//配位をupdate
void update(vector<double> &q, double T);

// eをget
double energy(vector<double> &q);

// 温度Tでの熱平衡状態を取得
void get_equilibrium(double T, int energy_num, vector<double> &q);

//比熱の解析式を計算
double C_analytic(double T);
double C_analytic2(double T);
//エネルギーの解析式を計算
double e_analytic(double T);
double e_analytic2(double T);
double sech(double x);
/***********************************************************************/

int main(){
  double tmp, dummy, dE;;
  int posi;
  double e_ave, e2_ave, e;//<e>と<e^2>の計算用
  double E_num, E2_num, E_ana, E_ana2;
  double C_num, C_ana, C_ana2;//比熱 
  double C_err, C_err2, E_err, E_err2;

  vector<double> vT;//計算する温度を入れる配列
  // 全てのスピンの状態を格納 q = {1,-1,1,-1,-1,...}
  vector<double> q(n);
  
  //初期化
  set_Tvector(vT);
  init_configuration(q);

  /*******************記録用のフォルダ、ファイル作成*******************/
  FILE *C_fp, *E_fp, *E_dist_fp;
  char C_filename[50], E_filename[50], E_dist_filename[50], dir[50];
  _mkdir("data");
  sprintf(dir, "data/spin%d", n);
  _mkdir(dir);
  sprintf(C_filename, "data/spin%d/heat_capacity.csv", n);
  C_fp = fopen(C_filename,"w");
  sprintf(E_filename, "data/spin%d/energy.csv", n);
  E_fp = fopen(E_filename,"w");
  /*****************************************************************/

  fprintf(C_fp, "T,numerical,analytic,analytic2,error,error2\n");
  fprintf(E_fp, "T,numerical,analytic,analytic2,error,error2,e2_numeric\n");
  printf("spin number:%d\n", n);

  //各温度について比熱を計算
  for(auto T : vT){
    /*******************記録用のフォルダ、ファイル生成*******************/
    sprintf(dir, "data/spin%d/energy_dist", n);
    _mkdir(dir);
    sprintf(E_dist_filename, "data/spin%d/energy_dist/spin%d_T%.3f.csv", n, n, T);
    E_dist_fp = fopen(E_dist_filename, "w");
    fprintf(E_dist_fp, "No,Energy,T=,%e\n",T);
    /*****************************************************************/

    // 温度Tでの熱平衡状態を取得
    get_equilibrium(T, energy_num, q);
    e_ave = 0.0;
    e2_ave = 0.0;

    // カノニカル分布にしたがったエネルギーをenergy_num個getする
    for(int i = 0; i < energy_num; i++) {
      //相関を消すために数個置きに配位を選ぶ -> なんか意味ない気がするなあ -> めっちゃ意味あった
      for(int j = 0; j < between; j++) {
        //配位をupdate
        update(q, T);      
      }
      e = energy(q);
      fprintf(E_dist_fp, "%d,%e\n", i+1, e);
      e_ave += e;
      e2_ave += e*e;
    }

    fclose(E_dist_fp);

    /************エネルギー分布からエネルギー・比熱を計算****************/
    e_ave /= energy_num;
    e2_ave /= energy_num;
    E_num = e_ave/N;
    E2_num = e2_ave/N;
    E_ana = e_analytic(T);
    E_ana2 = e_analytic2(T);
    E_err = abs(E_num-E_ana)/abs(E_ana) * 100.0;
    E_err2 = abs(E_num-E_ana2)/abs(E_ana2) * 100.0;
    C_num = (e2_ave - e_ave*e_ave)/(T*T)/N;
    C_ana = C_analytic(T);
    C_ana2 = C_analytic2(T);
    C_err = abs(C_num-C_ana)/C_ana * 100.0;
    C_err2 = abs(C_num-C_ana2)/C_ana2 * 100.0;
    /*****************************************************************/

    printf("T:%e C_numeric:%e C_analytic:%e C_analytic2:%e error:%e error2:%e\n", T, C_num, C_ana, C_ana2, C_err, C_err2);
    fprintf(C_fp,"%e,%e,%e,%e,%e,%e\n", T, C_num, C_ana, C_ana2, C_err, C_err2);
    fprintf(E_fp,"%e,%e,%e,%e,%e,%e,%e\n", T, E_num, E_ana, E_ana2, E_err, E_err2, E2_num);
  }
  fclose(E_fp);
  fclose(C_fp);
  return 0;
}

////////function/////////////

double C_analytic(double T){
  double bJ = 1.0/T;
  return bJ*bJ*sech(bJ)*sech(bJ);
}

double C_analytic2(double T){
  double bJ = 1.0/T;
  double mother = (1.0 + pow(tanh(bJ),N));
  double C = 1.0;
  C -= N*pow(sech(bJ),2.0)*pow(tanh(bJ),2.0*N-2.0)/(mother*mother);
  C += (N-1.0)*pow(sech(bJ),2.0)*pow(tanh(bJ),N-2.0)/mother;
  C -= 2.0*pow(tanh(bJ),N)/mother;
  
  return bJ*bJ*sech(bJ)*sech(bJ)*C;
}

double e_analytic(double T){
  double bJ = 1.0/T;
  return -tanh(bJ);
}
double e_analytic2(double T){
  double bJ = 1.0/T;
  double mother = (1.0 + pow(tanh(bJ),N));
  double tmp;
  tmp = tanh(bJ) + pow(sech(bJ),2.0)*pow(tanh(bJ),N-1.0)/mother;
  return -tmp;
}

double energy(vector<double> &q){
  double energy = 0.0;
  for(int j = 0; j < n; j++) energy += q[j]*q[(j+1)%n];
  return -energy;
}


void init_configuration(vector<double> &q){
  for(int i = 0; i < n; i++) {
    //randが偶数なら上向き、奇数なら下向き
    //if(rand()%2 == 0){
    if(rnd2(mt) == 0){
      q[i] = 1.0;
    }
    else{
      q[i] = -1.0;
    }
  }
}

void set_Tvector(vector<double> &vT){
  //0.1 ～ 2.0 はdT1刻みで、それ以降はdT2刻み
  double T = T_ini, dT = dT1;
  while(T < T_fin + 1e-9){
    vT.push_back(T);
    if(T > 1.0 - 1e-9) dT = dT2;
    if(T > 2.0 - 1e-9) dT = dT3;
    T += dT;
  }
}

void update(vector<double> &q, double T){
  int posi;
  double tmp, dummy, dE, q2;
  //配位が変わるまで繰り返す
  int accept = 0;
  while (accept == 0){
    /*****************************************/
  // i < n なら q[i] の中から一つだけflipさせる
  // i = n なら今の状態にとどまる
  posi = rndn(mt);
  if (posi == n)
  {
    accept = 1;
    continue;
    }
    /*****************************************/
    q2 = -q[posi];
    dE = 2.0 * q[posi] * (q[(posi -1 + n) % n] + q[(posi + 1) % n]);
    tmp = fmin(1.0, exp(-dE / T));
    //dummy = (double)rand()/RAND_MAX;
    dummy = rndd(mt);

    //受託する時はqをflip
    if (dummy < tmp){
      q[posi] = q2;
      accept = 1;
    }
    //棄却するときはそのまま
  }
}
void get_equilibrium(double T, int energy_num, vector<double> &q){

  double tmp, dummy, dE;
  int posi;
  for(int i = 0; i < n*n; i++){
    update(q, T);
    if(DEBUG_EQUILIBRIUM == 1){
      printf("num:%d energy:%e\n", i+1, energy(q));
    }
  }
}

double sech(double x){
  return 1.0/cosh(x);
}
