#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

using vd = std::vector<double>;
using vvd = std::vector<vd>;
using vvvd = std::vector<vvd>;
using vvvvd = std::vector<vvvd>;

const double Lx = 1.0;
const int NX = 10 + 1;
const double dx = Lx/(NX-1);
const double kappa = 1.0;
const double OMEGA = 1.8;
const double dt = 0.01;
const double tmax = 1.0 + 1e-9;
const int NT = 100+1;
const double lx = kappa*dt/dx/dx;
const double EPS = 1e-5;

void init(vvd &u);
void boundary(vd &nu);
void output(int ti, vvd &u);

int main(){
  vvd u(NT,vd(NX));
  vd nu(NX), a(NX);
  double t = 0.0;
  int ti = 0;
  double du, tmp;

  init(u);
  boundary(u[0]);
  boundary(nu);
  printf("CALICULATION START\n");
  output(ti, u);

  printf("lx:%f\n", lx);

  for(ti = 0; ti < NT-1; ti++) {
    printf("ti:%d\n", ti+1);
    // u[ti] から u[ti+1] を求める
    // まず仮に u[ti+1] を u[ti] とする
    u[ti+1] = u[ti];
    for(int i = 1; i < NX-1; i++) {
      a[i] = u[ti][i] + lx/2.0*(u[ti][i+1]+u[ti][i-1]-2.0*u[ti][i]);
    }
    du = 1.0;
    // u[ti+1] を求める
    while(du > EPS){
      // 仮の u[ti+1] から nu = (仮のu[ti+1]) を求める
      for(int i = 1; i < NX-1; i++) {
        nu[i] = OMEGA*( lx/2.0/(1.0+lx)*(u[ti+1][i+1]+nu[i-1]) + a[i]/(1.0+lx) );
        nu[i] += (1.0-OMEGA)*u[ti+1][i];
      }
      boundary(nu);
      // du を求めて判定
      du = 0.0;
      for(int i = 0; i < NX; i++) {
        du += std::abs(nu[i]-u[ti+1][i]);
      }
      du /= NX;
      u[ti+1] = nu;
      // du < EPS ならこの u[ti+1] をアクセプト
    }
    // u[ti+1] が求め終わったのでこれを出力
    output(ti+1, u);
  }
  return 0;
}

void init(vvd &u){
  for(int i = 0; i < NX; i++) {
    if(0.3*NX < i && i < 0.7*NX){
      u[0][i] = 1.0;
    }
    else{
      u[0][i] = 0.0;
    }
  }
}


void output(int ti, vvd &u){
  printf("time:%f\n", dt*ti);
  for(int i = 0; i < NX; i++) {
    if(i < NX-1) printf("%f ", u[ti][i]);
    else printf("%f\n", u[ti][i]);
  }
}

void boundary(vd &nu){
  nu[0] = nu[NX] = 0.0;
}