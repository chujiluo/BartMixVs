#include "heterbart.h"

//--------------------------------------------------
void heterbart::pr()
{
  Rcpp::Rcout << "+++++heterbart object:\n";
  bart::pr();
}
//--------------------------------------------------
void heterbart::draw(double *sigma, rn& gen)
{
  //cout << "Running heterbart::draw" << std::endl;
  for(size_t j=0;j<m;j++) {
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) {
      allfit[k] = allfit[k]-ftemp[k];
      //di.y is a pointer which is set to r in setdata()
      r[k] = y[k]-allfit[k];
    }
    heterbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen,mr_vec);
    heterdrmu(t[j],xi,di,pi,sigma,gen);
    fit(t[j],xi,p,n,x,ftemp);
    for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
  }
  if(dartOn) {
    draw_s(nv,lpv,theta,gen);
    draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
    for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
  }
}