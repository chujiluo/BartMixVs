/*
 *  This file is copied from:
 *  
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  and modified by:
 *  Chuji Luo and Michael J. Daniels
 */
#include "heterbd.h"

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, 
             std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen,
             std::vector<std::vector<double>>& mr_vec)
{
  tree::npv goodbots;  //nodes we could birth at (split on)
  double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x
  
  if(gen.uniform() < PBx) { //do birth or death
    
    //--------------------------------------------------
    //draw proposal
    tree::tree_p nx; //bottom node
    size_t v,c; //variable and cutpoint
    double pr; //part of metropolis ratio from proposal and prior
    bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);
    
    //--------------------------------------------------
    //compute sufficient statistics
    size_t nr,nl; //counts in proposed bots
    double bl,br; //sums of weights
    double Ml, Mr; //weighted sum of y in proposed bots
    hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma);
    
    //--------------------------------------------------
    //compute alpha
    double alpha=0.0, lalpha=0.0, lalpha0=0.0;
    double lhl, lhr, lht;
    if((nl>=5) && (nr>=5)) { //cludge?
      lhl = heterlh(bl,Ml,pi.tau);
      lhr = heterlh(br,Mr,pi.tau);
      lht = heterlh(bl+br,Ml+Mr,pi.tau);
      
      alpha=1.0;
      lalpha = log(pr) + (lhl+lhr-lht);
      //lalpha0 = log(pr) + (lhl+lhr-lht);
      lalpha = std::min(0.0,lalpha);
    }
    
    //--------------------------------------------------
    //try metrop
    double mul,mur; //means for new bottom nodes, left and right
    double temp_mr;
    //double temp_mr0;
    
    double uu = gen.uniform();
    bool dostep = (alpha > 0) && (log(uu) < lalpha);
    if(dostep) {
      mul = heterdrawnodemu(bl,Ml,pi.tau,gen);
      mur = heterdrawnodemu(br,Mr,pi.tau,gen);
      
      temp_mr = exp(lalpha);
      //temp_mr0 = exp(lalpha0);
      
      x.birthp(nx,v,c,mul,mur,temp_mr);
      
      nv[v]++;
      mr_vec[v].push_back(temp_mr);
      //mr0_vec[v].push_back(temp_mr0);
      
      return true;
    } else {
      return false;
    }
  } else {
    //--------------------------------------------------
    //draw proposal
    double pr;  //part of metropolis ratio from proposal and prior
    tree::tree_p nx; //nog node to death at
    dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);
    
    //--------------------------------------------------
    //compute sufficient statistics
    double br,bl; //sums of weights
    double Ml, Mr; //weighted sums of y
    hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigma);
    
    //--------------------------------------------------
    //compute alpha
    double lhl, lhr, lht;
    lhl = heterlh(bl,Ml,pi.tau);
    lhr = heterlh(br,Mr,pi.tau);
    lht = heterlh(bl+br,Ml+Mr,pi.tau);
    
    double lalpha = log(pr) + (lht - lhl - lhr);
    lalpha = std::min(0.0,lalpha);
    
    //--------------------------------------------------
    //try metrop
    //double a,b,s2,yb;
    double mu;
    double temp_mr;
    //double temp_mr0;
    
    if(log(gen.uniform()) < lalpha) {
      mu = heterdrawnodemu(bl+br,Ml+Mr,pi.tau,gen);
      
      size_t nx_getv = nx->getv();
      nv[nx->getv()]--;
      
      x.deathp(nx,mu,temp_mr);
      
      auto mr_it = find(mr_vec[nx_getv].begin(), mr_vec[nx_getv].end(), temp_mr);
      int mr_index = distance(mr_vec[nx_getv].begin(), mr_it);
      mr_vec[nx_getv].erase(mr_vec[nx_getv].begin() + mr_index);
      
      //auto mr0_it = find(mr0_vec[nx_getv].begin(), mr0_vec[nx_getv].end(), temp_mr0);
      //int mr0_index = distance(mr0_vec[nx_getv].begin(), mr0_it);
      //mr0_vec[nx_getv].erase(mr0_vec[nx_getv].begin() + mr0_index);
      
      return true;
    } else {
      return false;
    }
  }
}