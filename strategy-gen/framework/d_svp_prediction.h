#include "cost.h"

tuple<int,int,double,double,double> dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, pair<double,double> cum_GB_BKZ = {0.,0.});
tuple<int,int,double,double,double> fixed_dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, bool worst_case = false);
tuple<int,int,double,double,double> progressive_dsvp_predict(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, pair<double,double> cum_GB_BKZ = {0.,0.});
// tuple<int,int,double,double,double> progressive_dsvp_predict2(vector<double> l, double cumulated_proba, COST* cost, int cost_model=1, pair<double,double> cum_GB_BKZ = {0.,0.});