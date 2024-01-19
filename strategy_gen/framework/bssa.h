#include "bkz_with_jump_simulator.h"


class BSSA{
    public:
        //implement tours rounds of BKZ-beta-J
        struct strategy{
            int beta; 
            int jump;
            int tours;
        };
        struct blocksize_strategy{
            vector<strategy> S; //pnj-BKZ strategy: 
            vector<double> l;
            pair<double,double> cum_GB_BKZ;
            pair<double,double> cum_avg_GB_BKZ;
            double cum_pr; //cumulated probability     
        };
        

        map<int,blocksize_strategy> BS; //key: beta;  value: blocksize_strategy
        double G_min, G_tmp;
        blocksize_strategy bs0;
        int beta_tmp;
        BKZJSim* sim;
        COST* cost;
        Params params;
        double pre_pr = 0.;

        BSSA(Params params,int d){
            this->sim = new BKZJSim(params,d);
            cost = new COST(params);
            this->params = params;
        }


        void print_strategy(vector<strategy> S);
        void print_BS(map<int,BSSA::blocksize_strategy> BS);
        void print_bs(blocksize_strategy bs);

        BSSA::blocksize_strategy min_tour_to_each_goal_beta(BSSA::blocksize_strategy bs, int beta, int jump, pair<double,double> goal_quality); //Find the minial tours for a pnj-BKZ-beta-J to reach a basis quality of goal_dsvp.

        bool pnjbkz_beta_loop( vector<double> &l, pair<double,double> &cum_GB, pair<double,double> &cum_avg_GB, double &cum_pr, int beta, int jump, tuple<int,int,double,double,double> &dsvp_t_, double &slope);
        pair<double,double> max_tour_for_pnjbkz_beta(BSSA::blocksize_strategy bs, int beta); //maximal tours for one pnj-bkz-beta to reduce the basis.
        pair<double,double> max_tour_for_pnjbkz_beta(vector<double> l, int beta);

        void bssa_est(vector<double> l0, int sbeta, int gbeta); //bssa estimation
        // void bssa_est_mul_node(vector<double> l0, int sbeta, int gbeta); //bssa estimation with more node



        // pair<double,double> strategy_verification( vector<double> l,vector<strategy> strategy);

};