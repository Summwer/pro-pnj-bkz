#include "load_lwechal.h"

double log_gh_svp(double d, double delta_bkz, int svp_dim, int n, int q);

vector<tuple<int,int,int>> decoupler(LWEchal* lwechal);

tuple<int,int,int> find_min_complexity(vector<tuple<int,int,int>> params);

tuple<int,int,int> gsa_params(LWEchal* lwechal);

void primal_lattice_basis(LWEchal* lwechal);

LWEchal* gen_lwechal_instance(int n, double alpha);



