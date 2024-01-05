
#include "load_lwechal.h"

//Download files
size_t dl_req_reply(void *buffer, size_t size, size_t nmemb, void *user_p)
{
	FILE *fp = (FILE *)user_p;
	size_t return_size = fwrite(buffer, size, nmemb, fp);
	//cout << (char *)buffer << endl;
	return return_size;
}

//http GET request for file download.
CURLcode dl_curl_get_req(const std::string &url, std::string filename)
{

	const char* file_name = filename.c_str();
	char* pc = new char[1024];//long enough
	strcpy(pc, file_name);

	FILE *fp = fopen(pc, "wb");

	//curl initilization 
	CURL *curl = curl_easy_init();
	// curl return
	CURLcode res = curl_global_init(CURL_GLOBAL_ALL);
	if (curl)
	{
		//set header of curl
		struct curl_slist* header_list = NULL;
		header_list = curl_slist_append(header_list, "User-Agent: Mozilla/5.0 (Windows NT 10.0; WOW64; Trident/7.0; rv:11.0) like Gecko");
		curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header_list);

		//0: not accept 1: accept
		curl_easy_setopt(curl, CURLOPT_HEADER, 0);

		//URL address
		curl_easy_setopt(curl, CURLOPT_URL, url.c_str());

		//set ssl
		curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, false);
		curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, false);

		//CURLOPT_VERBOSE: 1, debug detail
		curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);

		curl_easy_setopt(curl, CURLOPT_READFUNCTION, NULL);

		//set data accept function
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &dl_req_reply);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);

		curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1);

		//set requst the longest duration
		// curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 6); // set transport and time out time  
		curl_easy_setopt(curl, CURLOPT_TIMEOUT, 6);

		// open request
		res = curl_easy_perform(curl);
	}
	// free curl 
	curl_easy_cleanup(curl);
	// free source
	fclose(fp);

	return res;
}

bool isFileExists_ifstream(string name) {
    ifstream f(name.c_str());
    return f.good();
}


LWEchal* load_lwe_challenge(int n, double alpha_){
    /*
    Load LWE challenge from file or website.

    :param n: LWE dimension
    :param alpha: the *standard deviation* of the secret is alpha*q

    */
    int alpha = int(round(alpha_ * 1000));
    ostringstream os, fname;
    const char* path = "lwechallenge/";
    fname << path << n << "-" << setw(3) << setfill('0')<< alpha <<".txt";

    if(not isFileExists_ifstream(fname.str())){
        os << "https://www.latticechallenge.org/lwe_challenge/challenges/LWE_" << n << "_" << setw(3) << setfill('0')<< alpha <<".txt";
        string dl_get_url = os.str();

        if(NULL==opendir(path))
            mkdir(path,0775);
        dl_curl_get_req(dl_get_url, fname.str());

        cout<<"Download \""<<fname.str()<<"\" successfully!"<<endl;
    }
    else{
        cout<<"\""<<fname.str()<<"\" exists."<<endl;
    }

	
    ifstream  fin;
    fin.open(fname.str(),ios::in);
    
	
    LWEchal* lwechal = new LWEchal;

	fin>>lwechal->n;
	fin>>lwechal->m;
	fin>>lwechal->q;
	fin>>lwechal->alpha;

	
	
	vector<Z_NR<ZT>> c;
	char ch;
	fin >> ch;
	if(ch == '['){
		while (fin >> ch && ch != ']'){
			fin.putback(ch);
			c.resize(c.size() + 1);
			// print_vector(c,0,c.size());
			if (!(fin >> c.back())){
				c.pop_back();
				break;
			}
		}
	}
	// print_vector(c,0,c.size());
	lwechal->c = c;
	
	

	vector<vector<Z_NR<ZT>>> matrix;
	fin>>ch;
	if(ch == '['){
		int i = 0;
		matrix.resize(0);
		// cout<<ch<<endl;
		while(fin >> ch && ch == '['){
			matrix.resize(i+1);
			// cout<<i<<endl;
			while (fin >> ch && ch != ']')
			{
				fin.putback(ch);
				matrix[i].resize(matrix[i].size() + 1);
				// cout<<lwechal->A[i].size()<<endl;
				if (!(fin >> matrix[i].back()))
				{
					matrix[i].pop_back();
					break;
				}
				// cout<<lwechal->A[i][0]<<endl;
			}
			i++;
		}
	}
	// print_matrix(matrix);
	


	lwechal->A.resize(lwechal->m,lwechal->n);
	for(int i = 0; i < lwechal->m ; i++){
		for(int j = 0; j < lwechal->n; j++){
			lwechal->A[i][j] = matrix[i][j]; 
		}
	}

	// lwechal->A.print(os);
	// cout<<os.str()<<endl;
	fin.close();
    return lwechal;
}


//min the minimal samples we should select for solving LWE through uSVP solver.
double find_m_in_gsa(int d, FP_NR<FT> dvol){
    /*
    input: d: dim of lattice 
           dvol: log(vol(L))
    output: BKZ-beta for solving LWE instance.
    */
    int bbeta = -1;
    double margin = -1, lhs, rhs, prev_margin = -1, pprev_margin = -1,beta_low;
    for(int beta = 2; beta < d; beta++){
        lhs = sqrt((double)beta);
        rhs = bkzgsa_gso_len(dvol, d - beta, d, beta).get_d();
        // cout<<"lhs = "<<lhs<< ", rhs = "<<bkzgsa_gso_len(dvol, d - beta, d, beta)<<endl;
        if(lhs < rhs and bbeta == -1){
            margin = rhs / lhs;
            prev_margin = pprev_margin;
            bbeta = beta;
            // cout<<"rhs<<","<<lhs<<endl;
        }

        if(lhs > rhs)
            bbeta = -1;
        pprev_margin = rhs / lhs; 
    }
    // ddelta = compute_delta(bbeta).get_d() * pow(margin, 1. / (double) d);
    if(prev_margin != - 1)
        beta_low = log(margin) / (log(margin) - log(prev_margin));
    else
        beta_low = 0;
    // if(d == 1004){
    //     cout<<"dvol:"<<dvol<<endl;
    //     cout<<margin<<", "<< prev_margin <<endl;
    //     cout<<"bbeta = "<<bbeta<<", beta_low = "<< beta_low<<endl;
    //     throw "";
    // }
    return (double) bbeta - beta_low;
}


//ZZ_mat<ZT> &A, ZZ_mat<ZT> &B, ZZ_mat<FT> &S,vector<Z_NR<ZT>> &b, vector<Z_NR<ZT>> &u, vector<double> &mu,
// LWEchal* gen_LWE_instance_with_input_distribution( int n, int q, int m, map<int,rational<int>> D_e,map<int,rational<int>> D_s, bool verbosity){
LWEchal* gen_LWE_instance_with_input_distribution( int n, int q, int m, map<int,double> D_e,map<int,double> D_s, bool verbosity){

    /*
    constructor that builds a LWE instance
    :n: (integer) size of the secret s
    :q: (integer) modulus
    :m: (integer) size of the error e
    :D_e: distribution of the error e (dictionnary form)
    :D_s: distribution of the secret s (dictionnary form)
    */
    
    if(verbosity){
       printf_red("     LWE instance generation    ");
       cout<<"n="<<n<<"\t"<<"m="<<m<<"\t"<<"q="<<q<<endl;
    }

	LWEchal* lwechal = new LWEchal;

    //define the mean and sigma of the instance
    // pair<rational<int>,rational<int>>
    pair<double,double> mu_s_e = average_variance(D_e), 
	                    mu_s_s = average_variance(D_s);
                
    // //rational<int> mu_e(mu_s_e.first), s_e(mu_s_e.second), mu_s(mu_s_s.first), s_s(mu_s_s.second);
    // double mu_e = rational_cast<double>(mu_s_e.first),
	// 	   s_e = rational_cast<double>(mu_s_e.second), //variance 
	// 	   mu_s = rational_cast<double>(mu_s_s.first),
    //        s_s = rational_cast<double>(mu_s_s.second);
    double //mu_e = mu_s_e.first, //average
           s_e = mu_s_e.second, //variance 
		   //mu_s = mu_s_s.first, 
           s_s = mu_s_s.second;
    // cout<<"s_e:"<<s_e<<endl;

    lwechal->sigma = sqrt(s_e) ; //Standard Deviation
	lwechal->alpha = lwechal->sigma / (double) q;
    //        
    // cout << "double=" << sigma << endl;
    // cout<<mu_s_e.first<<endl;
    
    // printf("mu_e = %.2lf, s_e = %.2lf, mu_s = %.2lf, s_s = %.2lf\n",mu_e,s_e,mu_s,s_s);
    
    
    // int d = m+n+1;
    
    // mu.resize(d);
    // for(int i = 0; i < d; i++){
    //     if(i<m)
    //         mu[i] = mu_e;
    //     else if(i<d-1)
    //         mu[i] = mu_s;
    //     else
    //         mu[i] = 1.;
    // }
    
	// vector<double> S;
    // S.resize(lwechal->dim);
    // for(int i = 0; i < lwechal->dim; i++){
    //     if(i<m)
    //         S[i] = log(abs(s_e));
    //     else if(i<m+n)
    //         S[i] = log(abs(s_s));
    //     else
    //         S[i] = log(1);
    // }
	// double  Svol = accumulate(S.begin(),S.end(),0.);


    //Find the best m.
    double min_beta = -1;
    for(int test_m = m; test_m > 0; test_m--){
        double dvol = test_m * log(q) - (test_m*log(abs(s_e))+n*log(abs(s_s))) / 2.; //  - log(lwechal->sigma)*lwechal->dim;
        // cout<<"s_e = "<<s_e<<", "<<"s_s = "<<s_s<<endl;
        // cout<<"vol: "<<test_m * log(q)<<","<<(test_m*log(abs(s_e))+n*log(abs(s_s))) / 2.<<endl;
        int d = test_m+n+1;
        double est_beta = find_m_in_gsa(d, dvol);
        if(dvol <= 0 or est_beta<=0)
            break;
        // else
        //     cout<<"dim = "<<d<<", est_beta = "<<est_beta<<endl;
        if(min_beta == - 1 or min_beta > est_beta ){
            min_beta = est_beta;
            lwechal->m = test_m;
        }
    }
    m = lwechal->m;
    lwechal->n = n;
	lwechal->q = q;
	lwechal->dim = m+n+1;
    lwechal->dvol = m * log(q) - (m*log(abs(s_e))+n*log(abs(s_s))) / 2.;
    printf("dim = %d, m = %d, dvol = %3.11f, β = %3.4f\n\n", lwechal->dim, lwechal->m, lwechal->dvol.get_d(), min_beta);



    // // draw the secrets
    vector<Z_NR<ZT>> s,e;
    s.resize(n);
    e.resize(m);
    for(int i =0; i<n; i++){    
        s[i] = draw_from_distribution(D_s,n);
    }
    for(int i =0; i<m; i++)
        e[i] = draw_from_distribution(D_e,m);
    lwechal->A.gen_zero(m,n); 
    gen_samples(lwechal->A,m,n,q); 
    // lwechal->A.print_comma(cout);
    // compute the public value t and build a target
    // b= (s*A.T + e)(mod q)
    lwechal->A.transpose(); 
    vector_matrix_product(lwechal->c, s, lwechal->A); 
    Z_NR<ZT> q2;
    q2 = q;
    for(int i =0; i < m; i++){
        lwechal->c[i].add(lwechal->c[i],e[i]);
        lwechal->c[i].mod(lwechal->c[i],q2);
    }   
    // print_vector(b);
    lwechal->A.transpose(); 




	// draw matrix B and define the lattice: B = [[qI,0],[-A^T,I]], A is lattice samples. Each col is one group of sample.
	// ZZ_mat<ZT> B;
    build_LWE_lattice(lwechal->B,lwechal->A,q);
    // // B.print_comma(cout);
    vector<Z_NR<ZT>> tmp = lwechal->c;
    tmp.resize(m+n);
    kannan_embedding(lwechal->B, tmp);
   

    return lwechal;
}
