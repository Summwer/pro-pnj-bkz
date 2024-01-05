
#include "utils.h"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cstddef>
#include<string>
#include<curl/curl.h>

using namespace std;

//Download files
size_t dl_req_reply(void *buffer, size_t size, size_t nmemb, void *user_p);

//http GET request for file download.
CURLcode dl_curl_get_req(const std::string &url, std::string filename);

bool isFileExists_ifstream(string name);

LWEchal* load_lwe_challenge(int n, double alpha_);

// LWEchal* gen_LWE_instance_with_input_distribution( int n, int q, int m, map<int,rational<int>> D_e,map<int,rational<int>> D_s, bool verbosity);
LWEchal* gen_LWE_instance_with_input_distribution( int n, int q, int m, map<int,double> D_e,map<int,double> D_s, bool verbosity);