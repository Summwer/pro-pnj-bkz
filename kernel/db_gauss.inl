#ifndef G6K_DB_INL
#define G6K_DB_INL

#ifndef G6K_SIEVER_H
#error Do not include siever.inl directly
#endif

// Insert e in db and in cdb, NOT thread-safe
inline void Siever::insert_in_db(Entry &&e)
{
    histo[histo_index(e.len)] ++;
    CompressedEntry ce;
    ce.len = e.len;
    if( ce.c.size()>0) ce.c = e.c;
    db.emplace_back(std::move(e));
    ce.i = db.size()-1;
    cdb.emplace_back(std::move(ce));
}

inline bool Siever::insert_in_db_and_uid(std::array<ZT,MAX_SIEVING_DIM> &x)
{
    auto uid = uid_hash_table.compute_uid(x);
    if (uid_hash_table.insert_uid(uid) == false)
    {
        return false;
    }
    db.emplace_back();
    Entry& e = db.back();
    e.x = std::move(x);
    e.uid = std::move(uid);
    recompute_data_for_entry<Recompute::recompute_all_and_consider_otf_lift & (~Recompute::recompute_uid)>(e);
    histo[histo_index(e.len)] ++;
    CompressedEntry ce;
    ce.len = e.len;
    if( ce.c.size()>0) ce.c = e.c;
    ce.i = db.size()-1;
    cdb.emplace_back(std::move(ce));
    return true;
}

inline bool Siever::verify_integrity(bool test_bijection)
{
    bool ret = true;
    size_t const db_size = db.size();
    if(cdb.size() != db_size) { std::cerr << "FATAL: cdb size and db size do not match.\n"; return false; }
    if(l > r) { std::cerr << "FATAL: l > r. \n"; return false; }
    if(l + n != r) { std::cerr << "FATAL: l + n != r.\n"; return false; }
    

    // Some integrity checks break when no yr-coordinates are saved
    /*
    std::vector<bool> indices_set; // Note that std::vector<bool> is special...
    if(test_bijection)
    {
        indices_set.resize(db_size,false);
    }
    for(size_t ind = 0; ind < db_size; ++ind)
    {
        if (cdb[ind].i >= db_size) { std::cerr << "FATAL: cdb index out of range.\n"; return false; }
        if (!NOYR and cdb[ind].c != db[cdb[ind].i].c) { std::cerr << "db and cdb disagree on c.\n"; ret = false; }
        if(test_bijection)
        {
            if (indices_set[cdb[ind].i] == true) { std::cerr << "2 cdb elements index the same db element\n."; ret = false; }
            indices_set[cdb[ind].i] = true;
        }
        if (cdb[ind].len < 0.) { std::cerr << "negative length in cdb.\n"; ret = false; } // cannot happen due to rounding.
        if (std::abs(cdb[ind].len - db[cdb[ind].i].len) > 0.001) { std::cerr << "cdb and db disagree about len.\n"; ret = false; }
    }
    for (size_t ind = 0; ind < db_size; ++ind)
    {
        // surjectivity is guaranteed by the pidgeonhole principle. So any failure here must be a failure of verify_integrity itself.
        assert(indices_set[ind]);
    }
    */

    Entry tmp;
    for (size_t ind = 0; ind < db_size; ++ind)
    {
        tmp = db[ind];
        if (tmp.len < 0.) { std::cerr << "negative length in db.\n"; ret = false; }
        recompute_data_for_entry<Recompute::recompute_all>(tmp);
        if (tmp.x != db[ind].x)  { std::cerr << "db entry is inconsistent @ x\n"; ret = false; }
        if(!NOYR) {
        for(unsigned int t = 0; t < n; ++t)
        {
            if (std::abs( tmp.yr[t] - db[ind].yr[t] )>0.001) { std::cerr << "db entry is inconstent @ yr\n"; ret = false; }
        }
        }
        if (POPCOUNT and tmp.c != db[ind].c)  { std::cerr << "db entry is inconsistent @ c \n"; ret = false; }
        if (tmp.uid == 0 ) { std::cerr << "0 uid in db\n"; ret = false; }
        if (tmp.uid != db[ind].uid) { std::cerr << "db entry is inconsistent @ uid" << ind << "\n"; ret = false; }
        if (std::abs(tmp.len - db[ind].len) > 0.01) { std::cerr << "db entry is inconsistent @ len" << ind << " " << db[ind].len << " " << tmp.len << " " << tmp.x[0] << "\n"; ret = false; }
        if (uid_hash_table.check_uid(tmp.uid) == false) { std::cerr << "uid of db entry is not in hash database\n"; ret =false; }
        // Note: if check_uid returns true, we increment the collision counter. We undo this as appropriate.
    }
    if (uid_hash_table.hash_table_size() != db_size) { std::cerr << "Number of uids in db and uid_db are inconsistent\n"; ret = false; }
    return ret;
}


// Sample a vector.
Entry Siever::sample(unsigned int large)
{
    ATOMIC_CPUCOUNT(213)
    Entry e;

    size_t fullS = cdb.size();

    // if (large || fullS < 200 + 10*n || !params.sample_by_sums)
    // {
    //     std::fill(e.x.begin(), e.x.end(), 0);
    //     size_t np = n / 2;
    //     np = (large > 3) ? 0 : np;
    //     for (size_t j = np; j < n; ++j){
    //         e.x[j] = (rng() % 5) - 2*(j < n-1); // Making sure the last cordinate is non-zero
    //         for (size_t l = 0; l < large; ++l)
    //         {
    //             e.x[j] += (rng() % 7) - 3*(j < n-1);
    //         }
    //     }

    //     recompute_data_for_entry_babai<Recompute::recompute_all_and_consider_otf_lift>(e,np);
    //     return e;
    // }

    //Gauss sampler
    // if (large || fullS < 200 + 10*n || !params.sample_by_sums){
    if(true){
        std::fill(e.x.begin(), e.x.end(), 0);
        size_t np = 0;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        float s = *(std::max_element(sqrt_rr.begin(),sqrt_rr.end()));
        // for(size_t i = 0; i<n; ++i){
        //     std::cerr<<sqrt_rr[i]<<" ";
        // }
        // std::cerr<<std::endl;
        // std::cerr<<sqrt_rr.end()-sqrt_rr.begin()<<std::endl;
        // auto maxposition = std::max_element(sqrt_rr.begin(),sqrt_rr.end());
        // std::cerr<<*maxposition<<std::endl;
        float sigma;
        for (size_t j = np; j<n; ++j){
            sigma = s/sqrt_rr[j];
            std::normal_distribution<double> distribution(0.0, n);
            //std::normal_distribution<double> distribution(0.0, sigma);
            e.x[j] = round(distribution(generator));

            // std::cerr<<sigma<<","<<e.x[j]<<" ";
        }

        // e.x[n-1] = 1; //Why the last one must be non-zero??
        // std::cerr<<std::endl;
        
        recompute_data_for_entry_babai<Recompute::recompute_all_and_consider_otf_lift>(e,np);


        return e;

    }
}





// x_full is a pointer to a (temporary) buffer, len is the non-normalized length.
inline void Siever::lift_and_compare(ZT * const x_full, FT len, LFT const * const helper)
{
    CPUCOUNT(212);
    if (UNLIKELY(len == 0.)) return;
    if (UNLIKELY(len >= params.lift_radius * gh)) return;
    if (UNLIKELY(len < lift_bounds[l])) lift_and_replace_best_lift(x_full, l);

    // i must be signed such that if ll == 0, we terminate properly.
    int i = static_cast<signed int>(l) - 1;
    const int llb = static_cast<signed int>(ll);

    if (helper)
    {
        //float prelen = len;
        for (int k = 0; k < int(OTF_LIFT_HELPER_DIM); ++k, --i)
        {
            if (UNLIKELY(i < llb)) return;

            FT yi = std::inner_product(x_full+l-k, x_full+l, full_muT[i].cbegin()+l-k,  static_cast<FT>(helper[k]));
            int const c = -std::floor(yi+0.5);
            x_full[i] = c;
            yi += c;
            len += yi * yi * full_rr[i];

            assert(i>=0);
            if (UNLIKELY(len < lift_bounds[i])) lift_and_replace_best_lift(x_full, static_cast<unsigned int>(i));
            if (UNLIKELY(len >= lift_max_bound)) return;
        }
        //std::cerr << len-prelen << std::endl;
    }

    for (; i >= llb; --i)
    {

        FT yi = std::inner_product(x_full+i+1, x_full+r, full_muT[i].cbegin()+i+1,  static_cast<FT>(0.));
        int const c = -std::floor(yi+0.5);
        x_full[i] = c;
        yi += c;
        len += yi * yi * full_rr[i];

        if (UNLIKELY(len < lift_bounds[i])) lift_and_replace_best_lift(x_full, static_cast<unsigned int>(i));
        if (UNLIKELY(len >= lift_max_bound)) return;
    }
}


inline void Siever::lift_and_compare(const Entry& e)
{
    CPUCOUNT(211);
    ZT x[r];
    std::fill (x,x+l,0);

    for (unsigned int i = 0; i < n; ++i)
    {
        x[i + l] = e.x[i];
    }

    if( e.otf_helper.size() > 0 ) 
        lift_and_compare(x, e.len * gh, e.otf_helper.data());
    else
        lift_and_compare(x, e.len * gh, nullptr);
}

//inline void Siever::clear_histo()
//{
//    std::fill(histo.begin(), histo.end(), 0);
//}

inline void Siever::recompute_histo()
{
    CPUCOUNT(210);
    if(histo_valid)
        return;
    histo_valid = true;
    // could be done more efficiently using the sorting of cdb, but this function does not
    // matter much anyway.
    std::fill(histo.begin(), histo.end(), 0);
    for (size_t i = 0; i < db_size(); ++i)
        ++histo[histo_index(cdb[i].len)];
    return;
}

inline void Siever::statistics_lift_hash( unsigned int N, float* liftlen, float* hashlen  ) {
    ZT x_full[r];
    for(size_t s = 0; s < N; s++ ) {
        Entry e = db[cdb[s].i];
        LFT* helper = e.otf_helper.data();
        std::fill( x_full, x_full+l, 0);
        for( size_t i =0; i < n; i++ )
            x_full[i+l] = e.x[i];

        // lift length
        double len = 0.;
        int i = int(l)-1;
        for (size_t k = 0; k < OTF_LIFT_HELPER_DIM; ++k, --i)
        {
            FT yi = std::inner_product(x_full+l-k, x_full+l, full_muT[i].cbegin()+l-k,  static_cast<FT>(helper[k]));
            int const c = -std::floor(yi+0.5);
            x_full[i] = c;
            yi += c;
            len += yi * yi * full_rr[i];

            assert(i>=0);
        }
        liftlen[s] = len;

        // hashlen
        hashlen[s] = dual_hashes.sqlen( e.dual_helper ); 
    }
}

#endif
