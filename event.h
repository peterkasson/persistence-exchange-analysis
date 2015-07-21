
class mol;

// 'anchor' (list of reference mols) for strain correction.
class anchor
{
public:
    void reinit(const mol* molp_, real weight_, const rvec origp_) 
    {
        molp=molp_;
        weight=weight_;
        origp[XX] = origp_[XX];
        origp[YY] = origp_[YY];
        origp[ZZ] = origp_[ZZ];
    }

    double get_weight() const { return weight; }
    const mol* get_molp() const { return molp; }

    const rvec &get_origp()  const { return origp; } 
protected:
    const mol* molp;
    real weight;
    rvec origp;
};

enum event_type
{
    no_event,
    persistence_event,
    exchange_event,
    initial_exchange_event // the first exchange time doesn't count.
};


// the directions (dimensions) in which to track the distance. 
enum event_direction 
{
    evdir_xyz,
    evdir_xy,
    evdir_z
};

// the data associated with determining an event
class eventpos
{
public:
    /* create a new possible event with a position and a
       list of lipid molecules. Makes references
       to the contents of the lipid_mols list.

       time = event start time
       r    =  event start pos
       mols_mes = the list molecules whose event times are measured
       mols_ref = the list of reference molecules 
       pbc = the periodic boundary conditions gmx object
       ref_trans = the reference transformation
    
       weight_sigma = the sigma of the weight gaussian for inclusion of 
                      reference molecules. Not used if < 1
       closest_n = the number of reference molecules to include (not used if <0)
       */
       
    void reinit(real time_, rvec r, 
                const std::deque<mol> &mols_meas, 
                const std::deque<mol> &mols_ref, 
                t_pbc *pbc, 
                matrix ref_trans,
                rvec com_meas_tot, rvec com_meas_tot_dx, 
                rvec com_ref_tot, rvec com_ref_tot_dx,
                real weight_sigma, int closest_n, event_type type_, 
                event_direction evdir_, bool calc_r_,
                std::list<anchor> &spare_anchor_list,
                real ref_dist_sq_);

    // check whether an event happened, based on a new
    // position and stored references to lipid_mols.
    bool check_event(rvec rdx, rvec com_meas_tot_dx, real event_dist);

    real weightfn(real dsq)
    {
        return exp( -(dsq/(2.*weight_sigmasq)) ) * weight_norm;
    }

    event_type get_type() const { return type; }
    real get_time() const { return time; }

    // return the number of anchors
    int get_N() const { return (int)(anchors.size()); }
    // return weighted average anchor dist
    real get_anchor_d() const { return anchor_d; }

    bool has_happened() const { return event_happened; }

    void write_all_pbc(FILE *out, t_pbc *pbc, rvec w);
    void write_all_nopbc(FILE *out, rvec w, t_pbc *pbc);

    // remember whether the event has already been checked
    bool been_checked() const { return checked; }
    void reset_checked() { checked=false; }

    void get_real_startpos(rvec rv) // get start r
    { 
        rv[XX]=real_start_pos[XX]; 
        rv[YY]=real_start_pos[YY]; 
        rv[ZZ]=real_start_pos[ZZ]; 
    }

    real get_r() const { return r_val; }
    real get_dr() const { return dr_val; }

    const rvec &get_refp() const { return refp; } // get the reference pos.

    real get_dist() const { return sqrt(dsq); }

    real get_ref_dist_sq() const { return ref_dist_sq; }
protected:
    std::list<anchor> anchors;
    real weight_sigmasq; 
    real weight_norm;

    //rvec startpos; // starting position of the event's particle.
    rvec dx; // the summed shift in position so far
    rvec com_dx; // total summed shift in com so far

    event_type type;

    real time;

    bool event_happened; // stores whether an event has already happened.
    event_direction evdir;

    bool use_weights; // whether to use weights at all.

    bool checked; // whether the event has already been checked

    rvec real_start_pos; // the starting position in real coords

    rvec refp; /* the reference position at the start relative to the ref pos.
                  and rotation-corrected if neccesary */

    real r_val; // the r value to return 
    real dr_val; // the dr value to return 

    real dsq; // the resulting distance squared

    real ref_dist_sq; // container for closest ref dist

    real anchor_d; // weighted average anchor distance at start
}; 

