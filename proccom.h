
#define BIG_REAL 1e10;

inline double sqr(double x)
{
    return x*x;
}



// the basic exchange event data for persistence/exchange times
class event
{
public:
    event(real t_, const rvec refp_, real r_, real dr_) : 
        t(t_), r(r_), dr(dr_)
    {
        refp[XX] = refp_[XX];
        refp[YY] = refp_[YY];
        refp[ZZ] = refp_[ZZ];
    }

    real get_t() const { return t; } // get the time

    real get_zb() const { return refp[ZZ]; }
    real get_r() const { return r; } // get the nearest nbr dist at ev start
    real get_dr() const { return dr; } // get the jump distance

    const rvec &get_refp() const { return refp; }

protected:
    real t;
    real r;
    real dr;
    rvec refp;
}; 

// list of exchange events for persistence/exchange times
class event_list
{
public:
    event_list() 
    { 
        //min_z = BIG_REAL; max_z = -BIG_REAL; 
        min_r = BIG_REAL; max_r = -BIG_REAL; 
        max_dr = -BIG_REAL; 
    }

    // add a single event.
    void add(real dt, 
             const rvec refp,
             //real zb, real za, real zm, real ze, 
             real r, real dr)
    {
        //event_t.push_back(event(dt, zb, za, ze, zm, r, dr, refp));
        event_t.push_back(event(dt, refp, r, dr));

        //printf("Event: z=%g, r=%g, dr=%g\n", z, r, dr);
        if (r < min_r)
            min_r = r;
        if (r > max_r)
            max_r = r;

        if (dr > max_dr)
            max_dr = dr;
    }

    /* post-process the event list into 
       h = histogram of times
       z_wh = a weighted histogram of <t> as a function of start z
       z_inv_wh = a weighted histogram of <1/t> as a function of start z
       z_rho_wh =  a weighted histogram of rho as a function of z
       rh = a histogram of nearest nbr dists (at the start of the ev)
       wrh = a weighted hist of nearest nbr dists with times
       wrinvh = a weighted hist of nearest nbr dists with <1/times>
       N = the total number
       av = the average time
       av_err = the error on the average time
       av_inv = the average inverse time
       av_inv_err = the error on the average inverse time
       startMoment = the lowest moment to calculate (can be 0 or negative)
       Nmoments = the number of moments to calculate 
       moments = the moments (output)
       */
    void proc(hist &h, /*hist &zh, */
              weighted_hist &z_wh, 
              weighted_hist &z_inv_wh, 
              weighted_hist *z_rho_wh,
              hist *rh, weighted_hist *wrh, weighted_hist *wrinvh,
              hist &drh,
              int &N, double &av, double &av_err,
              double &av_inv, double &av_inv_err,
              int startMoment, int Nmoments, 
              double *moments, double *moments_err, 
              grid *grd);

    /* get min, max z value */
    //real get_min_z() const { return min_z; }
    //real get_max_z() const { return max_z; }
    /* get min, max r value */
    real get_min_r() const { return min_r; }
    real get_max_r() const { return max_r; }
    /* get min, max dr value */
    real get_min_dr() const { return 0; }
    real get_max_dr() const { return max_dr; }


    size_t get_N() const { return event_t.size(); }

protected:
    std::list<event> event_t;
    //real min_z, max_z;
    real min_r, max_r;
    real max_dr;
};

// data associated with a molecule 
class mol
{
public:
    //mol() { N=0; prev_com_set=false; com_calced=false; }

    mol(int molblock_index_, int mol_index_, int molnr_,
        std::list<eventpos> *spare_event_list_,
        std::list<anchor> *spare_anchor_list_) : 
        molnr(molnr_),
        molblock_index(molblock_index_),
        mol_index(mol_index_)
    { 
        N=0; 
        prev_com_set=false;  
        com_calced=false; 
        spare_event_list=spare_event_list_;
        spare_anchor_list=spare_anchor_list_;
        startpos_event = NULL;
    }

    // add an initial atom 
    void add_atom(int molblock_index_, 
                  int mol_index_, int atom_index, int group_nr) 
    { 
        N++; 
        if (N>1 && group!=group_nr)
            printf("ERROR: group was %d, now %d, N=%d, atom_index=%d\n", 
                   group, group_nr, N, atom_index);

        if ( (mol_index_!=mol_index) || (molblock_index_ != molblock_index))
            printf("ERROR: mol was (%d,%d), now (%d,%d), N=%d\n", 
                   molblock_index, mol_index, 
                   molblock_index_, mol_index_ ,N);
        group=group_nr; // we assume it's the same for all.
    }
 
    void add_vec(const rvec v, double mass, t_pbc *pbc) 
    {
        if (com_N)
        {
            rvec vdx;

            pbc_dx(pbc, v, com_orig, vdx); 

            com[XX] += (vdx[XX]+com_orig[XX])*mass;
            com[YY] += (vdx[YY]+com_orig[YY])*mass;
            com[ZZ] += (vdx[ZZ]+com_orig[ZZ])*mass;
            com_m += mass;
        }
        else
        {
            com_orig[XX] = v[XX];
            com_orig[YY] = v[YY];
            com_orig[ZZ] = v[ZZ];
            com_m = mass;
            com[XX] = v[XX]*mass;
            com[YY] = v[YY]*mass;
            com[ZZ] = v[ZZ]*mass;
        }
        com_N++;
    }

    // get center of mass ASSUMING calc_com has been called!!!! 
    void get_com(rvec ret) const
    {
        ret[XX]=com[XX];
        ret[YY]=com[YY];
        ret[ZZ]=com[ZZ];
    }

    void reset_com()
    {
        if (com_calced)
        {
            prev_com[XX] = com[XX];
            prev_com[YY] = com[YY];
            prev_com[ZZ] = com[ZZ];
            prev_com_set=true;
        }
        com[XX] = 0;
        com[YY] = 0;
        com[ZZ] = 0;
        com_m = 0.;
        com_N = 0;
        com_calced=false;
        for(std::list<eventpos>::iterator ev=events.begin(); ev != events.end();
            ++ev)
        {
            ev->reset_checked();
        }
        if (startpos_event)
        {
            startpos_event->reset_checked();
        }
    }

    // calculate the new center of mass & displacement
    void calc_com(t_pbc *pbc) 
    {
        if (!com_calced)
        {
            com[XX]/=com_m;
            com[YY]/=com_m;
            com[ZZ]/=com_m;
            com_calced=true;
            /*com_m=1.;*/
            if (prev_com_set)
            {
                pbc_dx(pbc, com, prev_com, dx);
            }
            else
            {
                dx[XX] = dx[YY] = dx[ZZ] = 0;
            }
        }
    }

    // get center of mass displacement ASSUMING calc_com has been called!!!! 
    void get_dx(rvec ret) const
    {
        ret[XX] = dx[XX];
        ret[YY] = dx[YY];
        ret[ZZ] = dx[ZZ];
    }

    // get molecule group id
    int get_group() const { return group; }
    // get the number of atoms
    int get_N() const { return N; }
    // get the molecule index
    int get_mol_index() const { return mol_index; }
    // get the molecule block
    int get_molblock_index() const { return molblock_index; }
    // get the molecule nr
    int get_molnr() const { return molnr; }

    std::list<eventpos>::iterator get_first() { return events.begin(); }
    // checks for the next event. returns with the first event found
    // run this until it returns NULL.
    eventpos* check_events(std::list<eventpos>::iterator &next,
                           std::deque<mol> &mols_ref, t_pbc *pbc, 
                           rvec com_tot, rvec com_tot_dx, 
                           real jumpdist, real weightdist, int closest_n);

    /* allocate and generate an event based on 

       mols_meas = the list of coordinates of all the molecules of interest
       mols_ref = the list of coordinates of the reference molecules (for 
                    strain corr)
       pbc = the pbc struct
       time = the current time
       ref_trans = the matrix to correct the refernce group with
       com_tot = the system com
       com_tot_dx = the system com change wrt to last step
       com_ref_tot = the reference com
       com_ref_tot_dx = the reference com change wrt to last step
       jumpdist = the event distance to check for
       weightdist = the strain reference weight distance
       closest_n = the strain reference closest number of ref particles
       type = the type of event
       calc_r = whether to calculate the distance to ref grp
       */
    void make_event(std::deque<mol> &mols_meas, std::deque<mol> &mols_ref, 
                    t_pbc *pbc, real time, matrix ref_trans,
                    rvec com_tot, rvec com_tot_dx, 
                    rvec com_ref_tot, rvec com_ref_tot_dx,
                    real jumpdist, real weightdist, int closest_n,
                    event_type type, event_direction evdir,
                    bool calc_r,
                    real ref_dist_sq);

    void set_event(eventpos *ep,
                   std::deque<mol> &mols_meas, std::deque<mol> &mols_ref, 
                   t_pbc *pbc, real time, matrix ref_trans,
                   rvec com_tot, rvec com_tot_dx, 
                   rvec com_ref_tot, rvec com_ref_tot_dx,
                   real jumpdist, real weightdist,int closest_n,
                   event_type type, event_direction evdir,
                   bool calc_r,
                   real ref_dist_sq)
    {
        ep->reinit(time, com, mols_meas, mols_ref, pbc, ref_trans,
                   com_tot, com_tot_dx,
                   com_ref_tot, com_ref_tot_dx,
                   weightdist, closest_n, type, evdir, calc_r, 
                   *spare_anchor_list, ref_dist_sq);
    }


    // remove a found event.
    void remove_event(eventpos* ep);

    void count_events(int &pers, int &exch, int &preexch);

    void write_files(FILE *f1, FILE *f2, t_pbc *pbc);

    eventpos *startpos_event; /* event associated with the start position 
                                 (used when draw_pdb=true) */
protected:
    int com_N; // number of atoms already added 
    real com_m; // the molecule's mass
    rvec com; // the molecule's center of mass 
    rvec com_orig; // the first atom, for pbc reference

    rvec prev_com;// the previous center of mass
    bool prev_com_set; // whether the previous com was set 
    rvec dx; // distance in com w.r.t. previous step

    bool com_calced;

    int molnr; // a sequential number 
    int molblock_index; /* each molecule has a block and a mol index */
    int mol_index; 

    int N; // molecule's size in number of atoms
    int group; // group id

    std::list<eventpos> events; // the molecule's events, (if any).

    std::list<eventpos> *spare_event_list;
    std::list<anchor> *spare_anchor_list;

};

class comsystem
{
public:
    comsystem(int natom, real jumpdist_, real weightdist_, int closest_n_,
              real max_group_dist, 
              real max_group_atom_dist, 
              bool correct_group_rotation_, bool fit,
              bool z_absolute_,
              real z_min, real z_max, real z_spacing_, real r_spacing_, 
              real dr_spacing_, 
              event_direction evdir_,
              bool calc_r_, bool do_z_hist, bool do_z_rho_hist_,
              bool draw_pdb_, real grid_spacing_); 
   
    // initalize an atom (making a molecule if neccesary  
    void init_atom(int atom_index, int molblock_index, int mol_index, 
                   int group_index, 
                   const char *resname, const char *atomname);

    //  a new coordinate for an atom
    void add_atom_coord(int atom_index, rvec coord, real mass, t_pbc *pbc)
    {
        //printf("id=%d, mass=%g, (%g, %g, %g)\n", atom_index, mass, coord[XX], coord[YY], coord[ZZ]);
        atom_mols[atom_index]->add_vec(coord, mass, pbc);
    }

    void set_start_coord(int atom_index, rvec coord)
    {
        if (draw_pdb)
        {
            start_coords_x[atom_index] = coord[XX];
            start_coords_y[atom_index] = coord[YY];
            start_coords_z[atom_index] = coord[ZZ];
        }
    }

    // reset all centers of masses for a new frame
    void reset_coms();

    // handle atomic coordinates after reading a frame.
    void handle_atoms(t_pbc *pbc,
                      rvec x[],
                      int natoms,
                      int ngroups,
                      int size, 
                      atom_id *index,
                      t_atom *atom,
                      bool make_pers_events,
                      bool make_exch_events, 
                      real time, 
                      bool writing_out_coords);



    // frame is read; handle events
    void handle_mols(t_pbc *pbc, 
                     bool make_pers_events,
                     bool make_exch_events, 
                     real time);

  
    // post-process events. 
    void process_events(FILE *logf, size_t Nbin, 
                        const char *pers_exch_histfn,
                        const char *zhistfn, 
                        const char *gridfn,
                        const char *rhist, 
                        const char *drhist, 
                        const output_env_t oenv); 

    // whether there are events left
    int events_left()
    {
        return Nevents_left; 
    }

    // get the closest reference group atom's distance to a specific molecule
    real get_closest_ref_dist_sq(mol &m, t_pbc *pbc);

    // get com dist
    real get_com_dist_sq(mol &m, t_pbc *pbc);

    // check whether the event should be counted 
    bool event_relevant(eventpos *ep)
    {
        if (check_z_min_max)
        {
            real zb=ep->get_refp()[ZZ];
            if (zb < z_min_meas || zb >= z_max_meas)
                return false;
        }
#if 0
        if (max_group_dist > 0)
        {
            return (ep->get_ref_dist_sq() <  max_group_dist);
        }
        if (max_group_atom_dist > 0)
        {
            return (ep->get_ref_dist_sq() < max_group_atom_dist);
        }
#endif
        return true;
    }

    // check whether to even make the event.
    bool check_make_event(mol &m, t_pbc *pbc, real &ref_dist_sq)
    {
        if (max_group_dist_sq > 0)
        {
            ref_dist_sq = get_com_dist_sq(m, pbc);
            if (ref_dist_sq > max_group_dist_sq)
            {
                return false;
            }
        }
        else if (max_group_atom_dist_sq > 0)
        {
            ref_dist_sq = get_closest_ref_dist_sq(m, pbc);
            if (ref_dist_sq > max_group_atom_dist_sq)
            {
                return false;
            }
        }
        return true;
    }


    // write a PDB file with distances
    void write_pdb(const char *filename, t_pbc *pbc);
    // write a raw pdb from coordinates
    //void write_pdb_x(const char *filename, int natoms, rvec x[]);
    /*void open_pdb(const char *filename, matrix box);
    real write_pdb_atom(int atom_index, int res, const char *resname,
                        const char *atomname, rvec coord, real mass,
                        t_pbc *pbc);
    void close_pdb();*/
    
protected:
    // check for events
    void check_events(t_pbc *pbc,
                      bool make_pers_events,
                      bool make_exch_events,
                      real time);


    void write_ratio(FILE *outf, 
                     const weighted_hist &pers_hist, 
                     const weighted_hist &exch_hist);

    void proc_pers_rho(const weighted_hist &zb_pers, //input persistence times
                       const hist &zrho,    //input density times
                       weighted_hist &zb_rho_pers); // output <t_p>/rho

    void proc_inv_inv(const weighted_hist &inv,
                      weighted_hist &inv_inv);


    // process the first frame to extract r min and max values and allocate 
    //histogram if neccesary.
    void proc_first_frame(t_pbc *pbc);
    // process for z rho histogram
    void proc_z_rho_hist(t_pbc *pbc);


    int Natoms;
    std::deque<mol> mols_meas; // the meas molecules (water)
    std::deque<mol> mols_ref; // the ref molecules  (ref)
    std::list<eventpos> spare_event_list; // a 'spare' list
    std::list<anchor> spare_anchor_list; // a 'spare' list

    std::vector<mol*> atom_mols; // the molecule index for each atom

    /* data for pdb output */
    std::vector<real> start_coords_x;
    std::vector<real> start_coords_y;
    std::vector<real> start_coords_z;
    std::vector<const char *> atom_names;
    std::vector<const char *> res_names;
    std::vector<int> res_nrs;

    FILE *xyzout;
    FILE *zout;

    //std::list<eventpos> events; // all events 
    event_list persistence_events;
    event_list exchange_events;

    real first_t;
    real last_t;
    real delta_t; // time differential w.r.t. previous frame
    bool first_t_set;
    int frame_index;

    real jumpdist; // the displacement distance to check for 
    real weightdist; // the strain field size
    int closest_n; // the number of ref. particles to take

    real max_group_dist_sq; /* the maximum distance from the COM of a ref group
                               a mol can be to count as an event (or any 
                               distance if <0). */
    real max_group_atom_dist_sq; /* the maximum distance from an atom in the 
                                    ref group a mol can be to count as an event 
                                    (or any distance if <0). */

    bool correct_group_rotation; /* whether the ref group rotation should be 
                                    corrected */
    bool fit; /* whether the ref group rotation should be corrected with lsq 
                 fits*/
    rvec *fit_coords; /* the original coordinates to fit to. */
    real *fit_m; /* masses to weight the fit with. */
    matrix ref_trans; /* the rotation to correct the ref. group with for this 
                         step*/
    bool ref_trans_set; /* whether ref_trans has already been set in the prev. 
                           step. */

    bool z_absolute; // whether to use absolute coords for z
    bool check_z_min_max; // whether to check z min and max
    real z_min_meas, z_max_meas; // min and max z to measure for.
    real z_spacing; // spacing of z histogram (if > 0)
    real r_spacing; // spacing of r histogram (if > 0)
    real dr_spacing; // spacing of r histogram (if > 0)
    bool calc_r; /* whehter to calculate the distance to the ref grp com */
    bool do_z_hist; // whether to do a histogram in z
    bool do_z_rho_hist; // whether to do a density histogram in z

    real grid_spacing; // grid spacing (if >0)

    event_direction evdir; // the direction in which to do calcs

    real Nanchor_av; 
    real anchor_d; // weighted average anchor distance
    int NNanchor_av;

    rvec com_meas_tot; // total center of mass of meas molss
    rvec com_meas_tot_dx; // change in total center of mass of meas mols

    dvec com_ref_tot; // total center of mass for ref mols
    dvec com_ref_tot_dx; // change in com of ref mols
    //real com_ref_tot; // com z of the refs
    //real com_ref_tot_dx; // delta com z of the refs. 
    rvec com_ref_ev; // com of ref mols or zero, if z_absolute (for ev creation)
    rvec com_ref_ev_dx; // com of ref mols or zero, if z_absolute (for ev creation)

    bool atom_coms_for_ref; /* whether to use a com for each atom in the 
                               ref group (true), or a com for each molecule 
                               (false). True for max_group_(atom)_dist>0 */

    rvec min_r_found; // min coordinate values found
    rvec max_r_found; // max coordinate values found
    //real min_z_found, max_z_found; // minimum and maximum found values in z
    int Nzbin;
    hist *z_rho_hist_meas; // histogram in z for meas mols
    hist *z_rho_hist_ref; // histogram in z for ref mols
    static const int z_slack=25; // number of extra bins in both z dirs

    int Nevents_left; // number of outstanding events

    bool draw_pdb; // whether to collect data for drawing a pdb with movements
    FILE *pdbout; 

    rvec prcomp; // the principal components
    int pdbcounter;

    grid *exch_grid;
    grid *pers_grid;
    grid *ref_grid;
};



