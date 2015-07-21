
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <vector>
#include <list>
#include <deque>

#include <gromacs/typedefs.h>
#include <gromacs/pbc.h>
#include <gromacs/xvgr.h>
#include <gromacs/vec.h>
#include <gromacs/princ.h>
#include <gromacs/do_fit.h>

#include "hist.h"
#include "event.h"
#include "grid.h"
#include "proccom.h"



void event_list::proc(hist &h, /*hist &zh, */
                      weighted_hist &z_wh, 
                      weighted_hist &z_inv_wh, 
                      weighted_hist *z_rho_wh,
                      hist *rh, weighted_hist *wrh, weighted_hist *wrinvh,
                      hist &drh, 
                      int &N, double &av, double &av_err, 
                      double &av_inv, double &av_inv_err,
                      int startMoment, int Nmoments, double *moments, 
                      double *moments_err, grid *grd)
{
    av=0.;
    av_inv=0.;
    double block_av=0.;
    double block_av_sq=0.;
    double block_av_inv=0.;
    double block_av_inv_sq=0.;

    const unsigned long int Nblock=10;
    unsigned long int i=0;
    unsigned long int block_N=0;
    unsigned long int last_block_index=0;
    unsigned long int Nev=event_t.size();

    double moments_block[Nmoments];
    double moments_block_av_sq[Nmoments];

    for(int m=0;m<Nmoments;m++)
    {
        moments[m] = 0.;
        moments_block[m] = 0.; // current block avg
        moments_block_av_sq[m] = 0.; // average sq for all blocks
    }

    N=0;

    for(std::list<event>::iterator ev=event_t.begin(); ev!=event_t.end(); ++ev)
    {
        real t=ev->get_t();
        real zb=ev->get_zb();
        //real r=ev->get_r();
        real dr=ev->get_dr();

        rvec refp;
        refp[XX]=ev->get_refp()[XX];
        refp[YY]=ev->get_refp()[YY];
        refp[ZZ]=ev->get_refp()[ZZ];
        real r=sqrt(refp[XX]*refp[XX] + refp[YY]*refp[YY] + refp[ZZ]*refp[ZZ]);

        // this should work if size<MAXINT/nblock:
        unsigned long int block_index=(Nblock*i)/Nev;

        h.add(log(t));

        av += t;
        av_inv += 1./t;

        real mom=1.;
        for(int m=0;m<(Nmoments+startMoment);m++)
        {
            moments[m-startMoment] += mom;
            moments_block[m-startMoment] += mom;
            mom *= t;
        }
        mom=1.;
        for(int m=-1;m>=startMoment;m--)
        {
            mom /= t;
            moments[m-startMoment] += mom;
            moments_block[m-startMoment] += mom;
        }

        N++;

        block_av += t;
        block_av_inv += 1./t;
        block_N++;

        z_wh.add(zb, t);
        z_inv_wh.add(zb, 1./t);

        if (rh)
            rh->add(r);
        if (wrh)
            wrh->add(r, t);
        if (wrinvh)
            wrinvh->add(r, 1./t);

        drh.add(dr);

        if (grd)
        {
            grd->add(t, ev->get_refp());
        }

        if (last_block_index != block_index || ((i+1) == Nev))
        {
            //printf("Doing block %ld.\n", block_index);
            for(int m=0;m<Nmoments;m++)
            {
                moments_block[m] /= block_N;
                moments_block_av_sq[m] += moments_block[m]*moments_block[m];
                moments_block[m] = 0;
            }

            block_av /= block_N;
            block_av_sq += block_av*block_av;

            block_av_inv /= block_N;
            block_av_inv_sq += block_av_inv*block_av_inv;

            // now do z_rho_wh
            if (z_rho_wh)
            {
                for(size_t j=0; j<z_wh.get_N(); ++j)
                {
                    //if (z_wh.get_cur_block_count(j)>0)
                    z_rho_wh->add_bin(j, z_wh.get_cur_block_count(j));
                }
                z_rho_wh->finish_block();
            }

            // reset block averages
            block_av=0;
            block_av_inv=0;
            block_N=0;
            z_wh.finish_block();
            z_inv_wh.finish_block();
            if (wrh)
                wrh->finish_block();
            if (wrinvh)
                wrinvh->finish_block();
        }

        i++;
        last_block_index=block_index;
    }
    if (block_N > 0)
    {
        printf("ERROR block_N > 0!\n");
        exit(1);
    }

    av/=N;
    block_av_sq /= Nblock;
    av_err=sqrt( (block_av_sq - av*av)/(Nblock-1) );

    av_inv/=N;
    block_av_inv_sq /= Nblock;
    av_inv_err=sqrt( (block_av_inv_sq - av_inv*av_inv)/(Nblock-1) );

    for(int m=0;m<Nmoments;m++)
    {
        moments[m] /=N;
        moments_block_av_sq[m] /= Nblock;
        moments_err[m] = sqrt( (moments_block_av_sq[m]- moments[m]*moments[m])/
                               (Nblock - 1) );
    }
}




eventpos* mol::check_events(std::list<eventpos>::iterator &next, 
                            std::deque<mol> &mols_ref, t_pbc *pbc, 
                            rvec com_meas_tot, rvec com_meas_tot_dx, 
                            real jumpdist, real weightdist, int closest_n)
{
    calc_com(pbc);

    if (startpos_event)
    {
        startpos_event->check_event(dx, com_meas_tot_dx, jumpdist);
    }

    std::list<eventpos>::iterator ev=next;
    while(ev!=events.end()) 
    {
        ++next;
        if ( (!ev->been_checked()) &&
            ev->check_event(dx, com_meas_tot_dx, jumpdist))
        {
            // remove event
            spare_event_list->splice(spare_event_list->begin(), events, ev);
            return &(*ev);
        }
        ev=next; // we need to be sure to stay on the list
    }
    return NULL;
}

void mol::write_files(FILE *f1, FILE *f2, t_pbc *pbc)
{
    calc_com(pbc);
    for(std::list<eventpos>::iterator ev=events.begin(); ev != events.end();
        ++ev)
    {
        ev->write_all_pbc(f1, pbc, com);
        ev->write_all_nopbc(f2, com, pbc);
    }
}


void mol::make_event(std::deque<mol> &mols_meas, std::deque<mol> &mols_ref, 
                     t_pbc *pbc, real time, matrix ref_trans,
                     rvec com_tot, rvec com_tot_dx, 
                     rvec com_ref_tot, rvec com_ref_tot_dx,
                     real jumpdist, real weightdist, int closest_n,
                     event_type type, event_direction evdir, bool calc_r,
                     real ref_dist_sq)
{
    // allocate if needed
    if (spare_event_list->empty())
    {
        for(int i=0;i<1000;i++)
        {
            spare_event_list->push_back(eventpos());
        }
    }


    std::list<eventpos>::iterator newelem=spare_event_list->begin();
    set_event(&(*(newelem)), mols_meas, mols_ref, pbc, time,  ref_trans,
              com_tot, com_tot_dx, com_ref_tot, com_ref_tot_dx, 
              jumpdist, weightdist, closest_n, type, evdir, calc_r,
              ref_dist_sq);
    events.splice(events.end(), *spare_event_list, newelem);
}

void mol::remove_event(eventpos *ep)
{
    for(std::list<eventpos>::iterator ev=events.begin(); ev != events.end(); 
        ++ev)
    {
        if ( &(*ev) == ep) 
        {
            spare_event_list->splice(spare_event_list->begin(), events, ev);
            return;
        }
    }

    printf("ERROR: Could not remove event!!\n");
}

void mol::count_events(int &pers, int &exch, int &preexch)
{
    for(std::list<eventpos>::iterator ev=events.begin(); ev != events.end(); 
        ++ev)
    {
        switch (ev->get_type())
        {
            case persistence_event:
                pers++;
                break;
            case exchange_event:
                exch++;
                break;
            case initial_exchange_event:
                preexch++;
                break;
            default:
                printf("ERROR counting events!!");
                exit(1);
        }
    } 
}


comsystem::comsystem(int natoms, real jumpdist_, real weightdist_, 
                     int closest_n_,
                     real max_group_dist_, 
                     real max_group_atom_dist_, 
                     bool correct_group_rotation_, bool fit_,
                     bool z_absolute_, real z_min_, real z_max_,
                     real z_spacing_, real r_spacing_, 
                     real dr_spacing_, 
                     event_direction evdir_, 
                     bool calc_r_, 
                     bool do_z_hist_, bool do_z_rho_hist_,
                     bool draw_pdb_, real grid_spacing_) 
    : Natoms(natoms),
      atom_mols(natoms),
      start_coords_x(draw_pdb_ ? natoms : 0),
      start_coords_y(draw_pdb_ ? natoms : 0),
      start_coords_z(draw_pdb_ ? natoms : 0),
      atom_names(draw_pdb_ ? natoms : 0),
      res_names(draw_pdb_ ? natoms : 0),
      res_nrs(draw_pdb_ ? natoms : 0),
      jumpdist(jumpdist_), weightdist(weightdist_), closest_n(closest_n_),
      max_group_dist_sq(max_group_dist_ > 0 ? 
                        max_group_dist_*max_group_dist_ : -1.), 
      max_group_atom_dist_sq(max_group_atom_dist_ > 0 ?
                             max_group_atom_dist_*max_group_atom_dist_ : -1.), 
      correct_group_rotation(correct_group_rotation_), 
      fit(fit_), fit_coords(NULL), fit_m(NULL),
      z_absolute(z_absolute_), z_min_meas(z_min_), z_max_meas(z_max_),
      z_spacing(z_spacing_), r_spacing(r_spacing_), 
      dr_spacing(dr_spacing_), calc_r(calc_r_), 
      do_z_hist(do_z_hist_), do_z_rho_hist(do_z_rho_hist_),
      grid_spacing(grid_spacing_),
      evdir(evdir_), draw_pdb(draw_pdb_),
      exch_grid(NULL), pers_grid(NULL), ref_grid(NULL)
{
    pdbcounter=0;

    xyzout = NULL;
    zout = NULL;
    //xyzout = fopen("movie.xyz", "w");
    //zout = fopen("z.xvg", "w");

    printf("natoms=%d\n",natoms);

    first_t_set = false;
    frame_index=0;

    Nanchor_av=0;
    NNanchor_av=0;
    anchor_d=0;

    min_r_found[XX] = 0;
    min_r_found[YY] = 0;
    min_r_found[ZZ] = 0;
    max_r_found[XX] = 0;
    max_r_found[YY] = 0;
    max_r_found[ZZ] = 0;

    //min_z_found = 0;
    //max_z_found = 0;
    z_rho_hist_meas=NULL;
    z_rho_hist_ref=NULL;

    Nevents_left=0;

    if (z_min_meas > z_max_meas)
    {
        check_z_min_max = false;
    }
    else
    {
        check_z_min_max = ( (z_min_meas > -1000) && (z_max_meas < 1000) );
    }

    atom_coms_for_ref=false;
    if  (max_group_atom_dist_sq > 0.)
    {
        atom_coms_for_ref=true;
        fprintf(stderr, "Setting up for max. atom distance from group\n");
    }
    else if (max_group_dist_sq > 0.)
    {
        fprintf(stderr, "Setting up for max. distance from group\n");
    }

    // initialize to the unity matrix
    clear_mat(ref_trans);
    ref_trans[XX][XX]=1.;
    ref_trans[YY][YY]=1.;
    ref_trans[ZZ][ZZ]=1.;
    ref_trans_set=false;

    if ( (max_group_dist_sq > 0.) && (max_group_atom_dist_sq > 0.) )
    {
        fprintf(stderr, "both max_group_dist and max_group_atom_dist > 0\n");
        exit(1);
    }
}

void comsystem::init_atom(int atom_index, int molblock_index, int mol_index, 
                          int group_index, 
                          const char *resname, const char *atomname)
{
    bool found=false;
    bool newmol=false;
    // check if the molecule already exists

    std::deque<mol> *mols = (group_index==0) ? &mols_meas : &mols_ref ;
    std::deque<mol>::iterator ret=mols->begin();

    if (! atom_coms_for_ref)
    {
        if ( (ret->get_mol_index() == mol_index) &&
             (ret->get_molblock_index() == molblock_index) )
        {
            found=true;
        }
#if 0
        // IF COMMENTED, WE ASSUME THAT MOLECULES ARE ALWAYS SEQUENTIAL
        // IN ATOM NUMBER
        else
        {
            for(; ret!=mols->end(); ++ret)
            {
                if ( (ret->get_mol_index() == mol_index) && 
                     (ret->get_molblock_index() == molblock_index) )
                {
                    found=true;
                    fprintf(stderr, "ERROR ERROR mol out of order!\n");
                    break;
                }
            }
        }
#endif
    }

    if (!found)
    {
        mols->push_front(mol(molblock_index, mol_index, mols->size(),
                             &spare_event_list, &spare_anchor_list));
        ret = mols->begin();
        newmol=true;
    }

    ret->add_atom(molblock_index, mol_index, atom_index, group_index);
    atom_mols[atom_index]= &(*ret);
    if (draw_pdb)
    {
        atom_names[atom_index] = atomname;
        res_names[atom_index] = resname;
        res_nrs[atom_index] = mol_index;
    }
#if 0
    printf("Initizalizing atom %d %s mol (%d, %d, %d) (N=%d)\n", 
           atom_index,
           (newmol ? "new" : "   " ),
           atom_mols[atom_index]->get_molblock_index(), 
           atom_mols[atom_index]->get_mol_index(), 
           atom_mols[atom_index]->get_molnr(), 
           atom_mols[atom_index]->get_N());
#endif
}


void comsystem::reset_coms()
{
    for(std::deque<mol>::iterator i=mols_meas.begin();i!=mols_meas.end();++i)
    {
        i->reset_com();
    }
    for(std::deque<mol>::iterator i=mols_ref.begin();i!=mols_ref.end();++i)
    {
        i->reset_com();
    }
}

void comsystem::handle_atoms(t_pbc *pbc, 
                             rvec x[],
                             int natoms,
                             int ngroups,
                             int isize,
                             atom_id *index,
                             t_atom *atom,
                             bool make_pers_events,
                             bool make_exch_events, 
                             real time, 
                             bool writing_out_coords)
{
    if (fit)
    {
        if (!fit_m)
        {
            fit_m=new real[natoms];
            for(int i=0;i<natoms;i++)
            {
                // we set the weight to zero for all atoms, then set the
                // right weight for the atoms we're interested in
                fit_m[i]=0;
            }
            // set the initial weights.
            for(int i=0;i<isize;i++)
            {
                fit_m[index[i]] = atom[index[i]].m;
            }
        }

        // reset the com
        reset_x_ndim(3, isize, index, natoms, NULL, x, fit_m);

        if (fit_coords)
        {
            calc_fit_R(3, natoms, fit_m, fit_coords, x, ref_trans);
        }
        else 
        {
            // now set initial coordinates
            fit_coords=new rvec[natoms];
            for(int i=0;i<natoms;i++)
            {
                fit_coords[i][XX] = x[i][XX];
                fit_coords[i][YY] = x[i][YY];
                fit_coords[i][ZZ] = x[i][ZZ];
            }
        }
    }

    if (correct_group_rotation)
    {
        int     i,m;
        rvec    xcm; //,prcomp;
        matrix  ref_trans_new;

        // from gromacs' princ.c
        // calculate center of mass
        calc_xcm(x,isize,index,atom,xcm,FALSE);
        // and center the reference group at the COM
        for(i=0; i<natoms; i++)
            rvec_dec(x[i],xcm);
        principal_comp(isize,index,atom,x,ref_trans_new,prcomp);
        
        /* match up the new vectors with the old one. */
        if (ref_trans_set)
        {
            rvec old_vecs[DIM];
            rvec new_vecs[DIM];
            //rvec src;

            // decompose the matrix into 3 vectors
            for(int i=0;i<DIM;i++)
            {
                for(int j=0;j<DIM;j++)
                {
                    old_vecs[i][j] = ref_trans[i][j];
                    new_vecs[i][j] = ref_trans_new[i][j];
                }
            }
            // now find out the maximum dot product permutation
            int permutation[DIM];
            real multfact[DIM];
            for(int i=0;i<DIM;i++)
            {
                real biggest=-1e9;
                int biggest_index=-1;
                for(int j=0;j<DIM;j++)
                {
                    real dot=iprod(old_vecs[i], new_vecs[j]);
                    if (fabs(dot) > biggest)
                    {
                        biggest=fabs(dot);
                        permutation[i]=j;
                        if (dot > 0)
                            multfact[i]=1;
                        else 
                            multfact[i]=-1;
                    }
                }
                // check whether there's a duplicate permutation
                for(int j=0;j<i;j++)
                {
                    if (permutation[j] == permutation[i])
                    {
                        fprintf(stderr, "ERROR: permutation %d == %d\n", i, j);
                        fprintf(stdout, "ERROR: permutation %d == %d\n", i, j);
                        printf("old:\n");
                        printf("( %g, %g, %g )\n", ref_trans[XX][XX], 
                               ref_trans[XX][YY], ref_trans[XX][ZZ]);
                        printf("( %g, %g, %g )\n", ref_trans[YY][XX], 
                               ref_trans[YY][YY], ref_trans[YY][ZZ]);
                        printf("( %g, %g, %g )\n", ref_trans[ZZ][XX], 
                               ref_trans[ZZ][YY], ref_trans[ZZ][ZZ]);
                        printf("new:\n");
                        printf("( %g, %g, %g )\n", ref_trans_new[XX][XX], 
                               ref_trans_new[XX][YY], ref_trans_new[XX][ZZ]);
                        printf("( %g, %g, %g )\n", ref_trans_new[YY][XX], 
                               ref_trans_new[YY][YY], ref_trans_new[YY][ZZ]);
                        printf("( %g, %g, %g )\n", ref_trans_new[ZZ][XX], 
                               ref_trans_new[ZZ][YY], ref_trans_new[ZZ][ZZ]);

                       exit(1);
                    }
                }
            }        

            // now create the new matrix
            for(int i=0;i<DIM;i++)
            {
                for(int j=0;j<DIM;j++)
                {
                    ref_trans_new[i][j] = multfact[i]*
                                          new_vecs[permutation[i]][j];
                }
            }
        }

        if (fabs(det(ref_trans_new) - 1) > 1e-2 ) 
        {
            printf("Determinant != 1!! %g\n", det(ref_trans_new));
            exit(1);
        }


        copy_mat(ref_trans_new, ref_trans);
        ref_trans_set=true;
        /*we now have the transformation and can use it later on in the step*/
    }

    bool fill_ref_grid=(ref_grid && (frame_index % 100 == 0));
    if (writing_out_coords || fill_ref_grid )
    {
        /* do the actual rotation */
        rotate_atoms(natoms,NULL,x,ref_trans);
    }

    if (fill_ref_grid)
    {
        for(int i=0;i<isize;i++)
        {
            ref_grid->add(1, x[index[i]]);
        }
        ref_grid->add_sample();
        ref_grid->set_sample_size(isize);
    }
}


void comsystem::handle_mols(t_pbc *pbc, 
                            bool make_pers_events,
                            bool make_exch_events, 
                            real time)
{
    if (!first_t_set)
    {
        first_t_set=true;
        first_t=time;
        delta_t=0;
    }
    else
    {
        delta_t = time - last_t;
    }
    last_t=time;

    com_meas_tot_dx[XX] = com_meas_tot_dx[YY] = com_meas_tot_dx[ZZ] = 0.;
    /* first calculate the centers of mass */
    for(std::deque<mol>::iterator i=mols_meas.begin();i!=mols_meas.end();++i)
    {
        rvec dx;
        i->calc_com(pbc);
        i->get_dx(dx);

        com_meas_tot_dx[XX] += dx[XX];
        com_meas_tot_dx[YY] += dx[YY];
        com_meas_tot_dx[ZZ] += dx[ZZ];
    }
    com_ref_tot_dx[XX] = com_ref_tot_dx[YY] = com_ref_tot_dx[ZZ] = 0.;
    for(std::deque<mol>::iterator i=mols_ref.begin();i!=mols_ref.end();++i)
    {
        rvec dx;
        i->calc_com(pbc);
        i->get_dx(dx);

        com_ref_tot_dx[XX] += dx[XX];
        com_ref_tot_dx[YY] += dx[YY];
        com_ref_tot_dx[ZZ] += dx[ZZ];
    }

    // calculate new coms
    com_meas_tot_dx[XX]/=mols_meas.size();
    com_meas_tot_dx[YY]/=mols_meas.size();
    com_meas_tot_dx[ZZ]/=mols_meas.size();
    com_meas_tot[XX] += com_meas_tot_dx[XX];
    com_meas_tot[YY] += com_meas_tot_dx[YY];
    com_meas_tot[ZZ] += com_meas_tot_dx[ZZ];

    //printf("mols_ref.size()=%ld\n", mols_ref.size());
    size_t mrsize=mols_ref.size();
    if (mrsize == 0)
        mrsize = 1;

    com_ref_tot_dx[XX] /= mrsize;
    com_ref_tot_dx[YY] /= mrsize;
    com_ref_tot_dx[ZZ] /= mrsize;

    com_ref_tot[XX] += com_ref_tot_dx[XX];
    com_ref_tot[YY] += com_ref_tot_dx[YY];
    com_ref_tot[ZZ] += com_ref_tot_dx[ZZ];

    // for the first frame, calc absolute coms
    if (frame_index==0)
    {
        com_meas_tot[XX] = com_meas_tot[YY] = com_meas_tot[ZZ] = 0.;

        for(std::deque<mol>::iterator i=mols_meas.begin();i!=mols_meas.end();
            ++i)
        {
            rvec com;
            i->get_com(com);
            com_meas_tot[XX] += com[XX];
            com_meas_tot[YY] += com[YY];
            com_meas_tot[ZZ] += com[ZZ];
        }
        com_meas_tot[XX] /= mols_meas.size();
        com_meas_tot[YY] /= mols_meas.size();
        com_meas_tot[ZZ] /= mols_meas.size();

        com_ref_tot[XX] = com_ref_tot[YY] = com_ref_tot[ZZ] = 0.;

        if (mols_ref.size() > 0)
        {
            fprintf(stderr, "mols_ref.size=%ld\n", mols_ref.size());

            for(std::deque<mol>::iterator i=mols_ref.begin(); 
                i!=mols_ref.end(); 
                ++i)
            {
                /*if (i->get_N() != 3)
                    fprintf(stderr, "ERROR: size=%d\n", i->get_N());*/
                rvec com;
                i->get_com(com);
                com_ref_tot[XX] += com[XX];
                com_ref_tot[YY] += com[YY];
                com_ref_tot[ZZ] += com[ZZ];

                if (!std::isfinite(com[XX]))
                {
                    fprintf(stderr, 
                      "size: %d, molnr=%d, molblock_index=%d, mol_index=%d\n",
                      i->get_N(), i->get_molnr(), i->get_molblock_index(),
                      i->get_mol_index());
                }
            }

            com_ref_tot[XX] /= mols_ref.size();
            com_ref_tot[YY] /= mols_ref.size();
            com_ref_tot[ZZ] /= mols_ref.size();
        }
    }

    if (z_absolute)
    {
        com_ref_ev[XX] = com_ref_ev[YY] = com_ref_ev[ZZ] = 0;
        com_ref_ev_dx[XX] = com_ref_ev_dx[YY] = com_ref_ev_dx[ZZ] = 0;
    }
    else
    {
        com_ref_ev[XX]=com_ref_tot[XX]; 
        com_ref_ev[YY]=com_ref_tot[YY];
        com_ref_ev[ZZ]=com_ref_tot[ZZ];
        com_ref_ev_dx[XX]=com_ref_tot_dx[XX]; 
        com_ref_ev_dx[YY]=com_ref_tot_dx[YY];
        com_ref_ev_dx[ZZ]=com_ref_tot_dx[ZZ];
    }

    if (frame_index==0)
    {
        // calculate the z coord limits.
        proc_first_frame(pbc);
    }


    if (z_rho_hist_meas || z_rho_hist_ref)
        proc_z_rho_hist(pbc);


    if (zout)
    {
        fprintf(zout, "%g\t%g\t%g\t%g\n", time, 
                com_ref_tot[XX], com_ref_tot[YY], com_ref_tot[ZZ]);
    }
    if (xyzout)
    {
        //fprintf(xyzout, "%ld\n%g\n", mols_meas.size()+mols_ref.size(), time);
        fprintf(xyzout, "%ld\ntime = %g com = %g %g %g, com dx= %g %g %g\n", 
                mols_meas.size(), time, 
                com_meas_tot[XX], com_meas_tot[YY], com_meas_tot[ZZ],
                com_meas_tot_dx[XX], com_meas_tot_dx[YY], com_meas_tot_dx[ZZ]);
    }
    check_events(pbc, make_pers_events, make_exch_events, time);

    frame_index++;
}

real comsystem::get_closest_ref_dist_sq(mol &m, t_pbc *pbc)
{
    real min_dist_sq=1e10;
    rvec pos;
    m.get_com(pos);

    for(std::deque<mol>::iterator i=mols_ref.begin();i!=mols_ref.end();++i)
    {
        rvec dx;
        rvec rc;

        i->get_com(rc);
     
        pbc_dx(pbc, pos, rc, dx);

        real dist_sq = dx[XX]*dx[XX] + dx[YY]*dx[YY] + dx[ZZ]*dx[ZZ];
        if (dist_sq < min_dist_sq)
            min_dist_sq = dist_sq;
    }
    return min_dist_sq;
}

real comsystem::get_com_dist_sq(mol &m, t_pbc *pbc)
{
    rvec pos;
    m.get_com(pos);

    rvec dx;
    rvec com_ref_tot_r;
    com_ref_tot_r[XX] = com_ref_tot[XX];
    com_ref_tot_r[YY] = com_ref_tot[YY];
    com_ref_tot_r[ZZ] = com_ref_tot[ZZ];
    pbc_dx(pbc, pos, com_ref_tot_r, dx);

    return dx[XX]*dx[XX] + dx[YY]*dx[YY] + dx[ZZ]*dx[ZZ];
}

void comsystem::check_events(t_pbc *pbc,
                             bool make_pers_events,
                             bool make_exch_events, 
                             real time)
{
    int j=0,k=0;
    for(std::deque<mol>::iterator i=mols_meas.begin();i!=mols_meas.end();++i)
    {
        int ev=0;
        eventpos *ep;
        std::list<eventpos>::iterator  next=i->get_first();;
        while( (ep = i->check_events(next, mols_ref, pbc, com_meas_tot, 
                                     com_meas_tot_dx, jumpdist, weightdist,
                                     closest_n)) 
               != NULL)
        {
            switch(ep->get_type())
            {
                case persistence_event:
#if 0
                    printf("%d, %d, Persistence event: %g (N=%d)\n", 
                           frame_index, j, time-ep->get_time(), ep->get_N());
#endif
                    /* we subtract delta t/2 becasue otherwise we'd 
                       overestimate the time difference: we know exactly 
                       when the persistence time starts, but not when it stops,
                       only that is had stopped somewhere between the last
                       frames.*/
                    {
                        real dt= (time - ep->get_time()) - delta_t/2.;

                        if (event_relevant(ep))
                        {
                            persistence_events.add(dt, 
                                                   ep->get_refp(), 
                                                   ep->get_r(), 
                                                   ep->get_dr());
                        }
                    }
                    if (ev==0)
                        ev=1;
                    break;
                case exchange_event:
#if 0
                    printf("%d, %d, Exchange event: %g (N=%d)\n", 
                           frame_index, j, time-ep->get_time(), ep->get_N());
#endif
                    {
                        if (event_relevant(ep))
                        {
                            exchange_events.add(time - ep->get_time(), 
                                                ep->get_refp(),
                                                ep->get_r(), 
                                                ep->get_dr());
                        }
                    }
                    // NOTE: no break here
                case initial_exchange_event:
#if 0
                    printf("%d, %d: %s event: %g (N=%d, part=%d, r=%g)\n", 
                           frame_index, j,
                           (ep->get_type() == exchange_event ? 
                            "Exchange" : "Initial exchange" ),
                           time-ep->get_time(), ep->get_N(),j, ep->get_r());
#endif
                    if (make_exch_events)
                    {
                        real ref_dist_sq=0.;
                        bool keep = check_make_event(*i, pbc, ref_dist_sq);

                        if (keep)
                        {
                            i->make_event(mols_meas, mols_ref, pbc, time, 
                                          ref_trans,
                                          com_meas_tot, com_meas_tot_dx, 
                                          com_ref_ev, com_ref_ev_dx,
                                          jumpdist, weightdist, closest_n, 
                                          exchange_event, evdir,
                                          calc_r, ref_dist_sq);
                            Nevents_left++;
                        }
                    }
                    ev=2;
                    break;
                default:
                    printf("ERROR unknown event type\n");
                    exit(1);
            }
            Nevents_left--;
            Nanchor_av+=ep->get_N();
            anchor_d+=ep->get_anchor_d();
            NNanchor_av++;
            k++;
        }

        // in the first frame we should make initial_exchange_events
        if (frame_index==0)
        {
            real ref_dist_sq=0.;
            bool keep = check_make_event(*i, pbc, ref_dist_sq);

            if (keep)
            {
                i->make_event(mols_meas, mols_ref, pbc, time, ref_trans,
                              com_meas_tot, com_meas_tot_dx, com_ref_ev, 
                              com_ref_ev_dx,
                              jumpdist, weightdist, closest_n, 
                              initial_exchange_event, evdir,
                              calc_r, ref_dist_sq);
                Nevents_left++;
                if (draw_pdb)
                {
                    //i->make_initial_event(
                    i->startpos_event=new eventpos();
                    i->set_event(i->startpos_event,
                                 mols_meas, mols_ref, pbc, time, ref_trans,
                                 com_meas_tot, com_meas_tot_dx, com_ref_ev, 
                                 com_ref_ev_dx,
                                 jumpdist, weightdist, closest_n, 
                                 initial_exchange_event, evdir,
                                 calc_r, ref_dist_sq);
                }
            }
        }
        if (make_pers_events || frame_index == 0)
        {
            real ref_dist_sq=0.;
            bool keep = check_make_event(*i, pbc, ref_dist_sq);

            if (keep)
            {
                i->make_event(mols_meas, mols_ref, pbc, time, ref_trans,
                              com_meas_tot, com_meas_tot_dx, com_ref_ev, 
                              com_ref_ev_dx,
                              jumpdist, weightdist, closest_n, 
                              persistence_event, 
                              evdir,
                              calc_r, ref_dist_sq);
                Nevents_left++;
            }
        }
        j++;

        if (xyzout)
        {
            rvec v;
            rvec dx;
            i->get_com(v);
            i->get_dx(dx);
            
            char evs='C';
            if (ev==1)
                evs='N';
            else if (ev==2)
                evs='O';

            fprintf(xyzout, "%c  %g %g %g   %g %g %g\n", evs, 
                    v[XX], v[YY], v[ZZ], dx[XX], dx[YY], dx[ZZ]);
        }
    }
}

void comsystem::write_ratio(FILE *outf, 
                            const weighted_hist &pers_hist, 
                            const weighted_hist &exch_hist)
{
    for (size_t i=0; i<pers_hist.get_N(); ++i)
    {
        if (pers_hist.have_err(i) && exch_hist.have_err(i) && 
            pers_hist.get(i) > 1e-6 && exch_hist.get(i) > 1e-6)
        {
            double x = (pers_hist.get_x(i) + pers_hist.get_x(i+1))/2.;
            double y = (pers_hist.get(i) / exch_hist.get(i) );
            double relerra = pers_hist.get_err(i) / pers_hist.get(i);
            double relerrb = exch_hist.get_err(i) / exch_hist.get(i);
            double yerr = sqrt(relerra*relerra + relerrb*relerrb);
            /*printf("p=%g, e=%g, relerra=%g, relerrb=%g, yerr=%g\n", 
                   pers_hist.get(i), exch_hist.get(i),
                   relerra, relerrb, yerr);*/
            fprintf(outf, "%g %g %g\n", x, y, yerr);
        }
    }
}

void comsystem::proc_pers_rho(const weighted_hist &z_pers, 
                              const hist &zrho,    
                              weighted_hist &z_rho_pers) 
{
    if (0)
    {
    }
}

void comsystem::proc_inv_inv(const weighted_hist &inv,
                             weighted_hist &inv_inv)
{
    if (inv.get_N() != inv_inv.get_N())
    {
        fprintf(stderr, "ERROR inv and inv_inv unequal\n");
        exit(1);
    }
    for(unsigned int i=0; i<inv.get_N(); ++i)
    {
        if (inv.have_err(i))
        {
            double vo=inv.get(i);
            double so=inv.get_err(i);
            /* sigma_new^2 = v_new^2 * ( sigma_old/v_old  )^2 
               i.e. 
               sigma_new^2 = sigma_old^2 / v_old^4 
            */
            inv_inv.set(i, 1./vo, sqrt( so*so/(vo*vo*vo*vo) ));
        }
    }
}

void comsystem::process_events(FILE *logf, size_t Nbin, 
                               const char *pers_exch_histfn, 
                               const char *zhistfn, 
                               const char *gridfn, 
                               const char *rhist, 
                               const char *drhist, 
                               const output_env_t oenv)
{
    int pers=0,exch=0,initexch=0;
    for(std::deque<mol>::iterator i=mols_meas.begin();i!=mols_meas.end();++i)
    {
        i->count_events(pers,exch,initexch);
    }
    fprintf(logf,"%5ld meas molecules\n", mols_meas.size());
    fprintf(logf,"%5ld ref molecules\n", mols_ref.size());
    fprintf(logf,"Remaining pre-events:\n");
    fprintf(logf,"%5d        total pre-events\n", Nevents_left);
    fprintf(logf,"%5d  persistence pre-events\n", pers);
    fprintf(logf,"%5d     exchange pre-events\n", exch);
    fprintf(logf,"%5d pre-exchange pre-events\n", initexch);
    fprintf(logf,"\n");
    fprintf(logf,"                       Jump distance: %g\n", jumpdist);
    fprintf(logf,"               Weight sigma distance: %g\n", weightdist);
    fprintf(logf,"                    Weight closest N: %d\n", closest_n);
    switch(evdir)
    {
        case evdir_xyz:
        default:
            fprintf(logf,"                          Distances : 3D (xyz)\n");
            break;
        case evdir_xy:
            fprintf(logf,"                          Distances : 2D (xy)\n");
            break;
        case evdir_z:
            fprintf(logf,"                          Distances : 1D (z)\n");
            break;
    }
    fprintf(logf,"           Average number of anchors: %g\n", 
            Nanchor_av/NNanchor_av);
    fprintf(logf,"Average weighted distance of anchors: %g\n", 
            anchor_d/NNanchor_av);
    fprintf(logf,"\n");


    double delta_t = (last_t-first_t)/(frame_index-1);

    fprintf(logf,"\n");
    fprintf(logf,"      t: %g - %g\n", first_t, last_t);
    fprintf(logf," Frames: %d\n", frame_index);
    fprintf(logf,"Delta t: %g\n", delta_t);
    fprintf(logf,"\n");


    fprintf(logf,"\n");
    fprintf(logf, "N persistence events: %ld\n", persistence_events.get_N());
    fprintf(logf, "N    exchange events: %ld\n", exchange_events.get_N());
    fprintf(logf,"\n");

    if (persistence_events.get_N() == 0 || exchange_events.get_N()==0)
    {
        fprintf(logf, "ERROR: no events!\n");
        return;
    }

    /*real rmin=exchange_events.get_min_r(); 
    real rmax=exchange_events.get_max_r(); */
    real rmin=0.;
    real rmax_min=min_r_found[XX]*min_r_found[XX] + 
                  min_r_found[YY]*min_r_found[YY] +
                  min_r_found[ZZ]*min_r_found[ZZ];
    real rmax_max=max_r_found[XX]*max_r_found[XX] + 
                  max_r_found[YY]*max_r_found[YY] +
                  max_r_found[ZZ]*max_r_found[ZZ];

    real rmax;
    if (rmax_min > rmax_max)
        rmax=sqrt(rmax_min);
    else
        rmax=sqrt(rmax_max);

    int Nrbin=(rmax-rmin)/r_spacing;

    real drmin=exchange_events.get_min_dr();
    if (drmin>persistence_events.get_min_dr())
        drmin=persistence_events.get_min_dr();
    real drmax=exchange_events.get_max_dr();
    if (drmax<persistence_events.get_max_dr())
        drmax=persistence_events.get_max_dr();
    int Ndrbin=(drmax-drmin)/dr_spacing;



    hist exch_hist(Nbin, log(delta_t), log(last_t-first_t) );
    hist pers_hist(Nbin, log(delta_t), log(last_t-first_t) );


    // now determine spacing in z
    double min_z_hist=min_r_found[ZZ];
    double max_z_hist=max_r_found[ZZ];

    if (z_min_meas>1e-6)
        min_z_hist=z_min_meas;
    if (check_z_min_max)
        max_z_hist=z_max_meas;

    if (!zhistfn) //min_z >= max_z)
    {
        Nzbin = 100;
        min_z_hist=0;
        max_z_hist=1;
    }

    weighted_hist z_whist_exch(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_whist_pers(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_inv_whist_exch(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_inv_whist_pers(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_inv_inv_whist_exch(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_inv_inv_whist_pers(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_N_pers(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_N_exch(Nzbin, min_z_hist, max_z_hist);
    weighted_hist z_rho_pers(Nzbin, min_z_hist, max_z_hist);

    
    hist *r_hist_exch=NULL;
    hist *r_hist_pers=NULL;
    weighted_hist *r_whist_exch=NULL;
    weighted_hist *r_whist_pers=NULL;
    weighted_hist *r_inv_whist_exch=NULL;
    weighted_hist *r_inv_whist_pers=NULL;
    weighted_hist *r_inv_inv_whist_exch=NULL;
    weighted_hist *r_inv_inv_whist_pers=NULL;
   
    if (calc_r) 
    {
        r_hist_exch=new hist(Nrbin, rmin, rmax);
        r_hist_pers=new hist(Nrbin, rmin, rmax);
        r_whist_exch=new weighted_hist(Nrbin, rmin, rmax);
        r_whist_pers=new weighted_hist(Nrbin, rmin, rmax);
        r_inv_whist_exch=new weighted_hist(Nrbin, rmin, rmax);
        r_inv_whist_pers=new weighted_hist(Nrbin, rmin, rmax);
        r_inv_inv_whist_exch=new weighted_hist(Nrbin, rmin, rmax);
        r_inv_inv_whist_pers=new weighted_hist(Nrbin, rmin, rmax);
    }

    hist dr_hist_exch(Ndrbin, drmin, drmax);
    hist dr_hist_pers(Ndrbin, drmin, drmax);


    double exch_av;
    double exch_av_err;
    double exch_av_inv;
    double exch_av_inv_err;
    int exch_N;

    double pers_av;
    double pers_av_err;
    double pers_av_inv;
    double pers_av_inv_err;
    int pers_N;

    const int startMoment=-2;
    const int Nmoments=6;
    double exch_moments[Nmoments], exch_moments_err[Nmoments];
    double pers_moments[Nmoments], pers_moments_err[Nmoments];

    exchange_events.proc(exch_hist, 
                         z_whist_exch, 
                         z_inv_whist_exch, 
                         &z_N_exch,
                         r_hist_exch, r_whist_exch, r_inv_whist_exch, 
                         dr_hist_exch,
                         exch_N, exch_av, exch_av_err,
                         exch_av_inv, exch_av_inv_err,
                         startMoment, Nmoments, 
                         exch_moments, exch_moments_err, exch_grid);
    persistence_events.proc(pers_hist, 
                            z_whist_pers, 
                            z_inv_whist_pers, 
                            &z_N_pers,
                            r_hist_pers, r_whist_pers, r_inv_whist_pers,
                            dr_hist_pers,
                            pers_N, pers_av, pers_av_err,
                            pers_av_inv, pers_av_inv_err,
                            startMoment, Nmoments, 
                            pers_moments, pers_moments_err, pers_grid);

    if (gridfn && exch_grid && pers_grid)
    {
        fprintf(logf, "Writing grid data to %s\n\n", gridfn);
        write_grids(gridfn, ref_grid, pers_grid, exch_grid);
    }

    if (zhistfn)
    {
        proc_pers_rho(z_whist_pers, *z_rho_hist_meas, z_rho_pers);
    }
    proc_inv_inv(z_inv_whist_exch, z_inv_inv_whist_exch);
    proc_inv_inv(z_inv_whist_pers, z_inv_inv_whist_pers);
    if (r_inv_inv_whist_exch)
        proc_inv_inv(*r_inv_whist_exch, *r_inv_inv_whist_exch);
    if (r_inv_inv_whist_pers)
        proc_inv_inv(*r_inv_whist_pers, *r_inv_inv_whist_pers);

    fprintf(logf, "Moment            exch                           pers\n");
    for(int i=0;i<Nmoments;i++)
    {
        fprintf(logf,"moment %d    %10g +/- %10g      %10g +/- %10g\n", 
                i+startMoment, 
                exch_moments[i], exch_moments_err[i],
                pers_moments[i], pers_moments_err[i]);
    }
    fprintf(logf, "\n");

    {
        double est=exch_moments[2]/(2*exch_moments[1]);

        double est1=exch_moments_err[1]/exch_moments[1];
        double est2=exch_moments_err[2]/exch_moments[2];

        //double est_err=sqrt( est1*est1 + est2*est2 )*est;
        double est_err=sqrt(sqr(exch_moments_err[1])*
                            sqr(exch_moments[2]/
                                (2*exch_moments[1]*exch_moments[1]))
                            +
                            sqr(exch_moments_err[2])*
                            sqr(1./(2.*exch_moments[1]) ) );
        fprintf(logf, "Estimate for <tp> = %g +/- %g\n", est, est_err);
    }

    fprintf(logf,"\n");
    fprintf(logf, "   Exchange times: %d events, average time: %g +/- %g\n",
            exch_N, exch_av, exch_av_err);
    fprintf(logf, "Persistence times: %d events, average time: %g +/- %g\n",
            pers_N, pers_av, pers_av_err);

    double tot_err = sqrt( (pers_av_err*pers_av_err / (pers_av*pers_av)) +
                           (exch_av_err*exch_av_err / (exch_av*exch_av)) 
                         );
    fprintf(logf, "Ratio: %g +/- %g\n", pers_av/exch_av, tot_err);
    fprintf(logf, "\n");
    fprintf(logf, "Inverse: <1/t_x> = %g +/- %g , <1/t_p> = %g +/- %g\n",
            exch_av_inv, exch_av_inv_err,
            pers_av_inv, pers_av_inv_err);
    double tot_err_inv = sqrt( (pers_av_inv_err*pers_av_inv_err / 
                                    (pers_av_inv*pers_av_inv)) +
                               (exch_av_inv_err*exch_av_inv_err / 
                                    (exch_av_inv*exch_av_inv)) 
                         );
    fprintf(logf, "Inverse ratio: %g +/- %g\n", 
            exch_av_inv/pers_av_inv, tot_err_inv);
    fprintf(logf,"\n");



    if (pers_exch_histfn)
    {
        const char *sets[] = { "N_pers", "N_exch", NULL};
        fprintf(logf, "Writing %s\n", pers_exch_histfn);
        FILE *out=xvgropen_type(pers_exch_histfn, "Persistence/exchange time", 
                                "ln(t)", "N", exvggtXNY, oenv);
        xvgr_legend(out, 2, sets, oenv);
        pers_hist.writef(out, false);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        exch_hist.writef(out, false);
        fclose(out);
    }
    //fprintf(logf, "Writing exchange.hist\n");
    //exch_hist.write("exchange.hist");

    if (zhistfn)
    {
        const char *sets[] = { "density", "<t_pers>", "<t_exch>", 
                               "<1/t_pers>", "<1/t_exch>",
                               "1/<1/t_pers>", "1/<1/t_exch>",
                               "N_pers", "N_exch", 
                               "ratio",
                               z_absolute ?  NULL : "dens_ref",
                               NULL};
        int Nsets=0;
        // count the number of sets
        for(int i=0;i<100;i++)
        {
            if (sets[i] == NULL)
            {
                Nsets=i+1;
                break;
            }
        }
        fprintf(logf, "Writing %s\n", zhistfn);
        FILE *out=xvgropen_type(zhistfn, "Measurements in z",
                                "z", "", exvggtXYDY, oenv);
        xvgr_legend(out, Nsets, sets, oenv);
        //z_whist_ratio.writef(out, false, false);
        //xvgr_new_dataset(out, 0, 0, NULL, oenv);
        z_rho_hist_meas->writef(out, true, false, true);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);

        z_whist_pers.writef(out, false, false, true);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        z_whist_exch.writef(out, false, false, true);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        z_inv_whist_pers.writef(out, false, false, true);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        z_inv_whist_exch.writef(out, false, false, true);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);

        z_inv_inv_whist_pers.writef(out, false, false, true);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        z_inv_inv_whist_exch.writef(out, false, false, true);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);

        z_N_pers.writef(out, false, false, true, false);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        z_N_exch.writef(out, false, false, true, false);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);

        write_ratio(out, z_whist_pers, z_whist_exch);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);

        if (!z_absolute)
        {
            z_rho_hist_ref->writef(out,true, true);
        }
        /*z_hist_pers.writef(out);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        z_hist_exch.writef(out);*/

        fclose(out);
    }
    if (rhist && calc_r)
    {
        const char *sets[] = { "N_pers", "N_exch", "<t_pers>", "<t_exch>",
                               "<1/t_pers>", "<1/t_exch>",
                               "1/<1/t_pers>", "1/<1/t_exch>",
                               NULL };
        fprintf(logf, "Writing %s\n", rhist);
        FILE *out=xvgropen_type(rhist, 
                                //"Exchange times in nearest particle dist",
                                "Persistence/exchange times in ref. com dist",
                                "r (nm)", "N", exvggtXNY, oenv);
        xvgr_legend(out, 8, sets, oenv);
        r_hist_pers->writef(out);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        r_hist_exch->writef(out);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        r_whist_pers->writef(out, false, false);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        r_whist_exch->writef(out, false, false);

        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        r_inv_whist_pers->writef(out, false, false);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        r_inv_whist_exch->writef(out, false, false);

        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        r_inv_inv_whist_pers->writef(out, false, false);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        r_inv_inv_whist_exch->writef(out, false, false);

        fclose(out);
    }
    if (drhist)
    {
        const char *sets[] = { "N (pers)", "N (exch)", NULL };
        fprintf(logf, "Writing %s\n", drhist);
        FILE *out=xvgropen_type(drhist, 
                                "Jump distance",
                                "r (nm)", "N", exvggtXNY, oenv);
        xvgr_legend(out, 2, sets, oenv);
        dr_hist_pers.writef(out, true, false);
        xvgr_new_dataset(out, 0, 0, NULL, oenv);
        dr_hist_exch.writef(out, true, false);
        fclose(out);
    }

    /*if (exch_grid)
    {
        delete exch_grid;
    }
    if (pers_grid)
    {
        delete pers_grid;
    }*/
}

void comsystem::proc_first_frame(t_pbc *pbc)
{
    //min_z_found=0; max_z_found=0;

    clear_rvec(min_r_found);
    clear_rvec(max_r_found);


    for(std::deque<mol>::iterator i=mols_meas.begin();i!=mols_meas.end(); ++i)
    {
        rvec dxvec;
        rvec r;

        i->get_com(r);
        pbc_dx(pbc, r, com_ref_ev, dxvec);

        for(int j=0;j<DIM;j++)
        {
            if (dxvec[j] < min_r_found[j])
               min_r_found[j] = dxvec[j]; 
            if (dxvec[j] > max_r_found[j])
               max_r_found[j] = dxvec[j]; 
        }

        /*real z_val = dxvec[ZZ];

        if (z_val < min_z_found)
            min_z_found = z_val;
        if (z_val > max_z_found)
            max_z_found = z_val;*/
    }
    for(std::deque<mol>::iterator i=mols_ref.begin();i!=mols_ref.end(); ++i)
    {
        rvec dxvec;
        rvec r;

        i->get_com(r);
        pbc_dx(pbc, r, com_ref_ev, dxvec);

        for(int j=0;j<DIM;j++)
        {
            if (dxvec[j] < min_r_found[j])
               min_r_found[j] = dxvec[j]; 
            if (dxvec[j] > max_r_found[j])
               max_r_found[j] = dxvec[j]; 
        }
    }

    if (max_r_found[ZZ] > 0)
        max_r_found[ZZ] += z_slack*z_spacing;
    if (min_r_found[ZZ] < 0)
        min_r_found[ZZ] -= z_slack*z_spacing;
    Nzbin = (int)((max_r_found[ZZ] - min_r_found[ZZ])/z_spacing);

    // correct for z_min/max_meas 
    if (min_r_found[ZZ] < z_min_meas)
        min_r_found[ZZ]=z_min_meas;
    if (max_r_found[ZZ] > z_max_meas)
        max_r_found[ZZ]=z_max_meas;

    if (do_z_rho_hist)
    {
        z_rho_hist_meas = new hist(Nzbin, min_r_found[ZZ], max_r_found[ZZ]);
        z_rho_hist_ref = new hist(Nzbin, min_r_found[ZZ], max_r_found[ZZ]);
    }

    //grid *exch_grid=NULL;
    //grid *pers_grid=NULL;

    if (grid_spacing > 0)
    {
        fprintf(stdout, "Grid spacing: %g\n", grid_spacing);
        int nx,ny,nz;
        unsigned long int ntot;
        nx=(int)((max_r_found[XX] - min_r_found[XX])/grid_spacing);
        ny=(int)((max_r_found[YY] - min_r_found[YY])/grid_spacing);
        nz=(int)((max_r_found[ZZ] - min_r_found[ZZ])/grid_spacing);
        ntot = nx;
        ntot *= ny;
        ntot *= nz;
#if 0
        fprintf(stdout, "Allocating (%d, %d, %d)=%gM grid locations.\n\n",
                nx,ny,nz, ntot/(1024.*1024.));
        fflush(stdout);
#endif
        exch_grid = new grid(grid_spacing, min_r_found, max_r_found);
        pers_grid = new grid(grid_spacing, min_r_found, max_r_found);
        ref_grid = new grid(grid_spacing, min_r_found, max_r_found);
    }
}


void comsystem::proc_z_rho_hist(t_pbc *pbc)
{
    for(std::deque<mol>::iterator i=mols_meas.begin();i!=mols_meas.end(); ++i)
    {
        rvec dxvec;
        rvec r;

        i->get_com(r);
        pbc_dx(pbc, r, com_ref_ev, dxvec);
        real z_val = dxvec[ZZ];

        z_rho_hist_meas->add(z_val);
    }
    for(std::deque<mol>::iterator i=mols_ref.begin();i!=mols_ref.end(); ++i)
    {
        rvec dxvec;
        rvec r;

        i->get_com(r);
        pbc_dx(pbc, r, com_ref_ev, dxvec);
        real z_val = dxvec[ZZ];

        z_rho_hist_ref->add(z_val);
    }
}
#if 0
void comsystem::write_pdb_x(const char *filename, int natoms, rvec x[])
{
    //printf("Opening %s for PDB output\n", filename);
    FILE *pdbout=fopen(filename, "w");

    fprintf(pdbout, "REMARK    GENERATED BY TIMEPROC\n");
    fprintf(pdbout, "TITLE     SMP-ensv-03 t=   0.00000\n");
    fprintf(pdbout, "REMARK    THIS IS A SIMULATION BOX\n");
    //fprintf(out, "CRYST1   46.982   46.982   46.982  60.00  60.00  90.00 P 1           1\n");
    fprintf(pdbout, "MODEL        1\n");

    /*"ATOM  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n";
             ^ atom number
                 ^ atom name
                     ^ residue name
                          ^ residue number
                                 ^ xyz coordinates
                                                ^ occupancy
                                                     ^ temperature factor
    */

    real mind=1e10;
    real maxd=-1;

    for(int i=0;i<natoms;i++)
    {
        real d=atom_mols[i]->startpos_event->get_dist();
        fprintf(pdbout,
                "ATOM  %5i %4s %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",
                i,
                atom_names[i],
                res_names[i],
                res_nrs[i],
                start_coords_x[i]*10,
                start_coords_y[i]*10,
                start_coords_z[i]*10,
                1., d );
        if (d<mind)
            mind=d;
        if (d>maxd)
            maxd=d;
    }

    printf("\nmin d=%g maxd = %g\n", mind, maxd);
    fprintf(pdbout, "REMARK MIND=%g MAXD=%g\n", mind, maxd);
    fclose(pdbout);
}
#endif

void comsystem::write_pdb(const char *filename, t_pbc *pbc)
{
    //printf("Opening %s for PDB output\n", filename);
    FILE *pdbout=fopen(filename, "w");

    fprintf(pdbout, "REMARK    GENERATED BY TIMEPROC\n");
    fprintf(pdbout, "TITLE     SMP-ensv-03 t=   0.00000\n");
    fprintf(pdbout, "REMARK    THIS IS A SIMULATION BOX\n");
    //fprintf(out, "CRYST1   46.982   46.982   46.982  60.00  60.00  90.00 P 1           1\n");
    fprintf(pdbout, "MODEL        1\n");

    /*"ATOM  %5i %4s %3s  %4s    %8.3f%8.3f%8.3f%6.2f%6.2f    \n";
             ^ atom number
                 ^ atom name
                     ^ residue name
                          ^ residue number
                                 ^ xyz coordinates
                                                ^ occupancy
                                                     ^ temperature factor
    */

    real mind=1e10;
    real maxd=-1;

    for(int i=0;i<Natoms;i++)
    {
        //real d=atom_mols[i]->startpos_event->get_dist();
        real d2,d;
        mol *m=&(*((atom_mols[i])));
        if (m)
        {
            d=check_make_event(*m, pbc, d2);
        }
        else
        {
            d=0;
        }


        fprintf(pdbout,
                "ATOM  %5i %4s %3s  %4i    %8.3f%8.3f%8.3f%6.2f%6.2f    \n",
                i,
                atom_names[i],
                res_names[i],
                res_nrs[i],
                start_coords_x[i]*10,
                start_coords_y[i]*10,
                start_coords_z[i]*10,
                1., d );
        if (d<mind)
            mind=d;
        if (d>maxd)
            maxd=d;
    }

    printf("\nmin d=%g maxd = %g\n", mind, maxd);
    fprintf(pdbout, "REMARK MIND=%g MAXD=%g\n", mind, maxd);
    fclose(pdbout);


}

