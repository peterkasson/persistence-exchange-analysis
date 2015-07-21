
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include <vector>
#include <list>
#include <deque>

#include <gromacs/typedefs.h>
#include <gromacs/pbc.h>
#include <gromacs/vec.h>

using namespace std;

#include "hist.h"
#include "event.h"
#include "grid.h"
#include "proccom.h"



void eventpos::reinit(real time_, rvec r, 
                      const deque<mol> &meas_mols, 
                      const deque<mol> &ref_mols, t_pbc *pbc,
                      matrix ref_trans,
                      rvec com_meas_tot, rvec com_meas_tot_dx,
                      rvec com_ref_tot, rvec com_ref_tot_dx,
                      real weight_sigma, int closest_n, event_type type_,
                      event_direction evdir_,
                      bool calc_r_, std::list<anchor> &spare_anchor_list,
                      real ref_dist_sq_)
{
    event_happened=false;
    checked=true; // because we don't want cumulative dxes to start at t=0
    time=time_;
    type=type_;

    ref_dist_sq = ref_dist_sq_;

    // remove all the existing anchors
    spare_anchor_list.splice(spare_anchor_list.begin(), anchors);

    if (weight_sigma > 0)
    {
        weight_sigmasq = weight_sigma*weight_sigma;
        use_weights=true;
        if (closest_n > 0)
        {
            gmx_fatal(FARGS, "ERRROR: both closest_n >0 and weight_sigma >0");
        }
    }
    else if (closest_n > 0)
    {
        use_weights=true;
    }
    else
    {
        weight_sigmasq=0;
        use_weights=false;
    }
        
    weight_norm = 1./sqrt(2.*M_PI*weight_sigmasq);

    real com_weight=0.;

    evdir=evdir_;

    anchor_d=0.;

    if (use_weights)
    {
        if (closest_n > 0)
        {
            const mol* closest[closest_n];
            double dsq_max[closest_n];
            for(int i=0;i<closest_n;i++)
            {
                closest[i] = NULL;
                dsq_max[i] = 1e20;
            }
            int k=0;

            // now find the closest n mols
            for(deque<mol>::const_iterator i=ref_mols.begin(); 
                i!=ref_mols.end(); ++i)
            {
                rvec dxvec;
                rvec lcom;
                i->get_com(lcom);

                // calculate distance vector
                pbc_dx(pbc, r, lcom, dxvec); 
                // and calculate distance to base weights on
                real dsq;
                switch(evdir)
                {
                    case evdir_xyz:
                    default:
                        dsq = (dxvec[XX]*dxvec[XX] + dxvec[YY]*dxvec[YY] + 
                               dxvec[ZZ]*dxvec[ZZ]);
                        break;
                    case evdir_xy:
                        dsq = dxvec[XX]*dxvec[XX] + dxvec[YY]*dxvec[YY];
                        break;
                    case evdir_z:
                        dsq = dxvec[ZZ]*dxvec[ZZ];
                        break;
                }
                k++;
                /* check whether it's closer than any before */
                if (dsq < dsq_max[closest_n-1])
                {
                    int j;
                    for(j=closest_n - 1; j >= 0; j--)
                    {
                        if (dsq < dsq_max[j])
                        {
                            // shift everything one place up
                            if (j < closest_n - 1)
                            {
                                dsq_max[j+1] = dsq_max[j];
                                closest[j+1] = closest[j];
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                    // we always overshoot by 1
                    j++;
                    dsq_max[j] = dsq;
                    closest[j] = &(*i);
                }
            }
            for(int i=0;i<closest_n;i++)
            {
                /*printf("closest[%d], dsq=%g, p=%p\n", i, sqrt(dsq_max[i]),
                       closest[i]);
                fflush(stdout);*/
                // add the anchors
                if (spare_anchor_list.empty())
                {
                    for(int i=0;i<1000;i++)
                    {
                        spare_anchor_list.push_back(anchor());
                    }
                }

                rvec lcom;
                closest[i]->get_com(lcom);
                std::list<anchor>::iterator newelem=spare_anchor_list.begin();
                newelem->reinit(closest[i], 1., lcom);
                anchors.splice(anchors.begin(), spare_anchor_list, newelem);
                com_weight += 1.;
                anchor_d += sqrt(dsq_max[i]);
            }
            anchor_d /= com_weight;
        }
        else
        {
            for(deque<mol>::const_iterator i=ref_mols.begin(); 
                i!=ref_mols.end(); ++i)
            {
                rvec dxvec;
                rvec lcom;
                i->get_com(lcom);

                // calculate distance vector
                pbc_dx(pbc, r, lcom, dxvec); 
                // and calculate distance to base weights on
                real dsq;
                switch(evdir)
                {
                    case evdir_xyz:
                    default:
                        dsq = (dxvec[XX]*dxvec[XX] + dxvec[YY]*dxvec[YY] + 
                               dxvec[ZZ]*dxvec[ZZ]);
                        break;
                    case evdir_xy:
                        dsq = dxvec[XX]*dxvec[XX] + dxvec[YY]*dxvec[YY];
                        break;
                    case evdir_z:
                        dsq = dxvec[ZZ]*dxvec[ZZ];
                        break;
                }

                real weight=weightfn( dsq );

                if (weight > 0.01)
                {
                    if (spare_anchor_list.empty())
                    {
                        for(int i=0;i<1000;i++)
                        {
                            spare_anchor_list.push_back(anchor());
                        }
                    }
                    // add the anchors
                    std::list<anchor>::iterator newelem=spare_anchor_list.begin();
                    newelem->reinit(&(*i), weight, lcom);
                    anchors.splice(anchors.begin(), spare_anchor_list, newelem);
                    //anchors.push_back(anchor(&(*i), weight, lcom));
                    com_weight += weight;
                    anchor_d += sqrt(dsq) * weight;
                }
            }
            anchor_d /= com_weight;
        }
        if (com_weight == 0.)
        {
            /*printf("ERRRORRREEEERROOORRR no weight\n");*/
            printf("WARNING no weight. Using total  COM as ref point\n");
            com_weight=1.;
        }
    }
    else
    {
        com_weight=1.;
    }

    { // calculate position relative to reference pos.
        rvec dxvec;

        pbc_dx(pbc, r, com_ref_tot, dxvec); 
        // now do reference transformation
        mvmul(ref_trans, dxvec, refp);

        //z_val_start=dxvec[ZZ];
        real_start_pos[XX] = r[XX];
        real_start_pos[YY] = r[YY];
        real_start_pos[ZZ] = r[ZZ];
    }

    // calculate r_val
    if (calc_r_)
    {
        r_val=0.;
#if 0
        r_val=FLT_MAX;
        bool found0=false;
        for(list<mol>::iterator i=meas_mols.begin(); i!=meas_mols.end(); ++i)
        {
            // calcualte the closest distance > 0
            rvec dxvec;
            rvec lcom;
            i->get_com(lcom);
            if (lcom[XX] != r[XX] || lcom[YY]!=r[YY] || lcom[ZZ]!=r[ZZ] )
            {

                // calculate distance vector
                pbc_dx(pbc, r, lcom, dxvec); 
                // and calculate distance to base weights on
                real dsq;
                switch(evdir)
                {
                    case evdir_xyz:
                    default:
                        dsq = (dxvec[XX]*dxvec[XX] + dxvec[YY]*dxvec[YY] + 
                               dxvec[ZZ]*dxvec[ZZ]);
                        break;
                    case evdir_xy:
                        dsq = dxvec[XX]*dxvec[XX] + dxvec[YY]*dxvec[YY];
                        break;
                    case evdir_z:
                        dsq = dxvec[ZZ]*dxvec[ZZ];
                        break;
                }
                if (dsq > 0 && dsq < r_val)
                {
#ifdef GMX_DOUBLE
                    r_val = sqrt(dsq);
#else
                    r_val = sqrtf(dsq);
#endif
                }
            }
            else
                found0=true;
        }
        //printf("r_val = %g\n", r_val);
        if (!found0)
        {
            printf("0 not found!\n");
            exit(1);
        }
        if (r_val > 1e5)
        {
            printf("r_val=%g!\n",r_val);
            exit(1);
        }
#endif
    }
    else
    {
        r_val=0.;
    }

    com_dx[XX] = com_dx[YY] = com_dx[ZZ] = 0;
    dx[XX] = dx[YY] = dx[ZZ] = 0;
#if 0
    printf(" r=(%g,%g,%g)\n", r[XX], r[YY], r[ZZ]);
#endif
}

void eventpos::write_all_pbc(FILE *out, t_pbc *pbc, rvec w)
{
    fprintf(out, "%ld\n \n", anchors.size());
    for(list<anchor>::iterator i=anchors.begin(); i!=anchors.end(); ++i)
    {
        rvec dxvec;
        rvec lcom;

        i->get_molp()->get_com(lcom); // get the anchor's com
        pbc_dx(pbc, w, lcom, dxvec);
#if 1
        dxvec[XX] += w[XX];
        dxvec[YY] += w[YY];
        dxvec[ZZ] += w[ZZ];
#endif
        fprintf(out, "C\t%g\t%g\t%g\n", dxvec[XX], dxvec[YY], dxvec[ZZ]);
    }
}

void eventpos::write_all_nopbc(FILE *out, rvec w, t_pbc *pbc)
{
    fprintf(out, "%ld\n \n", anchors.size());
    for(list<anchor>::iterator i=anchors.begin(); i!=anchors.end(); ++i)
    {
        rvec lcom;

        i->get_molp()->get_com(lcom); // get the anchor's com

        fprintf(out, "C\t%g\t%g\t%g\n", lcom[XX], lcom[YY], lcom[ZZ]);
    }
}


bool eventpos::check_event(rvec rdx, rvec com_meas_tot_dx, real event_dist)
{
    if (checked)
        return event_happened;

    // first calculate the change in position
    dx[XX] += rdx[XX];
    dx[YY] += rdx[YY];
    dx[ZZ] += rdx[ZZ];

    if (use_weights && anchors.size()>0 )
    {    
        rvec com_dx_sum;

        real ncom_weight=0.;

        com_dx_sum[XX]=0;
        com_dx_sum[YY]=0;
        com_dx_sum[ZZ]=0;

        /* loop over the anchors and calculate their delta x w.r.t. 
            the event's starting time */
        for(list<anchor>::iterator i=anchors.begin(); i!=anchors.end(); ++i)
        {
            rvec dx;

            real weight=i->get_weight();
            i->get_molp()->get_dx(dx);

            com_dx_sum[XX] += dx[XX]*weight;
            com_dx_sum[YY] += dx[YY]*weight;
            com_dx_sum[ZZ] += dx[ZZ]*weight;

            ncom_weight += weight;
        }
        // add the difference in com
        com_dx[XX] += com_dx_sum[XX]/ncom_weight;
        com_dx[YY] += com_dx_sum[YY]/ncom_weight;
        com_dx[ZZ] += com_dx_sum[ZZ]/ncom_weight;
    }
    else
    {
        // add the difference in com.
        com_dx[XX] += com_meas_tot_dx[XX];
        com_dx[YY] += com_meas_tot_dx[YY];
        com_dx[ZZ] += com_meas_tot_dx[ZZ];

        //pbc_dx(pbc, r, startpos, distvec);
    }

    rvec distvec;
    // calcualte the distance as the position difference minus the com movement
    distvec[XX] = dx[XX] - com_dx[XX];
    distvec[YY] = dx[YY] - com_dx[YY];
    distvec[ZZ] = dx[ZZ] - com_dx[ZZ];

    // only in the actual distance calculation is
    // the direction setting taken into effect.
    switch(evdir)
    {
        case evdir_xyz:
        default:
            dsq = (distvec[XX]*distvec[XX] + distvec[YY]*distvec[YY] +
                   distvec[ZZ]*distvec[ZZ]);
            break;
        case evdir_xy:
            dsq = distvec[XX]*distvec[XX] + distvec[YY]*distvec[YY];
            break;
        case evdir_z:
            dsq = distvec[ZZ]*distvec[ZZ];
            break;
    }
    /*dsq = distvec[XX]*distvec[XX] + distvec[YY]*distvec[YY];
    if (!do2d)
        dsq += distvec[ZZ]*distvec[ZZ];*/


    if (dsq > event_dist*event_dist)
    {
#ifdef GMX_DOUBLE
        dr_val = sqrt(dsq);
#else
        dr_val = sqrtf(dsq);
#endif
        event_happened=true;
        return true;
    }
    return false;
}



