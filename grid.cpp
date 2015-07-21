
#include <iostream>
#include <cstdlib>
#include <cmath>

#include <vector>
#include <list>

#include <gromacs/typedefs.h>
#include <gromacs/pbc.h>
#include <gromacs/xvgr.h>
#include <gromacs/vec.h>
#include <gromacs/princ.h>


#include "grid.h"

grid::grid(real spacing_, rvec minr_, rvec maxr_) :
    spacing(spacing_), //minr(minr_), maxr(maxr_), 
    nx((size_t)( (maxr_[XX] - minr_[XX])/spacing_)),
    ny((size_t)( (maxr_[YY] - minr_[YY])/spacing_)),
    nz((size_t)( (maxr_[ZZ] - minr_[ZZ])/spacing_)),
    grid_N(nx*ny*nz),
    grid_sum(nx*ny*nz),
    grid_inv_sum(nx*ny*nz)
{
    //copy_rvec(minr_, minr);
    //copy_rvec(maxr_, maxr);
    //zero_r[XX] = 

    // the grid indices of the zero coordinates
    int zx=int(0.5-minr_[XX]/spacing);
    if (zx < 1)
        zx=1;
    int zy=int(0.5-minr_[YY]/spacing);
    if (zy < 1)
        zy=1;
    int zz=int(0.5-minr_[ZZ]/spacing);
    if (zz < 1)
        zz=1;

    minr[XX] = (-zx+0.5)*spacing;
    minr[YY] = (-zy+0.5)*spacing;
    minr[ZZ] = (-zz+0.5)*spacing;
    maxr[XX] = minr[XX] + nx*spacing;
    maxr[YY] = minr[YY] + ny*spacing;
    maxr[ZZ] = minr[ZZ] + nz*spacing;

    printf("Allocating grid:   (%g -- %g, %g -- %g, %g -- %g)\n",
           minr_[XX], maxr_[XX],
           minr_[YY], maxr_[YY],
           minr_[ZZ], maxr_[ZZ]);
    printf("Allowed positions: (%g -- %g, %g -- %g, %g -- %g)\n",
           minr[XX], maxr[XX],
           minr[YY], maxr[YY],
           minr[ZZ], maxr[ZZ]);
    printf("    Grid size:       (%lu, %lu, %lu)\n", nx, ny, nz);
    //printf("    Grid zero point: (%ld, %ld, %ld)\n", zx, zy, zz);
    printf("    %gM grid points\n", (nx*ny*nz)/(1024.*1024.));
    printf("\n");


    for(size_t i=0; i<grid_N.size(); ++i)
    {
        grid_N[i] = 0;
        grid_sum[i] = 0.;
    }

    sample_count=0;
}

void grid::avg(real &arithmetic, real &harmonic, size_t &Na, size_t &Nh)
{
    real asum=0;
    real hsum=0;
    Na=0;
    Nh=0;

    for(size_t i=0;i<grid_N.size();i++)
    {
        int N=grid_N[i];
        real gsum=grid_sum[i];

        hsum += N/gsum;
        Nh++;
        if (N>0)
        {
            asum += gsum/N; 
            Na++;
        }
    }

    arithmetic=asum/Na;
    harmonic=Nh/hsum;
}


void write_grids(const char *filename, const grid *ref_grid,
                 const grid *pers_grid, 
                 const grid *exch_grid)
{
    if ( (exch_grid->get_nx() != pers_grid->get_nx()) ||
         (exch_grid->get_ny() != pers_grid->get_ny()) ||
         (exch_grid->get_nz() != pers_grid->get_nz()) )
    {
        fprintf(stderr, "ERROR ERROR grids unequal in write_grids!!!!");
        fprintf(stdout, "ERROR ERROR grids unequal in write_grids!!!!");
        exit(1);
    }

    size_t nx=exch_grid->get_nx();
    size_t ny=exch_grid->get_ny();
    size_t nz=exch_grid->get_nz();
    const rvec &minr(exch_grid->get_minr());
    double spacing=exch_grid->get_spacing();
    double zx=(0 - minr[XX])/spacing;
    double zy=(0 - minr[YY])/spacing;
    double zz=(0 - minr[ZZ])/spacing;
    size_t sample_count=ref_grid->get_sample_count();
    size_t sample_size=ref_grid->get_sample_size();

    //double ref_norm=((double)sample_count * (double)sample_size);
    double ref_norm=((double)sample_count)*spacing*spacing*spacing;

    FILE *outf=fopen(filename, "w");
    
    fprintf(outf, "# Nx Ny Nz spacing zero_x zero_y zero_z\n");
    fprintf(outf, "%lu %lu %lu %g %g %g %g\n", nx, ny, nz, 
            spacing, zx, zy, zz);
    fprintf(outf, "# ix iy iz N_ref N_pers sum<t_pers> sum<1/t_pers> N_exch sum<t_exch> sum<1/t_exch>\n");

    for(size_t iz=0; iz<nz; ++iz)
    {
        for(size_t iy=0; iy<ny; ++iy)
        {
            for(size_t ix=0; ix<nx; ++ix)
            {
                size_t index=pers_grid->get_index(ix, iy, iz);

                unsigned int N_ref =    ref_grid->get_N(index);
                double N_ref_norm =     N_ref/ref_norm;

                unsigned int N_pers =   pers_grid->get_N(index);
                double sumt_pers =      pers_grid->get_sum(index);
                double suminv_pers =    pers_grid->get_inv_sum(index);

                unsigned int N_exch =   exch_grid->get_N(index);
                double sumt_exch =      exch_grid->get_sum(index);
                double suminv_exch =    exch_grid->get_inv_sum(index);


                fprintf(outf, "%lu %lu %lu %g %u %g %g %u %g %g\n", 
                        ix, iy, iz, N_ref_norm, 
                        N_pers, sumt_pers, suminv_pers, 
                        N_exch, sumt_exch, suminv_exch);
            }
        }
    }
    fclose(outf);
}

