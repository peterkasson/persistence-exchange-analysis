

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <gromacs/typedefs.h>

#include "hist.h"

void hist_base::write(const std::string &filename, bool normalize, 
                      bool histtype, bool error, bool write_zeroes)
{
    FILE *out=fopen(filename.c_str(), "w");
    if (!out)
    {
        fprintf(stderr,"ERROR: opening file %s\n", filename.c_str());
        exit(1);
        //throw std::iostream::failure("Error opening file");
    }
    writef(out, normalize, histtype, write_zeroes);
    fclose(out);
}


void hist::writef(FILE *out, bool normalize, bool histtype, bool error, 
                  bool write_zeroes)
{
    double normfact=1.;

    if (normalize)
        normfact=1./(double(sum)*bin_width);

    for(unsigned int i=0;i<N;i++)
    {
        double x1=bin_width*i + minx;
        double x2=bin_width*(i+1) + minx;

        double y=normfact*double(count[i]);

        if (histtype)
        {
            if (write_zeroes || count[i]!=0)
            {
                if (error)
                    fprintf(out, "%g %g 0\n%g %g 0\n", x1, y, x2, y);
                else
                    fprintf(out, "%g %g\n%g %g\n", x1, y, x2, y);
            }
        }
        else
        {
            if (write_zeroes || count[i]!=0)
            {
                if (error)
                    fprintf(out, "%g %g 0\n", (x1+x2)/2., y);
                else
                    fprintf(out, "%g %g\n", (x1+x2)/2., y);
            }
        }
    }
}

void weighted_hist::writef(FILE *out, bool normalize, bool histtype, 
                           bool error, bool write_zeroes)
{
    double normfact=1.;

    for(unsigned int i=0;i<N;i++)
    {
        if (count[i]>0)
        {
            bool doprint=true;
            bool have_err=false;
            double x1=bin_width*i + minx;
            double x2=bin_width*(i+1) + minx;

            double w=sum_weight[i]/count[i];

            double y=w;
            double yerr=0;
            if (nblock>1)
            {
                if (blktot_count[i] > 1)
                {
                    /* we calculate the average again because the block's
                       average may be weighted differently */
                    double wsq = blktot_sq_weight[i]/blktot_count[i];
                    double wn = blktot_weight[i]/blktot_count[i];
                    yerr = sqrt( wsq - wn*wn )/(blktot_count[i]-1);
                    have_err=true;
                    if (!write_zeroes && count[i]==0)
                        doprint=false;
                }
                else if (blktot_count[i]==-1)
                {
                    /* a hack: we use the sq_weight directly */
                    yerr = blktot_sq_weight[i]; 
                    if (!write_zeroes && count[i]==0)
                        doprint=false;
                    have_err=true;
                }
                else
                {
                    doprint=false;
                }
            }
            else
            {
                if (blktot_count[i]==-1)
                {
                    /* a hack: we use the sq_weight directly */
                    yerr = blktot_sq_weight[i]; 
                    have_err=true;
                }
            }


            if (histtype)
            {
                if (have_err) 
                {
                    if (doprint)
                    {
                        if (error)
                            fprintf(out, "%g %g %g\n%g %g %g\n", x1, y, yerr, 
                                    x2, y, yerr);
                        else
                            fprintf(out, "%g %g\n%g %g\n", x1, y, x2, y);

                    }
                }
                else
                {
                    if (error)
                        fprintf(out, "%g %g 0\n%g %g 0\n", x1, y, x2, y);
                    else
                        fprintf(out, "%g %g\n%g %g\n", x1, y, x2, y);

                }
            }
            else
            {
                if (have_err) 
                {
                    if (doprint)
                    {
                        if (error)
                            fprintf(out, "%g %g %g\n", (x1+x2)/2., y, yerr);
                        else
                            fprintf(out, "%g %g\n", (x1+x2)/2., y);

                    }
                }
                else
                {
                    if (error)
                        fprintf(out, "%g %g 0\n", (x1+x2)/2., y);
                    else
                        fprintf(out, "%g %g\n", (x1+x2)/2., y);
                }
            }
        }
    }
}


void weighted_hist::finish_block()
{
    for(unsigned int i=0;i<N;i++)
    {
        if (blksum_count[i] > 0)
        {
            double bs=blksum_weight[i]/blksum_count[i];
            blktot_weight[i] += bs;
            blktot_sq_weight[i] += bs*bs;
            blktot_count[i]++;
        }
        blksum_weight[i]=0;
        blksum_count[i]=0;
    }
    nblock++;
}


