


class hist_base
{
public:
    hist_base(size_t N_, double minx_, double maxx_) 
              : N(N_), minx(minx_), maxx(maxx_)
    {
        bin_width=(maxx-minx)/N;
    }
    size_t get_bin(double x) const
    {
        if (x>=minx && x<maxx)
            return (int)((x-minx)/bin_width);
        else
            return N+1;
    }

    double get_x(size_t bin) const
    {
        return minx + bin_width*bin;
    }

    size_t get_N(void) const { return N; }

    // write histogram to file name
    virtual void write(const std::string &filename, bool normalize=false,
                       bool histtype=true, bool error=false, 
                       bool write_zeroes=true);
    // write histogram to open file
    virtual void writef(FILE *out, bool normalize=false, bool histtype=true,
                        bool error=false, bool write_zeroes=true)=0;
protected:
    size_t N; // number of bins
    double minx, maxx; // min and max x values
    double bin_width; // the bin width
}; 


class hist : public hist_base
{
public:
    hist(size_t N_, double minx_, double maxx_) 
              : hist_base(N_, minx_, maxx_), count(N_, 0), sum(0)
    {
    }

    void add(double x)
    {
        size_t bin=get_bin(x);
        if (bin<N)
        {
            count[bin]++;
        }
        sum++;
    }

    // write histogram to open file
    virtual void writef(FILE *out, bool normalize=false, bool histtype=true,
                        bool error=false, bool write_zeroes=true);

    // reset the histogram
    void reset()
    {
        for(unsigned int i=0;i<N;i++)
            count[i]=0;
        sum=0;
    }

    unsigned long get(size_t bin) const
    {
        return count[bin];
    }

protected:
    std::vector<unsigned long> count;
    unsigned long sum;
};    


class weighted_hist : public hist_base
{
public:
    weighted_hist(size_t N_, double minx_, double maxx_) 
              : hist_base(N_, minx_, maxx_), 
                count(N_, 0), sum_weight(N_, 0.), 
                blksum_weight(N_, 0),
                blksum_count(N_, 0),
                blktot_sq_weight(N_, 0),
                blktot_weight(N_, 0),
                blktot_count(N_, 0),
                sum(0), sumw(0.), sumwsq(0.), nblock(0)
    {
    }

    // add to a specific bin
    void add_bin(size_t bin, double w)
    {
        if (bin<N)
        {
            sum_weight[bin] += w;
            count[bin]++;

            blksum_weight[bin] += w;
            blksum_count[bin]++;
        }
        sum++;
        sumw+=w;
        sumwsq+=w*w;
    }

    void add(double x, double w)
    {
        size_t bin=get_bin(x);
        add_bin(bin, w);
    }

    void reset()
    {
        for(unsigned int i=0;i<N;i++)
        {
            count[i]=0;
            sum_weight[i]=0;
            blksum_weight[i]=0;
            blksum_count[i]=0;
            blktot_sq_weight[i]=0;
            blktot_weight[i]=0;
            blktot_count[i]=0;
        }
        sum=0;
        sumw=0; 
        sumwsq=0; 
        nblock=0;
    }


    double get(size_t bin) const
    {
        return sum_weight[bin]/count[bin];
    }
    double get_err(size_t bin) const
    {
        if (nblock > 1 && blktot_count[bin]>1)
        {
            /* we calculate the average again because the block's
               average may be weighted differently */
            double wsq = blktot_sq_weight[bin]/blktot_count[bin];
            double wn = blktot_weight[bin]/blktot_count[bin];
            double yerr = sqrt( wsq - wn*wn )/(blktot_count[bin]-1);
            return yerr;
        }
        else if (blktot_count[bin] < 0)
        {
            return blktot_sq_weight[bin];
        }
        else
            return 0;
    }
    bool have_err(size_t bin) const
    {
        return (nblock > 1 && blktot_count[bin]>1);
    }

    void finish_block();

    // write histogram to open file
    virtual void writef(FILE *out, bool normalize=false, bool histtype=true,
                        bool error=false, bool write_zeroes=true);


    // return the current block's count
    unsigned long get_cur_block_count(size_t bin) const
    {
        return blksum_count[bin];
    }

    // return the current block's average
    double get_cur_block_avg(size_t bin) const
    {
        return blksum_weight[bin]/blksum_count[bin];
    }

    // force a specific bin to have a specific value and error
    void set(size_t bin, double avg, double err)
    {
        count[bin]=1;
        sum_weight[bin]=avg;
        blktot_weight[bin]=avg;
        blktot_sq_weight[bin]=err;
        blktot_count[bin]=-1;
    }
protected:
    std::vector<unsigned long> count;
    std::vector<double> sum_weight;

    std::vector<double> blksum_weight; // the sum of the block so far
    std::vector<unsigned int> blksum_count; // the count for blocksum_weight

    std::vector<double> blktot_sq_weight; // the squared sum of all blocks
    std::vector<double> blktot_weight; // the sum of all blocks
    std::vector<int> blktot_count; // the count of all blocks

    unsigned long sum;
    double sumw;
    double sumwsq;
    unsigned int nblock; // number of blocks encountered
};


