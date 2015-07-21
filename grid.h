// 


class grid
{
public:
    grid(real spacing_, rvec minr_, rvec maxr_); 

    void add(real t, const rvec r)
    {
        if ( (r[XX] >= minr[XX]) &&
             (r[YY] >= minr[YY]) &&
             (r[ZZ] >= minr[ZZ]) )
        {
            size_t ix = (size_t) ((r[XX] - minr[XX])/spacing);
            size_t iy = (size_t) ((r[YY] - minr[YY])/spacing);
            size_t iz = (size_t) ((r[ZZ] - minr[ZZ])/spacing);

            // check whether the coordinates fit in the grid
            if ( (ix >= nx) || (iy >= ny) || (iz >= nz) )
            {
                // do nothing if they don't
                return;
            }
            size_t index=get_index(ix, iy, iz); 

            grid_N[index] ++ ;
            grid_sum[index] += t;
            grid_inv_sum[index] += 1./t;
        }
    }

    // get an index out of an int (x,y,z) index
    size_t get_index(size_t ix, size_t iy, size_t iz) const
    {
        return ix + nx*iy + nx*ny*iz;
    }

    /* calculate the arithmetic and harmonic grid mean

       arithmetic = the arithmetic mean over all non-zero bins
       harmonic = the harmonic mean over all bins */
    void avg(real &arithmetic, real &harmonic, size_t &Na, size_t &Nh); 


    size_t get_nx() const { return nx; }
    size_t get_ny() const { return ny; }
    size_t get_nz() const { return nz; }
    const rvec &get_minr() const { return minr; }
    const rvec &get_maxr() const { return maxr; }
    float get_spacing() const { return spacing; }

    unsigned int get_N(size_t index) const
    {
        return grid_N[index];
    }
    
    real get_sum(size_t index) const
    {
        return grid_sum[index];
    }

    real get_inv_sum(size_t index) const
    {
        return grid_inv_sum[index];
    }

    void add_sample() { sample_count++; }
    size_t get_sample_count() const { return sample_count; }

    void set_sample_size(size_t size) { sample_size = size; }
    size_t get_sample_size() const { return sample_size; }
protected:
    real spacing; // spacing

    //rvec minr, maxr; // min&max coords
    size_t nx, ny, nz;

    //long int zx, zy, zz; // the zero grid coordinates

    rvec minr, maxr; // min & max coordinates

    std::vector<unsigned int> grid_N; // exchange event number
    std::vector<real> grid_sum; // exchange time sum
    std::vector<real> grid_inv_sum; // exchange time sum

    size_t sample_count; // the sample count for the reference grid
    size_t sample_size; // the sample size for the reference grid
};



void write_grids(const char *filename, const grid *ref_grid,
                 const grid *pers_grid, const grid *exch_grid);

