#!/usr/bin/env python

import numpy
import scipy.integrate
import math
import sys


#constants:
# Avogadro's number
Na=6.0221e23
# Absolute viscosity of water at 298K:
eta0=1.002e-3
# Temperature
T=293
# boltzmann's constant
k=1.38065e-23

class grid:
    def __init__(self, filename):
        inf=open(filename, 'r')
        # read comment line
        inf.readline()
        sizelinesp=inf.readline().split()
        nx=int(sizelinesp[0])
        ny=int(sizelinesp[1])
        nz=int(sizelinesp[2])
        spacing=float(sizelinesp[3])
        zx=float(sizelinesp[4])
        zy=float(sizelinesp[5])
        zz=float(sizelinesp[6])
        print "nx=%d, ny=%d, nz=%d"%(nx, ny, nz)
        print "zx=%g, zy=%g, zz=%g"%(zx, zy, zz)
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.zx=zx
        self.zy=zy
        self.zz=zz
        self.spacing=spacing
        # read comment line
        inf.readline()
        # now allocate the arrays
        #self.N_ref = numpy.zeros((nx,ny,nz), dtype=numpy.float)

        #self.N_pers = numpy.zeros((nx,ny,nz), dtype=numpy.int)
        #self.sum_t_pers = numpy.zeros((nx,ny,nz), dtype=numpy.float)

        #self.N_exch = numpy.zeros((nx,ny,nz), dtype=numpy.int)
        #self.sum_t_exch = numpy.zeros((nx,ny,nz), dtype=numpy.float)

        #self.use = numpy.zeros((nx,ny,nz), dtype=numpy.bool)
        #self.grid=numpy.array
        
        # The points as functions of r
        rmax=self.spacing*min(nx,ny,nz)/2
        print "rmax=%g"%rmax
        self.r_bin_size = self.spacing
        self.Nr = int(rmax/self.r_bin_size)
        print "Nr=%d"%self.Nr
        self.pers_av_r = numpy.zeros( (self.Nr), dtype=numpy.float)
        self.pers_av_N = numpy.zeros( (self.Nr), dtype=numpy.int)
        for line in inf:
            spl=line.split()
            ix=int(spl[0])
            iy=int(spl[1])
            iz=int(spl[2])
            N_ref=float(spl[3])
            N_pers=int(spl[4])
            sum_t_pers=float(spl[5])
            sum_t_pers_inv=float(spl[6])
            N_exch=int(spl[7])
            sum_t_exch=float(spl[8])
            sum_t_exch_inv=float(spl[9])
            #self.N_ref[ix,iy,iz] = N_ref
            #self.N_pers[ix,iy,iz] = N_pers
            #self.sum_t_pers[ix,iy,iz] = sum_t_pers
            #self.N_exch[ix,iy,iz] = N_exch
            #self.sum_t_exch[ix,iy,iz] = sum_t_exch
            #print "%d  %d  %d  -> %d %s"%(ix,iy,iz, N_ref, spl)
            # calculate the distance to the zero coordinate
            rx=(ix+0.5) - zx
            ry=(iy+0.5) - zy
            rz=(iz+0.5) - zz
            r = math.sqrt(rx*rx + ry*ry + rz*rz)*self.spacing
            #print r
            rbin = int(r/self.r_bin_size)
            if rbin < self.Nr:
                #if N_pers > 0:
                #    self.pers_av_r [rbin] += N_pers/sum_t_pers
                if N_pers > 0:
                    self.pers_av_r [rbin] += sum_t_pers_inv/N_pers
                self.pers_av_N [rbin] += 1
                #self.mark_ref_close(ix,iy,iz, grid_dist_check)

    #def mark_ref_close(self, ix, iy, iz, grid_dist_check):
    #    #print "%d  %d  %d"%(ix,iy,iz)
    #    zmin=max(iz-grid_dist_check, 0)
    #    zmax=min(iz+grid_dist_check+1, self.nz)
    #    ymin=max(iy-grid_dist_check, 0)
    #    ymax=min(iy+grid_dist_check+1, self.ny)
    #    xmin=max(ix-grid_dist_check, 0)
    #    xmax=min(ix+grid_dist_check+1, self.nx)
    #    for diz in range(zmin, zmax):
    #        for diy in range(ymin, ymax):
    #            for dix in range(xmin, xmax):
    #                self.use[dix, diy, diz] = True
    #                #if self.N_ref[dix, diy, diz] >= lower_ref_threshold:
    #                #    return True
    #    #return False
 
    #def have_ref_close(self, ix, iy, iz, grid_dist_check, lower_ref_threshold):
    #    zmin=max(iz-grid_dist_check, 0)
    #    zmax=min(iz+grid_dist_check+1, self.nz)
    #    ymin=max(iy-grid_dist_check, 0)
    #    ymax=min(iy+grid_dist_check+1, self.ny)
    #    xmin=max(ix-grid_dist_check, 0)
    #    xmax=min(ix+grid_dist_check+1, self.nx)
    #    for diz in range(zmin, zmax):
    #        for diy in range(ymin, ymax):
    #            for dix in range(xmin, xmax):
    #                if self.N_ref[dix, diy, diz] >= lower_ref_threshold:
    #                    return True
    #    return False
                                

    def write_avg(self, filename, tp0):
        outf=open(filename, 'w')
        for i in range(self.Nr):
            r=(i+0.5)*self.r_bin_size
            outf.write("%g %g %d\n"%(r, tp0*self.pers_av_r[i]/self.pers_av_N[i],
                                     self.pers_av_N[i]))
        outf.close()


    def write_int(self, filename, tp0, lower_n_threshold):
        outf=open(filename, 'w')
        x=[]
        y=[]
        for i in range(self.Nr):
            if self.pers_av_r[i] > 0 and self.pers_av_N[i]>=lower_n_threshold:
                r=(i+0.5)*self.r_bin_size
                ri=(tp0*self.pers_av_r[i]/self.pers_av_N[i])/(r*r*r*r)
                x.append(r)
                y.append(ri)
                outf.write("%g %g %d\n"%(r, ri, self.pers_av_N[i]))
        intg=scipy.integrate.simps(y, x=x)
        print("integral [%g-%g] = %g"%(x[0], x[-1], intg))
        outf.close()

    def calc_avgs(self, grid_spacing, lower_ref_threshold, lower_avg_threshold,
                  grid_dist_check,tp0):
        pers_sum=0.
        N_tot=0
        N_cosol=0
        N_sol=0
        N_pers_max=0
        for iz in range(self.nz):
            for iy in range(self.ny):
                for ix in range(self.nx):
                    # first check whether we need to include it into our count
                    N_pers=self.N_pers[ix, iy, iz]
                    if self.use[ix,iy,iz]:
                        # we always count
                        N_tot += 1
                        if N_pers >= lower_avg_threshold:
                            pers_sum += N_pers/self.sum_t_pers[ix, iy, iz]
                            N_cosol += 1
                    if N_pers > N_pers_max:
                        N_pers_max=N_pers
        pers_avg_tot = N_tot/pers_sum
        pers_avg_cosol = N_cosol/pers_sum
        cell_size=(grid_spacing*grid_spacing*grid_spacing)
        #N_sol = N_tot - N_cosol
        N_sol = self.Nsol
        vol_tot = N_tot * cell_size
        vol_sol = N_sol * cell_size
        vol_cosol = N_cosol * cell_size
        Drel=1. + (tp0*pers_sum - float(N_cosol))/float(N_tot)
        # the radius from the volume
        a3=(3*vol_sol*1e-27 / (4*math.pi))
        D_SE=k*T/(8*math.pi*eta0*a3)/1e6
        D_abs=D_SE*Drel
        print("Total average persistence:       %g"%pers_avg_tot)
        print("Co-solvent average persistence:  %g"%pers_avg_cosol)
        print("Total rel. average persistence:       %g"%(tp0/pers_avg_tot))
        print("Co-solvent rel. average persistence:  %g"%(tp0/pers_avg_cosol))
        #print("<t^0>*<1/t>:  %g"%(tp0*pers_sum/N_cosol))
        print("N_tot    %d"%N_tot)
        print("N_sol    %d"%N_sol)
        print("N_co-sol %d"%N_cosol)
        print("V_tot    %g"%vol_tot)
        print("V_sol    %g"%vol_sol)
        print("V_co-sol %g"%vol_cosol)
        print("N_co-sol/N_tot = %g"%(float(N_cosol)/float(N_tot)))
        print("N_pers_max %d"%N_pers_max)
        print("D/D^SE = %g"%Drel)
        print("D^SE = %g"%D_SE)
        print("D = %g"%D_abs)
        

if len(sys.argv) != 4:
    print "usage: "
    print "gridavg grid.dat lower-n-threshold t_p0"
    print
    print "calculates the average exch and pers times over a grid,"
    print "for grid points that have a grid point nearby (within "
    print "grid-dist-check that has at least lower-ref-threshold reference "
    print "counts."
    print
    print "Counts points with fewer than lower-avg-threshold as infinite"
    print "t_p0 = the average fluid persistence time."
    sys.exit(1)

filename=sys.argv[1]
#grid_spacing=float(sys.argv[2])
lower_n_threshold=int(sys.argv[2])
#lower_avg_threshold=int(sys.argv[3])
#grid_dist_check=int(sys.argv[5])
tp0=float(sys.argv[3])

g=grid(filename)

g.write_avg('eta-1.dat', tp0)
g.write_int('eta-integrate.dat', tp0, lower_n_threshold)

#g.calc_avgs(grid_spacing, lower_ref_threshold, lower_avg_threshold, 
#            grid_dist_check,tp0)

