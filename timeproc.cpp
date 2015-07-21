
#include <cstdio>
#include <cmath>
#include <list>
#include <vector>
#include <string>
#include <deque>

#define GMX_THREADS


#if 0
#include <gromacs/typedefs.h>
#if !defined(GMX_API_VERSION) || (GMX_API_VERSION < 40600)
#error 4.5 
#elif (GMX_API_VERSION < 40700)
#error 4.6 incompatible
#endif
#endif

#if 0
#include <gromacs/version.h>
#if (GMX_API_VERSION < 50000)
#error pre-5.0
#else
#error post-5.0
#endif
#endif

#include <gromacs/typedefs.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/statutil.h>
#include <gromacs/copyrite.h>
#include <gromacs/macros.h>
#include <gromacs/sysstuff.h>
#include <gromacs/string2.h>
#include <gromacs/smalloc.h>
#include <gromacs/tpxio.h>
#include <gromacs/index.h>
#include <gromacs/confio.h>
#include <gromacs/pbc.h>
#include <gromacs/vec.h>
#include <gromacs/mtop_util.h>

#include "hist.h"
#include "event.h"
#include "grid.h"
#include "proccom.h"



static void proc_traj(const char *conffn, 
                      const char *tprfn, 
                      const char *ndxfn,
                      const char *logfn,
                      const char *pers_exch_histfn,
                      const char *zhistfn,
                      const char *rhist,
                      const char *drhist,
                      const char *distpdb,
                      const char *gridfn,
                      real jumpdist,
                      real weightsigma,
                      int closest_n,
                      real max_group_dist,
                      real max_group_atom_dist,
                      bool correct_group_rotation,
                      bool fit,
                      real grid_spacing,
                      gmx_bool z_absolute,
                      real z_min,
                      real z_max,
                      real z_spacing,
                      real r_spacing,
                      real dr_spacing,
                      real first_event_t,
                      real last_event_t,
                      int nstevent,
                      event_direction evdir,
                      output_env_t oenv)
{
    int natoms;
    //matrix box;
    real t;
    rvec *x;
    t_trxstatus *status;
    int framenr=0;
    int ePBC=-1;
    int *isize;
    char **grpname;
    atom_id **index;
    int i,j;
    gmx_mtop_t mtop;
    t_tpxheader tpxh;
    int version, generation;
    t_pbc pbc;
    int ngroups=1;
    t_trxframe fr;
    bool term_msg_sent=false;

    // read the topology
    read_tpxheader(tprfn, &tpxh, TRUE, &version , &generation);
    ePBC = read_tpx(tprfn,NULL,NULL,&natoms,NULL,NULL,NULL,&mtop);

    // read the first frame of the trajectory
    //natoms = read_first_x(oenv, &status, conffn, &t, &x, box);
    int ret=read_first_frame(oenv, &status, conffn, &fr, TRX_NEED_X);
    if (!ret)
        gmx_fatal(FARGS, "Couldn't open %s\n",conffn);



    gmx_mtop_finalize(&mtop);
    t_atoms atoms=gmx_mtop_global_atoms(&mtop);

    printf("Topology read. n=%d\n", fr.natoms);    

    printf("Select the group to track");
    if (weightsigma > 0 || closest_n > 0)
    {
        printf(" and the strain reference group");
        ngroups=2;
    }
    else if (max_group_atom_dist > 0)
    {
        printf(" and the group to base max. atom distance on");
        ngroups=2;
    }
    else if (max_group_dist > 0)
    {
        printf(" and the group to base max. distance on");
        ngroups=2;
    }
    else if ( (zhistfn) && (!z_absolute) )
    {
        printf(" and the z-coordinate reference group");
        ngroups=2;
    }
    else if (grid_spacing > 0)
    {
        printf(" and the reference group");
        ngroups=2;
    }
    printf(":\n");
    snew(index,ngroups);
    snew(isize,ngroups);
    snew(grpname,ngroups);
    //rd_index(ndxfn, ngroups, isize, index, grpname);
    get_index(&atoms, ndxfn, ngroups, isize, index, grpname);
#if 1
    if (fr.natoms != tpxh.natoms)
    {
        gmx_fatal(FARGS,"Number of atoms in the run input file (%d) is different from the trajectory (%d)", tpxh.natoms, fr.natoms );
    }

    for(i=0; (i<ngroups); i++)
    {
        for(j=0; (j<isize[i]); j++)      
        {
            if (index[i][j] > fr.natoms)
            gmx_fatal(FARGS,"An atom number in group %s is larger than the number of atoms in the trajectory");
        }
    }
#endif

    bool do_rcalc=(rhist != NULL) && (r_spacing>0);
    bool do_z_hist=(zhistfn) && (z_spacing>0);
    bool do_z_rho_hist=(zhistfn ) && (z_spacing>0);
    //bool correct_group_rotation= (max_group_dist > 0);
    //printf("do_rcalc=%d\n", do_rcalc);
    comsystem cs(fr.natoms, jumpdist, weightsigma, closest_n, 
                 max_group_dist,
                 max_group_atom_dist,
                 correct_group_rotation, fit,
                 z_absolute!=FALSE, z_min, z_max, z_spacing, r_spacing, 
                 dr_spacing, evdir,
                 do_rcalc, do_z_hist, do_z_rho_hist, distpdb != NULL,
                 grid_spacing);

    fprintf(stderr, "Categorizing molecules\n");
    for(j=0;j<ngroups;j++)
    {
        fprintf(stderr, "Categorizing molecules: group %d\n", j);
        for(i=0;i<isize[j];i++)
        {
            int molb,molnr, atnr_mol;
            int id=index[j][i];

            int resnr; 
            char *atomname; 
            char *resname;
 

            gmx_mtop_atomnr_to_molblock_ind(&mtop, index[j][i], 
                                            &molb, &molnr, &atnr_mol);
            gmx_mtop_atominfo_global(&mtop, index[j][i], &atomname, &resnr, 
                                     &resname);

#if 0
            printf("Atom %d [%d][%d] in molblock %d, molnr %d, mass %g\n",
                   id, j, i, molb, molnr, atoms.atom[index[j][i]].m);
            if (i%1000==0)
            {
                fprintf(stderr, "group %d, N=%d\n", j, i);
            }
#endif
            cs.init_atom(id, molb, molnr, j, resname, atomname);
        }
    }

    if (fr.natoms == 0)
        gmx_fatal(FARGS,"No atoms in trajectory!");

    do
    {
        bool make_exch_events = (last_event_t<0) || (fr.time<last_event_t);
        bool make_pers_events = ((framenr%nstevent)==0) && make_exch_events;
        bool last_frame=false;

        if ( !make_exch_events )
        {
            if (!term_msg_sent)
            {
                printf("\nLast event start time reached. Handling remaining %d events.\n", 
                       cs.events_left());
                term_msg_sent=true;
                last_frame=true;
            }
            if (cs.events_left()==0)
            {
                printf("\nStopping: all events handled.\n");
                break;
            }
        }

        if (fr.time >= first_event_t)
        {

            /*proc_frame(&mtop, &pbc, ePBC, fr.natoms, ngroups, fr.x, fr.box, 
                       fr.time, isize, index, &atoms, cs, make_pers_events,
                       make_exch_events, framenr);*/

            // process the frame
            cs.reset_coms();
            set_pbc(&pbc, ePBC, fr.box);
            // first calculate the COMs 
            for(int j=0;j<ngroups;j++)
            {
                for(int i=0;i<isize[j];i++)
                {
                    int id=index[j][i]; 
                    cs.add_atom_coord(id, fr.x[id], atoms.atom[id].m, &pbc);
                }
            }
            // for the first frame, set the start coords for pbc jumps
            if (framenr == 0)
            {
                for(int j=0;j<ngroups;j++)
                {
                    for(int i=0;i<isize[j];i++)
                    {
                        int id=index[j][i]; 

                        cs.set_start_coord(id, fr.x[id]);
                    }
                }
            }

            bool writing_out_coords=((framenr % 1000) == 0) && 
                      (correct_group_rotation||fit);
            

            /* now the coordinates are used for COM calculations. We can 
               safely change them for principal axis orientation, 
               discretization, etc. */
            cs.handle_atoms(&pbc, fr.x, fr.natoms, ngroups, isize[ngroups-1], 
                            index[ngroups-1],
                            atoms.atom, make_pers_events, make_exch_events, 
                            fr.time, writing_out_coords);
            if (writing_out_coords)
            {
                char filename[256];
                sprintf(filename, "out%06d.gro", framenr);
                write_sto_conf(filename, "bla", &atoms, fr.x, NULL, ePBC, 
                               fr.box);
            }
            // handle the COMs. this is where most of the work takes place.
            cs.handle_mols(&pbc, make_pers_events, make_exch_events, fr.time);
            
            framenr++;
        }

        //if (last_frame && distpdb)
        if (distpdb)
        {
            char filename[512];
            sprintf(filename, "%05d-%s", framenr, distpdb);
            cs.write_pdb(filename, &pbc);
        }
    } while(read_next_frame(oenv, status, &fr));

  

    close_trj(status);

    printf("       Atoms: %d\n", natoms);
    printf("Frame number: %d\n", framenr);

    cs.process_events(stdout, 100, pers_exch_histfn, 
                      zhistfn, 
                      gridfn,
                      rhist, 
                      drhist, 
                      oenv);
    if (logfn)
    {
        FILE *log=fopen(logfn, "w");
        cs.process_events(log, 100, 
                          NULL, 
                          NULL, 
                          NULL, 
                          NULL, 
                          NULL, oenv);
        fclose(log);
    }
}

int main(int argc, char *argv[])
{
    const char *desc[] = {
        "Persistence vs exchange time calculation tool: uses a cutoff ",
        "distance between events, set with -j.\n",
        "\n",
        "Displacements can be made relative to a reference set of molecules",
        "with either a 'sigma weight' -w that, if postive, sets the decay",
        "length of a weight function with which to calculate the",
        "displacement of a particle relative to its neighboring particles.",
        "Alternatively, a fixed number of neighbouring reference particles",
        "can be set with -cn'. If neither is set, absolute displacements are",
        "used.\n",
        "\n",
        "If -zt is set, the persistence & exchange times are ",
        "plotted w.r.t. the z coordinate. If -zn is set, the number of events",
        "is plotted w.r.t. the z coordinate.\n",
        "\n",
        "If -zrho is set, the density is plotted for both the tracked atoms'",
        "COMs and the reference atoms' COMs.\n",
        "\n",
        "With -jumph, the actual measured jump distances (at exchange events)",
        "are binned\n" 
    };

    // the jumping distance. Based on the mean distance between water oxygens 
    real jumpdist=0.28;
    // this is just to give a number of anchors close to 10
    real weightsigma=-1;
    // The nubmer of particles to have in the reference group
    int closest_n=0;

    // directions to check
    const char *dirstr[]={NULL,"xyz", "xy", "z", NULL };
    event_direction evdir=evdir_xyz;
    // whether to use absolute coordinates a z values
    gmx_bool z_absolute=FALSE;
    // spacing of the z histogram (or no histogram if <=0)
    real z_spacing=0.1;
    // spacing of the r histogram 
    real r_spacing=0.01;
    // spacing of the dr histogram 
    real dr_spacing=0.001;
    // steps between persistence events
    int nstevent=50;
    // last time to allow starts of events
    real first_event_t=0;
    // last time to allow starts of events
    real last_event_t=-1;
    // minimum and maximum z-value to look at
    real z_min=-1e10;
    real z_max=1e10;

    // grid spacing
    real grid_spacing=-1;
    // maximum distance from index group to measure
    real max_group_dist=-1;
    // maximum distance from any atom in the index group to measure
    real max_group_atom_dist=-1;
    // whether to correct for the index group rotation
    gmx_bool correct_rotation=FALSE;
    // whether to correct for the index group shifts with lsq fitting
    gmx_bool fit=FALSE;

    t_pargs pa[] =
    {
        {"-j", FALSE, etREAL, {&jumpdist}, "Jump distance" },
        {"-w", FALSE, etREAL, {&weightsigma}, 
            "Weight sigma (negative = no weights)" },
        {"-cn", FALSE, etINT, {&closest_n}, 
            "Closest N in ref group (0 = no ref group or use weightsigma)" },
        {"-ddir", FALSE, etENUM, {dirstr}, 
            "Direction (xyz/xy/z) to check distances in"},
        {"-ne", FALSE, etINT, {&nstevent}, 
            "Steps between persistence time samples"},
        {"-zabs", FALSE, etBOOL, {&z_absolute}, 
            "Whether to use absolute coordinates a z values"},
        {"-zmin", FALSE, etREAL, {&z_min}, 
            "min. z value for event measurement. (no minimum if min>max)"},
        {"-zmax", FALSE, etREAL, {&z_max}, 
            "max. z value for event measurement. (no maximum if min>max)"},
        {"-group-dist", FALSE, etREAL, {&max_group_dist}, 
            "max distance to the COM of a given index group for event measurement (<0 for no max dist)"},
        {"-group-atom-dist", FALSE, etREAL, {&max_group_atom_dist}, 
            "max distance to the closest atom of a given index group for event measurement (<0 for no max dist)"},
        {"-correct-group-rot", FALSE, etBOOL, {&correct_rotation},
            "whether to correct for reference group rotation"},
        {"-fit", FALSE, etBOOL, {&fit},
            "whether to correct for reference group shifts by fitting with lsq"},
        {"-grid", FALSE, etREAL, {&grid_spacing}, 
            "grid spacing. Use >0 to assign events to grid points"},
        {"-zd", FALSE, etREAL, {&z_spacing}, 
            "z spacing for z histograms"},
        {"-rd", FALSE, etREAL, {&r_spacing}, 
            "r spacing for r histograms."},
        {"-rdr", FALSE, etREAL, {&dr_spacing}, 
            "dr spacing for jumpdist histograms."},
        {"-fet", FALSE, etREAL, {&first_event_t}, "first event start time."},
        {"-let", FALSE, etREAL, {&last_event_t}, "last event start time."}
    };
    t_filenm fnm[] = {
        { efTRX, "-f", "traj", ffREAD },
        { efTPS, "-s", NULL, ffREAD },
        { efNDX, "-n", NULL, ffOPTRD },
        { efLOG, "-l", NULL, ffOPTRD },
        { efXVG, "-pxh", "pers_exch", ffWRITE },
        { efXVG, "-zh", "zhist", ffOPTWR },
        { efXVG, "-rh", "rh", ffOPTWR },
        { efXVG, "-jumph", "jumps", ffOPTWR },
        { efPDB, "-dpdb", "dists", ffOPTWR },
        { efDAT, "-g", "grid", ffWRITE }
    };
#define NFILE asize(fnm)
    output_env_t oenv;
    int npargs;

    npargs=asize(pa);
    /*ppa = add_acf_pargs(&npargs, pa);*/

    parse_common_args(&argc,argv,PCA_CAN_TIME,NFILE,fnm,npargs,
                      pa,asize(desc),desc,0, NULL,&oenv);

    evdir=(event_direction)(nenum(dirstr)-1);
    //printf("direction=%s, %d\n", dirstr[0], (int)evdir);

    if (fit && correct_rotation)
    {
        gmx_fatal(FARGS, 
                  "correct_group_rotation cannot be used together with fit");
    }

    proc_traj(opt2fn("-f", NFILE, fnm), 
              opt2fn("-s", NFILE, fnm),
              ftp2fn_null(efNDX, NFILE, fnm),
              ftp2fn_null(efLOG, NFILE, fnm),
              opt2fn("-pxh", NFILE, fnm),
              opt2fn_null("-zh", NFILE, fnm),
              opt2fn_null("-rh", NFILE, fnm),
              opt2fn_null("-jumph", NFILE, fnm),
              opt2fn_null("-dpdb", NFILE, fnm),
              opt2fn("-g", NFILE, fnm),
              jumpdist, 
              weightsigma,
              closest_n,
              max_group_dist,
              max_group_atom_dist,
              correct_rotation,
              fit,
              grid_spacing,
              z_absolute,
              z_min,
              z_max,
              z_spacing,
              r_spacing,
              dr_spacing,
              first_event_t,
              last_event_t,
              nstevent,
              evdir,
              oenv);

    return 0;
}


