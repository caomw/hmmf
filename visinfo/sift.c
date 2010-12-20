#include <vl/generic.h>
#include <vl/stringop.h>
#include <vl/pgm.h>
#include <vl/sift.h>
#include <vl/getopt_long.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define ERRF(msg, arg) {                                        \
    err = VL_ERR_BAD_ARG ;                                      \
    snprintf(err_msg, sizeof(err_msg), msg, arg) ;              \
    break ;                                                     \
}

#define ERR(msg) {                                              \
    err = VL_ERR_BAD_ARG ;                                      \
    snprintf(err_msg, sizeof(err_msg), msg) ;                   \
    break ;                                                     \
}


/* ----------------------------------------------------------------- */
/** @brief Keypoint ordering
 ** @internal
 **/
int
korder (void const* a, void const* b) {
    double x = ((double*) a) [2] - ((double*) b) [2] ;
    if (x < 0) return -1 ;
    if (x > 0) return +1 ;
    return 0 ;
}

/* ---------------------------------------------------------------- */
/** @brief SIFT driver entry point
 **/
int main(int argc, char **argv)
{
    /* algorithm parameters */
    double   edge_thresh  = 10 ;
    double   peak_thresh  = 10 ;
    double   magnif       = -1 ;
    int      O = 3, S = 3, omin = 0;

    vl_bool  err    = VL_ERR_OK ;
    char     err_msg [1024] ;
    int      exit_code          = 0 ;
    int      verbose            = 4 ;
    vl_bool  force_orientations = 0 ;


    while (argc--) 
    {
        char const      *name = *argv++ ;

        FILE            *in    = 0 ;
        vl_uint8        *data  = 0 ;
        vl_sift_pix     *fdata = 0 ;
        VlPgmImage       pim ;

        VlSiftFilt      *filt = 0 ;
        vl_size          q ;
        int              i ;
        vl_bool          first ;

        double           *ikeys = 0 ;
        int              nikeys = 0, ikeys_size = 0 ;

        /* open input file */
        in = fopen (name, "rb") ;
        if (!in) 
        {
            err = VL_ERR_IO ;
            snprintf(err_msg, sizeof(err_msg),
                    "Could not open '%s' for reading.", name) ;
            goto done ;
        }

        /* read PGM header */
        err = vl_pgm_extract_head (in, &pim) ;

        if (err) 
        {
            switch (vl_get_last_error()) {
                case  VL_ERR_PGM_IO :
                    snprintf(err_msg, sizeof(err_msg),
                            "Cannot read from '%s'.", name) ;
                    err = VL_ERR_IO ;
                    break ;

                case VL_ERR_PGM_INV_HEAD :
                    snprintf(err_msg, sizeof(err_msg),
                            "'%s' contains a malformed PGM header.", name) ;
                    err = VL_ERR_IO ;
                    goto done ;
            }
        }

        if (verbose)
            printf ("sift: image is %d by %d pixels\n",
                    pim. width,
                    pim. height) ;

        // allocate buffer
        data  = malloc(vl_pgm_get_npixels (&pim) *
                vl_pgm_get_bpp       (&pim) * sizeof (vl_uint8)   ) ;
        fdata = malloc(vl_pgm_get_npixels (&pim) *
                vl_pgm_get_bpp       (&pim) * sizeof (vl_sift_pix)) ;

        if (!data || !fdata) 
        {
            err = VL_ERR_ALLOC ;
            snprintf(err_msg, sizeof(err_msg),
                    "Could not allocate enough memory.") ;
            goto done ;
        }

        // read PGM body 
        err  = vl_pgm_extract_data (in, &pim, data) ;

        if (err) {
            snprintf(err_msg, sizeof(err_msg), "PGM body malformed.") ;
            err = VL_ERR_IO ;
            goto done ;
        }

        // convert data type 
        for (q = 0 ; q < (unsigned) (pim.width * pim.height) ; ++q) {
            fdata [q] = data [q] ;
        }

        // make filter
        filt = vl_sift_new (pim.width, pim.height, O, S, omin) ;

        if (edge_thresh >= 0) vl_sift_set_edge_thresh (filt, edge_thresh) ;
        if (peak_thresh >= 0) vl_sift_set_peak_thresh (filt, peak_thresh) ;
        if (magnif      >= 0) vl_sift_set_magnif      (filt, magnif) ;

        if (!filt) {
            snprintf (err_msg, sizeof(err_msg),
                    "Could not create SIFT filter.") ;
            err = VL_ERR_ALLOC ;
            goto done ;
        }

        if (verbose > 1) {
            printf ("sift: filter settings:\n") ;
            printf ("sift:   octaves      (O)     = %d\n",
                    vl_sift_get_noctaves     (filt)) ;
            printf ("sift:   levels       (S)     = %d\n",
                    vl_sift_get_nlevels      (filt)) ;
            printf ("sift:   first octave (o_min) = %d\n",
                    vl_sift_get_octave_first (filt)) ;
            printf ("sift:   edge thresh           = %g\n",
                    vl_sift_get_edge_thresh  (filt)) ;
            printf ("sift:   peak thresh           = %g\n",
                    vl_sift_get_peak_thresh  (filt)) ;
            printf ("sift:   magnif                = %g\n",
                    vl_sift_get_magnif       (filt)) ;
            printf ("sift: will source frames? %s\n",
                    ikeys ? "yes" : "no") ;
            printf ("sift: will force orientations? %s\n",
                    force_orientations ? "yes" : "no") ;
        }

        /* ...............................................................
         *                                             Process each octave
         * ............................................................ */
        i     = 0 ;
        first = 1 ;
        while (1) 
        {
            VlSiftKeypoint const *keys = 0 ;
            int                   nkeys ;

            /* calculate the GSS for the next octave .................... */
            if (first) 
            {
                first = 0 ;
                err = vl_sift_process_first_octave (filt, fdata) ;
            } 
            else 
                err = vl_sift_process_next_octave  (filt) ;

            if (err) 
            {
                err = VL_ERR_OK ;
                break ;
            }

            if (verbose > 1) 
            {
                printf("sift: GSS octave %d computed\n",
                        vl_sift_get_octave_index (filt));
            }

            /* run detector ............................................. */
            if (ikeys == 0) 
            {
                vl_sift_detect (filt) ;

                keys  = vl_sift_get_keypoints(filt) ;
                nkeys = vl_sift_get_nkeypoints (filt) ;
                i     = 0 ;

                if (verbose > 1) {
                    printf ("sift: detected %d (unoriented) keypoints\n", nkeys) ;
                }
            } 
            else 
                nkeys = nikeys ;

            /* for each keypoint ........................................ */
            for (; i < nkeys ; ++i) 
            {
                double                angles [4] ;
                int                   nangles ;
                VlSiftKeypoint        ik ;
                VlSiftKeypoint const *k ;

                /* obtain keypoint orientations ........................... */
                if (ikeys) 
                {
                    vl_sift_keypoint_init (filt, &ik,
                            ikeys [4 * i + 0],
                            ikeys [4 * i + 1],
                            ikeys [4 * i + 2]) ;

                    if (ik.o != vl_sift_get_octave_index (filt)) 
                    {
                        break ;
                    }
                    k = &ik ;

                    /* optionally compute orientations too */
                    if (force_orientations) {
                        nangles = vl_sift_calc_keypoint_orientations
                            (filt, angles, k) ;
                    } 
                    else 
                    {
                        angles [0] = ikeys [4 * i + 3] ;
                        nangles    = 1 ;
                    }
                } 
                else 
                {
                    k = keys + i ;
                    nangles = vl_sift_calc_keypoint_orientations
                        (filt, angles, k) ;
                }

                /* for each orientation ................................... */
                for (q = 0 ; q < (unsigned) nangles ; ++q) 
                {
                    vl_sift_pix descr [128] ;

                    /* compute descriptor (if necessary) */
                    vl_sift_calc_keypoint_descriptor(filt, descr, k, angles [q]);
                    VL_PRINT("[%.2f %.2f %.2f]\n", k->x, k->y, angles[q]);
                    int l;
                    for(l=0; l<128; l++)
                    {
                        //VL_PRINT("%d ", (vl_uint8)(512.0*descr[l]));
                    }
                    //VL_PRINT("\n");
                }
            }
        }

        /* ...............................................................
         *                                                       Finish up
         * ............................................................ */

done :
        /* release input keys buffer */
        if (ikeys) {
            free (ikeys) ;
            ikeys_size = nikeys = 0 ;
            ikeys = 0 ;
        }

        /* release filter */
        if (filt) {
            vl_sift_delete (filt) ;
            filt = 0 ;
        }

        /* release image data */
        if (fdata) {
            free (fdata) ;
            fdata = 0 ;
        }

        /* release image data */
        if (data) {
            free (data) ;
            data = 0 ;
        }

        /* close files */
        if (in) 
        {
            fclose (in) ;
            in = 0 ;
        }
    }

    /* quit */
    return exit_code ;
}
