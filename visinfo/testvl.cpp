extern "C"
{
#include <vl/generic.h>
#include <vl/quickshift.h>
#include <vl/pgm.h>
}
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <stdio.h>


IplImage * imseg(IplImage *im, int * flatmap, int size)
{
    /********** Mean Color **********/
    double * meancolor = (double *) calloc(size*3, sizeof(double)) ;
    double * counts    = (double *) calloc(size, sizeof(double)) ;
    
    for (int p = 0; p < size; p++)
    {
        counts[flatmap[p]]++;
        for (int k = 0; k < 3; k++)
            meancolor[flatmap[p] + k*size] += im->imageData[(p + k*size)*4];
    }

    int roots = 0;
    for (int p = 0; p < size; p++)
    {
        if (flatmap[p] == p)
            roots++;
    }
    printf("Roots: %d\n", roots);
    int nonzero = 0;
    for (int p = 0; p < size; p++)
    {
        if (counts[p] > 0)
        {
            nonzero++;
            for (int k = 0; k < 3; k++)
                meancolor[p + k*size] /= counts[p];
        }
    }
    if (roots != nonzero)
        printf("Nonzero: %d\n", nonzero);
    assert(roots == nonzero);

    /*
    image_t imout = im;
    imout.I = (float *) calloc(im.N1*im.N2*im.K, sizeof(float));
    for (int p = 0; p < im.N1*im.N2; p++)
        for (int k = 0; k < im.K; k++)
            imout.I[p + k*im.N1*im.N2] = meancolor[flatmap[p] + k*im.N1*im.N2];
    */
    IplImage *tmp = cvCloneImage(im);
    for (int p = 0; p < size; p++)
        for (int k = 0; k < 3; k++)
            tmp->imageData[p + k*size] = meancolor[flatmap[p] + k*size];
     
    free(meancolor);
    free(counts);
    
    return tmp;
}

int * map_to_flatmap(int * map, unsigned int size)
{
    /********** Flatmap **********/
    int *flatmap = (int *) malloc(size*sizeof(int)) ;
    for (unsigned int p = 0; p < size; p++)
    {
        flatmap[p] = map[p];
    }

    bool changed = true;
    while (changed)
    {
        changed = false;
        for (unsigned int p = 0; p < size; p++)
        {
            changed = changed || (flatmap[p] != flatmap[flatmap[p]]);
            flatmap[p] = flatmap[flatmap[p]];
        }
    }

    /* Consistency check */
    for (unsigned int p = 0; p < size; p++)
        assert(flatmap[p] == flatmap[flatmap[p]]);

    return flatmap;
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        VL_PRINTF("Usage: ./testvl <img.pgm>");
        return 0;
    }
    
    if (argc == 2)
    {
        IplImage *tmp = cvLoadImage(argv[1], 1);        // color
        int height = tmp->height;
        int width = tmp->width;
        VL_PRINTF("Depth: %d Size: %d\n", tmp->depth, tmp->imageSize);

        cvCvtColor(tmp, tmp, CV_BGR2Lab);
        cvNamedWindow("opencv");
        cvShowImage("opencv", tmp);
        cvWaitKey();

        IplImage *newimg = cvCreateImage(cvGetSize(tmp), IPL_DEPTH_32F, 3);
        cvConvertScale(tmp, newimg, 1/255.0);
        cvShowImage("opencv", newimg);
        cvWaitKey();
    
        VL_PRINTF("Start quickshift\n");
        
        VlQS *q = vl_quickshift_new( (double const *)(newimg->imageData), height, width, 1);
        vl_quickshift_set_kernel_size(q, 2.0);
        vl_quickshift_set_max_dist(q, 10.0);
        vl_quickshift_process(q);
        VL_PRINTF("End quickshift\n");
        
        vl_qs_type *gap = vl_quickshift_get_dists(q);
        int *map = vl_quickshift_get_parents(q);
        vl_qs_type *density = vl_quickshift_get_density(q);

        int* flatmap = map_to_flatmap(map, height*width);
        IplImage *seg = imseg(newimg, flatmap, height*width);
        cvShowImage("opencv", seg);
        cvWaitKey();

        vl_quickshift_delete(q);
        cvReleaseImage(&tmp);
        cvReleaseImage(&newimg);
        cvReleaseImage(&seg);
    }
    return 0;
}
