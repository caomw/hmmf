#include "Image.h"
#include "Exception.h"
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "quickshift_common.h"
#include <sys/time.h>
#include <stdlib.h>

double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

void write_image(image_t im, const char * filename)
{
    /********** Copy from matlab style **********/
    Image IMGOUT(im.K > 1 ? Image::RGB : Image::L, im.N2, im.N1);
    for(int k = 0; k < im.K; k++)
        for(int col = 0; col < im.N2; col++)
            for(int row = 0; row < im.N1; row++)
            {
                /* Row transpose */
                unsigned char * pt = IMGOUT.getPixelPt(col, im.N1-1-row);
                /* scale 0-255 */
                pt[k] = (unsigned char) (im.I[row + col*im.N1 + k*im.N1*im.N2]/32*255);
            }

    /********** Write image **********/
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        throw Exception("Could not open the file");
    }
    ofs<<IMGOUT;
}

image_t imseg(image_t im, int * flatmap)
{
    /********** Mean Color **********/
    float * meancolor = (float *) calloc(im.N1*im.N2*im.K, sizeof(float)) ;
    float * counts    = (float *) calloc(im.N1*im.N2, sizeof(float)) ;

    int roots = 0;
    for (int p = 0; p < im.N1*im.N2; p++)
    {
        counts[flatmap[p]]++;
        for (int k = 0; k < im.K; k++)
            meancolor[flatmap[p] + k*im.N1*im.N2] += im.I[p + k*im.N1*im.N2];
        
        if (flatmap[p] == p)
            roots++;
    }
    printf("Roots: %d\n", roots);

    int nonzero = 0;
    for (int p = 0; p < im.N1*im.N2; p++)
    {
        if (counts[p] > 0)
        {
            nonzero++;
            for (int k = 0; k < im.K; k++)
                meancolor[p + k*im.N1*im.N2] /= counts[p];
        }
    }
    if (roots != nonzero)
        printf("Nonzero: %d\n", nonzero);
    assert(roots == nonzero);


    /********** Create output image **********/
    image_t imout = im;
    imout.I = (float *) calloc(im.N1*im.N2*im.K, sizeof(float));
    for (int p = 0; p < im.N1*im.N2; p++)
        for (int k = 0; k < im.K; k++)
            imout.I[p + k*im.N1*im.N2] = meancolor[flatmap[p] + k*im.N1*im.N2];

    free(meancolor);
    free(counts);

    return imout;
}

int * map_to_flatmap(float * map, unsigned int size)
{
    /********** Flatmap **********/
    int *flatmap      = (int *) malloc(size*sizeof(int)) ;
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

void image_to_matlab(Image & IMG, image_t & im)
{
    /********** Convert image to MATLAB style representation **********/
    im.N1 = IMG.getHeight();
    im.N2 = IMG.getWidth();
    im.K  = IMG.getPixelSize();
    im.I = (float *) calloc(im.N1*im.N2*im.K, sizeof(float));
    for(int k = 0; k < im.K; k++)
        for(int col = 0; col < im.N2; col++)
            for(int row = 0; row < im.N1; row++)
            {
                unsigned char * pt = IMG.getPixelPt(col, im.N1-1-row);
                im.I[row + col*im.N1 + k*im.N1*im.N2] = 32. * pt[k] / 255.; // Scale 0-32
            }
}

int main(int argc, char ** argv)
{
    //Use command-line specified CUDA device, otherwise use device with highest Gflops/s
    float sigma = 2, tau = 10;
    char *file = "flowers2.pnm";
    if(argc == 4)
    {
        file = argv[1];
        sigma = atoi(argv[2]);
        tau = atoi(argv[3]);
    }

    /********** Read image **********/
    Image IMG;
    char outfile[1024];

    std::ifstream ifs(file, std::ios::binary);
    if (!ifs) 
    {
        throw Exception("Could not open the file");
    }
    ifs>>IMG;
    image_t im;

    image_to_matlab(IMG, im);

    float *map, *E, *gaps;
    int * flatmap;
    image_t imout;

    map          = (float *) calloc(im.N1*im.N2, sizeof(float)) ;
    gaps         = (float *) calloc(im.N1*im.N2, sizeof(float)) ;
    E            = (float *) calloc(im.N1*im.N2, sizeof(float)) ;

    /********** Quick shift **********/
    double start = get_msec();
    quickshift(im, sigma, tau, map, gaps, E);
    printf("quickshift: [%3.3f ms]\n", get_msec() - start);

    /* Consistency check */
    for(int p = 0; p < im.N1*im.N2; p++)
        if(map[p] == p) 
            assert(gaps[p] == INF);

    flatmap = map_to_flatmap(map, im.N1*im.N2);
    imout = imseg(im, flatmap);

    sprintf(outfile, "%s", file);
    char * c = strrchr(outfile, '.');
    if(c) *c = '\0';
    sprintf(outfile, "%s-cpu.pnm", outfile); 
    write_image(imout, outfile);
    

    /********** Cleanup **********/
    free(flatmap);
    free(imout.I);
    free(im.I);

    free(map);
    free(E);
    free(gaps);
};
