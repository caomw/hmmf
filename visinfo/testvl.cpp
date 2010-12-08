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

        //cvCvtColor(tmp, tmp, CV_BGR2Lab);
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
        
        vl_quickshift_delete(q);
        cvReleaseImage(&tmp);
        cvReleaseImage(&newimg);
    }
    return 0;
}
