extern "C"
{
#include <vl/generic.h>
#include <vl/quickshift.h>
#include <vl/pgm.h>
#include <vl/sift.h>
}
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <stdio.h>

int height;
int width;

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        VL_PRINTF("Usage: ./testvl <img.pgm>");
        return 0;
    }
    
    if (argc == 2)
    {
        IplImage *img = cvLoadImage(argv[1], 1);        // color
        height = img->height;
        width = img->width;
        VL_PRINTF("Depth: %d Size: %d\n", img->depth, img->imageSize);

        //cvCvtColor(img, img, CV_BGR2Lab);
        cvNamedWindow("opencv");
        cvShowImage("opencv", img);
        cvWaitKey();
        
        VlSiftFilt *sift = vl_sift_new(width, height, 4, 5, 0);
        int i=0, err = ~VL_ERR_EOF;
        while( err != VL_ERR_EOF)
        {
            err = vl_sift_process_next_octave(sift);
            printf("process octave %d\n", i++);
        }
        vl_sift_detect(sift);
        VlSiftKeypoint const *keypoints = vl_sift_keypoints(sift);
        


        cvReleaseImage(&img);
    }
    return 0;
}
