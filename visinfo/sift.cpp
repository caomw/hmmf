#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/core/internal.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc_c.h"

#include <iostream>
#include <vector>
#include <sys/time.h>
#include <time.h>

using namespace cv;
using namespace std;

int threshold = 50;

double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

IplImage* process(IplImage *img, vector<KeyPoint> keypoints) 
{
    IplImage *newimg = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
    IplImage *toshow = cvCloneImage(img);
    cvConvertImage(img, newimg, CV_BGR2GRAY);

    int count = 0;
    for(int octave=0; octave<1; octave++)
    {
        IplImage *new_octave = cvCreateImage( cvSize(newimg->width/pow(2, octave), newimg->height/pow(2,octave) ), IPL_DEPTH_8U, 1);
        cvResize(newimg, new_octave);

        Mat frame(new_octave);
        keypoints.clear();
        FAST(frame, keypoints, threshold, true);

        for(int i=0; i< keypoints.size(); i++)
        {
            count++;
            keypoints[i].pt.x *= pow(2, octave);
            keypoints[i].pt.y *= pow(2, octave);
            keypoints[i].octave = pow(2, octave);

            //cout<<i<<" "<<keypoints[i].pt.x<<" "<<keypoints[i].pt.y<<" "<<endl;
            cvCircle(toshow, cvPoint(keypoints[i].pt.x, keypoints[i].pt.y), 10*keypoints[i].octave, Scalar(0, 255, 0), 1);
        }
        
        Mat descriptors;
        SIFT sift(3.0);
        sift(frame, Mat(), keypoints, descriptors, true);
        FileStorage fs("descp.xml", FileStorage::WRITE);
        fs<<"descp"<<descriptors;
    }
    cout<<count<<endl;

    return toshow;
}

int main(int argc, char** argv) 
{
    if(argc != 2)
        printf("Usage: ./sift <thresh>\n");
    else
    {
        threshold = atoi(argv[1]);
        
        string window_name = "video | q or esc to quit";
        cout << "press q or esc to quit" << endl;
        namedWindow(window_name); //resizable window;
#if 0
        CvCapture *capture = cvCaptureFromCAM(0);

        for (;;) 
        {
            IplImage *img = 0;

            if (!cvGrabFrame(capture))
            {
                printf("Couldn't capture frame\n");
                return 1;
            }
            img = cvRetrieveFrame(capture);
            
            vector<KeyPoint> keypoints;
            IplImage *toshow = process(img, keypoints); 

            imshow(window_name, toshow);
            char key = (char)waitKey(5); //delay N millis, usually long enough to display and capture input
            switch (key) 
            {
                case 'q':
                case 'Q':
                case 27: //escape key
                    return 0;
                default:
                    break;
            }
        }
#else
        IplImage *img = cvLoadImage("a.jpg", 1);
        vector<KeyPoint> keypoints;
        IplImage *toshow = process(img, keypoints); 
        imshow(window_name, toshow);
        waitKey();
#endif
    }
    return 0;
}
