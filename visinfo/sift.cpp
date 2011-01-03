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

vector<KeyPoint> keypoints;
Mat descriptors;
cv::flann::Index *kdtree;
bool kdtree_empty = true;

double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

IplImage* process(IplImage *img) 
{
    double start = 0;
    start = get_msec();
    IplImage *newimg = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
    IplImage *toshow = cvCloneImage(img);
    cvConvertImage(img, newimg, CV_BGR2GRAY);
    
    Mat frame(newimg);

    vector<float> descp;
    keypoints.clear();
    SURF surf(5.0e3);
    surf(frame, Mat(), keypoints, descp, false);

    for(int i=0; i< keypoints.size(); i++)
    {
        cvCircle(toshow, cvPoint(keypoints[i].pt.x, keypoints[i].pt.y), 10*keypoints[i].octave, Scalar(0, 255, 0), 1);
    }
    if(keypoints.size() > 2)
    {
        descriptors = Mat(keypoints.size(), 64, CV_32F, &descp[0], 64*sizeof(float));
        printf("keypoints size: %d, descp size: %d\n", keypoints.size(), descp.size());
        
        kdtree = new cv::flann::Index(descriptors, cv::flann::KDTreeIndexParams(4));
        
        // find closest vertex to query
        Mat indices(1, 1, CV_32SC1), dists(1, 1, CV_32FC1);
        kdtree->knnSearch( descriptors.row(1), indices, dists, 1, cv::flann::SearchParams(32) );
        printf("Searching for first row...%d, %.3f\n", indices.at<int>(0, 0), dists.at<float>(0, 0));
        
    }
    return toshow;
}

int main(int argc, char** argv) 
{

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

        IplImage *toshow = process(img); 
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
    IplImage *toshow = process(img);
    imshow(window_name, toshow);
    waitKey();
#endif
    return 0;
}
