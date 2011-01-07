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

// each element stores all the keypoints found in that scene
vector< vector<KeyPoint> > keypoints;
// correspondingly the descriptors of keypoints[i]
vector<Mat> descriptors;
cv::flann::Index *kdtree;
bool empty_tree = true, note_this = false;

Mat indices(1, 2, CV_32SC1), dists(1, 2, CV_32FC1);


double get_msec()
{
    struct timeval start;
    gettimeofday(&start, NULL);
    return start.tv_sec*1000 + start.tv_usec/1000.0;
}

IplImage* process(IplImage *img) 
{
    double start = 0;
    int match_count = 0;
    start = get_msec();
    IplImage *newimg = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
    IplImage *toshow = cvCloneImage(img);
    cvConvertImage(img, newimg, CV_BGR2GRAY);
    
    Mat frame(newimg);

    vector<float> descpVec;
    vector<KeyPoint> key;
    SURF surf(5000, 4, 2, true);
    surf(frame, Mat(), key, descpVec, false);
    Mat descpMat = Mat(key.size(), 128, CV_32F, &descpVec[0], 128*sizeof(float))*1000;

    //save the keypoints and descriptors
    if(note_this == true)
    {
        note_this = false;
        keypoints.push_back(key);
        descriptors.push_back(descpMat);    
    }
    else
    {
        cout<<"Matching...\n";
        BruteForceMatcher<L2<float> > matcher;
        for(int i=0; i<descriptors.size(); i++)
        {
            vector<DMatch> matches;
            matcher.match(descriptors[i], descpMat, matches);
            for(int j=0; j<matches.size(); j++)
                cout<<matches[j].distance<<endl;
        }
    }
    for(int i=0; i< key.size(); i++)
    {
        cvCircle(toshow, cvPoint(key[i].pt.x, key[i].pt.y), 10*key[i].octave, Scalar(0, 255, 0), 1);
    }
    
    /*
    if( (keypoints.size() > 0) && (empty_tree == true) )
    {
        empty_tree = false;
        printf("1: keypoints size: %d, descp size: %d\n", keypoints.size(), descp.size());
        
        kdtree = new cv::flann::Index(descriptors, cv::flann::KDTreeIndexParams(4));
        
    }
    else if ( (keypoints.size() > 0) && (empty_tree != true) )
    {
        //use old tree
        Mat descriptors_curr = Mat(keypoints.size(), 128, CV_32F, &descp[0], 128*sizeof(float))*1000;
        printf("------------------------\nKeypoints size: %d \n", keypoints.size());

        for(int i=0; i<descriptors_curr.rows; i++)
        {
            // find closest vertex to row(i)
            kdtree->knnSearch( descriptors_curr.row(i), indices, dists, 2, cv::flann::SearchParams(64) );
            //cout<<"Indices :"<<indices<<endl;
            //cout<<"Dists: "<<dists<<endl;

            if( dists.at<float>(0, 0) > 100000.0 )
            {
                if(first_few > 0)
                {
                    //printf("Adding new row: %d\n", i);
                    Mat new_descp_mat = Mat( descriptors.rows + 1, 128, CV_32F, Scalar(0.0) );
                    new_descp_mat.rowRange(0, descriptors.rows) = descriptors.rowRange(0, descriptors.rows);        // copy the descriptors matrix fully
                    new_descp_mat.row(descriptors.rows) = descriptors_curr.row(i);

                    descriptors = new_descp_mat;                                                                    // copy the appended matrix to descriptors
                    kdtree = new cv::flann::Index(descriptors, cv::flann::KDTreeIndexParams(4));                    // create new kdtree
                }
            }
            else
            {
                //printf("Matching row: %d\n", i);
                match_count++;
            }
        }
        printf("Match stats -- Descriptor size: %d Match count: %d Scene match: %.2f\n\n", descriptors.rows, match_count, ((double)match_count/keypoints.size()) );
    }
    */
    return toshow;
}

int main(int argc, char** argv) 
{

    string window_name = "video | q or esc to quit";
    cout << "press q or esc to quit" << endl;
    namedWindow(window_name); //resizable window;
#if 0
    CvCapture *capture = cvCaptureFromCAM(0);
    waitKey(2000);

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

        char key = (char)waitKey(10); //delay N millis, usually long enough to display and capture input
        switch (key) 
        {
            case 'w':
                note_this = true;
            case 'q':
            case 'Q':
            case 27: //escape key
                return 0;
            default:
                break;
        }
    }
#else
    IplImage *img = cvLoadImage("box.pgm", 1);
    note_this = true;
    IplImage *toshow = process(img);
    IplImage *img1 = cvLoadImage("scene.pgm", 1);
    IplImage *toshow1 = process(img1);
    imshow(window_name, toshow);
    waitKey();
    imshow(window_name, toshow1);
    waitKey();
#endif

    /*
     * cv::Mat test
    double data[3][3] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Mat M = Mat(3, 3, CV_64F, data);
    Mat r = M.row(0);
    printf("r: [%d %d]\n", r.rows, r.cols);
    for(int i=0; i<r.rows; i++)
    {
        for(int j=0; j<r.cols; j++)
        {
            printf("%.1f ", r.at<double>(i, j));
        }
        printf("\n");
    }
*/

    return 0;
}
