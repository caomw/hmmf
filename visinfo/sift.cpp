#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc_c.h"

#include <iostream>
#include <vector>

using namespace cv;
using namespace std;

int threshold = 50;

int process(CvCapture *capture) 
{
	string window_name = "video | q or esc to quit";
	cout << "press q or esc to quit" << endl;
	namedWindow(window_name, CV_WINDOW_KEEPRATIO); //resizable window;

	for (;;) 
    {
        IplImage *img = 0;

		if (!cvGrabFrame(capture))
        {
            printf("Couldn't capture frame\n");
			return 1;
        }
        img = cvRetrieveFrame(capture);
        IplImage *newimg = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
        cvConvertImage(img, newimg, CV_BGR2GRAY);
        Mat frame(newimg);

        vector<KeyPoint> keypoints;
        FAST(frame, keypoints, threshold, true);
        IplImage toshow = frame;

        int count = 0;
        for(int i=0; i< keypoints.size(); i++)
        {
            count++;
            //cout<<i<<" "<<keypoints[i].pt.x<<" "<<keypoints[i].pt.y<<" "<<endl;
            cvCircle(&toshow, cvPoint(keypoints[i].pt.x, keypoints[i].pt.y), 5, Scalar(0, 0, 255), 2);
        }

		imshow(window_name, &toshow);
        cout<<count<<endl;
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
	return 0;
}

int main(int argc, char** argv) 
{
    if(argc != 2)
        printf("Usage: ./sift <thresh>\n");
    else
    {
        threshold = atoi(argv[1]);
        CvCapture *capture = cvCaptureFromCAM(0);
        return process(capture);
    }
}
