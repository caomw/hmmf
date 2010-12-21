#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <vector>

using namespace cv;
using namespace std;


int process(VideoCapture& capture) 
{
	string window_name = "video | q or esc to quit";
	cout << "press q or esc to quit" << endl;
	namedWindow(window_name, CV_WINDOW_KEEPRATIO); //resizable window;
	Mat frame;
	for (;;) 
    {
		capture >> frame;
		if (frame.empty())
			continue;
		imshow(window_name, frame);
        cout<<frame.cols<<" "<<frame.rows<<endl;
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

int main(int ac, char** av) 
{
	if (ac != 2) 
	{
	    cout<< "Pass video input\n";
		return 1;
	}
	std::string arg = av[1];
	VideoCapture capture(arg);
	if (!capture.isOpened())
		capture.open(atoi(arg.c_str()));
	if (!capture.isOpened()) 
	{
		cerr << "Failed to open a video device or video file!\n" << endl;
		return 1;
	}
	return process(capture);
}
