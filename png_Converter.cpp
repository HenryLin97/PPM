// Copyright @2018 Pony AI Inc. All rights reserved.
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include<string.h>
#include<iostream>

using namespace cv;

/*This program is used to convert ppm to png file*/
int main() {
  cv::Mat image;
  // ATTENTION!!! : please use absolute path for reading the data file.
  std::string s;
  std::cout<<"Please enter the name of the ppm file: (without suffix)\n";
  std::cin>>s;
  image = imread(s + ".ppm", CV_LOAD_IMAGE_COLOR);
  namedWindow("Gray", CV_WINDOW_NORMAL);
  imwrite(s + ".png",image);
  imshow("Gray",image);
  waitKey(0);
  return 0;
}