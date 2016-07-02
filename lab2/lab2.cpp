#include <iostream>
#include <ctime>
#include <cstdlib>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

Mat addSaltPepperNoise(Mat image){
	Mat noiseImage = image.clone();
	unsigned int r;
	srand((unsigned int)time(NULL));
	int nr = image.rows;
	int nc = image.cols;
	for (int i = 0; i < nr; i++){
		uchar* data = noiseImage.ptr<uchar>(i);
		for (int j = 0; j < nc; j++){
			r = rand();
			if (r < RAND_MAX*0.03){
				if (r & 0x01){
					data[j] = 255;
				}
				else{
					data[j] = 0;
				}
			}
		}
	}
	return noiseImage;
}
Mat meanFilter(Mat noiseImage){
	Mat meanImage = noiseImage.clone();
	int nr = noiseImage.rows;
	int nc = noiseImage.cols;
	for (int i = 1; i < nr - 1; i++){
		for (int j = 1; j < nc - 1; j++){
			meanImage.at<uchar>(i, j) = (
				meanImage.at<uchar>(i - 1, j - 1) + meanImage.at<uchar>(i - 1, j) + meanImage.at<uchar>(i - 1, j + 1) +
				meanImage.at<uchar>(i, j - 1) + meanImage.at<uchar>(i, j) + meanImage.at<uchar>(i, j + 1) +
				meanImage.at<uchar>(i + 1, j - 1) + meanImage.at<uchar>(i + 1, j) + meanImage.at<uchar>(i + 1, j + 1)
				) / 9;
		}
	}
	for (int j = 1; j < nc - 1; j++){
		meanImage.at<uchar>(0, j) = (meanImage.at<uchar>(0, j - 1) + meanImage.at<uchar>(0, j) + meanImage.at<uchar>(0, j + 1)
			+ meanImage.at<uchar>(1, j - 1) + meanImage.at<uchar>(1, j) + meanImage.at<uchar>(1, j + 1)) / 6.0;
		meanImage.at<uchar>(nr - 1, j) = (meanImage.at<uchar>(nr - 1, j - 1) + meanImage.at<uchar>(nr - 1, j) + meanImage.at<uchar>(nr - 1, j + 1)
			+ meanImage.at<uchar>(nr - 2, j - 1) + meanImage.at<uchar>(nr - 2, j) + meanImage.at<uchar>(nr - 2, j + 1)) / 6.0;
	}
	for (int i = 1; i < nr - 1; i++){
		meanImage.at<uchar>(i, 0) = (meanImage.at<uchar>(i - 1, 0) + meanImage.at<uchar>(i, 0) + meanImage.at<uchar>(i + 1, 0)
			+ meanImage.at<uchar>(i - 1, 1) + meanImage.at<uchar>(i, 1) + meanImage.at<uchar>(i + 1, 1)) / 6.0;
		meanImage.at<uchar>(i, nc - 1) = (meanImage.at<uchar>(i - 1, nc - 1) + meanImage.at<uchar>(i, nc - 1) + meanImage.at<uchar>(i + 1, nc - 1)
			+ meanImage.at<uchar>(i - 1, nc - 2) + meanImage.at<uchar>(i, nc - 2) + meanImage.at<uchar>(i + 1, nc - 2)) / 6.0;
	}
	meanImage.at<uchar>(0, 0) = (meanImage.at<uchar>(0, 0) + meanImage.at<uchar>(0, 1) + meanImage.at<uchar>(1, 1) + meanImage.at<uchar>(1, 0)) / 4.0;
	meanImage.at<uchar>(0, nc - 1) = (meanImage.at<uchar>(0, nc - 1) + meanImage.at<uchar>(0, nc - 2) + meanImage.at<uchar>(1, nc - 2) + meanImage.at<uchar>(1, nc - 1)) / 4.0;
	meanImage.at<uchar>(nr - 1, 0) = (meanImage.at<uchar>(nr - 1, 0) + meanImage.at<uchar>(nr - 1, 1) + meanImage.at<uchar>(nr - 1, 1) + meanImage.at<uchar>(nr - 1, 0)) / 4.0;
	meanImage.at<uchar>(nr - 1, nc - 1) = (meanImage.at<uchar>(nr - 1, nc - 1) + meanImage.at<uchar>(nr - 1, nc - 2) + meanImage.at<uchar>(nr - 2, nc - 2) + meanImage.at<uchar>(nr - 2, nc - 1)) / 4.0;
	return meanImage;
}
void swap(int &a, int &b){
	if (&a != &b){
		a ^= b ^= a ^= b;
	}
}
uchar select_middle(uchar *data, int start, int last, int nth){
	int i = start + 1, j = last - 1, pivot = data[start];
	int k;
	if (last - start < 2)
		return data[start];
	while (i <= j){
		if (data[i] < pivot){
			++i;
			continue;
		}
		if (data[j] >= pivot){
			--j;
			continue;
		}
		swap(data[i], data[j]);
	}
	swap(data[i - 1], data[start]);
	k = i - start;
	if (k == nth)  return data[i - 1];
	else if (k > nth)  return select_middle(data, start, i, nth);
	else return select_middle(data, i, last, nth - k);
}
Mat medianFilter(Mat noiseImage){
	Mat medianImage = noiseImage.clone();
	int nr = noiseImage.rows;
	int nc = noiseImage.cols;
	uchar arr[9];
	for (int i = 1; i < nr - 1; i++){
		for (int j = 1; j < nc - 1; j++){
			arr[0] = medianImage.at<uchar>(i - 1, j - 1);
			arr[1] = medianImage.at<uchar>(i - 1, j);
			arr[2] = medianImage.at<uchar>(i - 1, j + 1);
			arr[3] = medianImage.at<uchar>(i, j - 1);
			arr[4] = medianImage.at<uchar>(i, j);
			arr[5] = medianImage.at<uchar>(i, j + 1);
			arr[6] = medianImage.at<uchar>(i + 1, j - 1);
			arr[7] = medianImage.at<uchar>(i + 1, j);
			arr[8] = medianImage.at<uchar>(i + 1, j + 1);
			medianImage.at<uchar>(i, j) = select_middle(arr, 0, 9, 5);
		}
	}
	for (int j = 1; j < nc - 1; j++){
		arr[0] = medianImage.at<uchar>(0, j - 1);
		arr[1] = medianImage.at<uchar>(0, j);
		arr[2] = medianImage.at<uchar>(0, j + 1);
		arr[3] = medianImage.at<uchar>(1, j - 1);
		arr[4] = medianImage.at<uchar>(1, j);
		arr[5] = medianImage.at<uchar>(1, j + 1);
		medianImage.at<uchar>(0, j) = select_middle(arr, 0, 6, 3);
		arr[0] = medianImage.at<uchar>(nr - 1, j - 1);
		arr[1] = medianImage.at<uchar>(nr - 1, j);
		arr[2] = medianImage.at<uchar>(nr - 1, j + 1);
		arr[3] = medianImage.at<uchar>(nr - 2, j - 1);
		arr[4] = medianImage.at<uchar>(nr - 2, j);
		arr[5] = medianImage.at<uchar>(nr - 2, j + 1);
		medianImage.at<uchar>(nr - 1, j) = select_middle(arr, 0, 6, 3);
	}
	for (int i = 1; i < nr - 1; i++){
		arr[0] = medianImage.at<uchar>(i - 1, 0);
		arr[1] = medianImage.at<uchar>(i, 0);
		arr[2] = medianImage.at<uchar>(i + 1, 0);
		arr[3] = medianImage.at<uchar>(i - 1, 1);
		arr[4] = medianImage.at<uchar>(i, 1);
		arr[5] = medianImage.at<uchar>(i + 1, 1);
		medianImage.at<uchar>(i, 0) = select_middle(arr, 0, 6, 3);
		arr[0] = medianImage.at<uchar>(i - 1, nc - 1);
		arr[1] = medianImage.at<uchar>(i, nc - 1);
		arr[2] = medianImage.at<uchar>(i + 1, nc - 1);
		arr[3] = medianImage.at<uchar>(i - 1, nc - 2);
		arr[4] = medianImage.at<uchar>(i, nc - 2);
		arr[5] = medianImage.at<uchar>(i + 1, nc - 2);
		medianImage.at<uchar>(i, nc - 1) = select_middle(arr, 0, 6, 3);
	}
	arr[0] = medianImage.at<uchar>(0, 0);
	arr[1] = medianImage.at<uchar>(0, 1);
	arr[2] = medianImage.at<uchar>(1, 1);
	arr[3] = medianImage.at<uchar>(1, 0);
	medianImage.at<uchar>(0, 0) = select_middle(arr, 0, 4, 2);
	arr[0] = medianImage.at<uchar>(0, nc - 1);
	arr[1] = medianImage.at<uchar>(0, nc - 2);
	arr[2] = medianImage.at<uchar>(1, nc - 2);
	arr[3] = medianImage.at<uchar>(1, nc - 1);
	medianImage.at<uchar>(0, nc - 1) = select_middle(arr, 0, 4, 2);
	arr[0] = medianImage.at<uchar>(nr - 1, 0);
	arr[1] = medianImage.at<uchar>(nr - 1, 1);
	arr[2] = medianImage.at<uchar>(nr - 2, 1);
	arr[3] = medianImage.at<uchar>(nr - 2, 0);
	medianImage.at<uchar>(nr - 1, 0) = select_middle(arr, 0, 4, 2);
	arr[0] = medianImage.at<uchar>(nr - 1, nc - 1);
	arr[1] = medianImage.at<uchar>(nr - 1, nc - 2);
	arr[2] = medianImage.at<uchar>(nr - 2, nc - 2);
	arr[3] = medianImage.at<uchar>(nr - 2, nc - 1);
	medianImage.at<uchar>(nr - 1, nc - 1) = select_middle(arr, 0, 4, 2);
	return medianImage;
}
int main(int argc,char** argv){
	String s;
	if (argc > 1) {
		s = argv[1];
	}
	else {
		s = "../source/lena.bmp";
	}
	Mat image = imread(s, CV_8U);
	if (image.empty()) {
		cout << "Can't open image " << s << endl;
		return -1;
	}
	namedWindow("origin image");
	imshow("origin image", image);
	Mat noiseImage = addSaltPepperNoise(image);
	namedWindow("noise image");
	imshow("noise image", noiseImage);
	Mat meanImage = meanFilter(noiseImage);
	namedWindow("mean filter image");
	imshow("mean filter image", meanImage);
	Mat medianImage = medianFilter(noiseImage);
	namedWindow("median filter image");
	imshow("median filter image", medianImage);
	waitKey();
	return 0;
}