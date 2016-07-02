#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

void lab11(Mat image){
	double a, b, t;
	Mat outImage;
	cout << "input a b" << endl;
	cin >> a >> b;
	outImage.create(image.size(), image.type());
	int nr = image.rows;
	int nl = image.cols*image.channels();
	for (int k = 0; k < nr; k++){
		const uchar* inData = image.ptr<uchar>(k);
		uchar* outData = outImage.ptr<uchar>(k);
		for (int i = 0; i < nl; i++){
			t = a*inData[i] + b + 0.5;
			if (t < 0)
				outData[i] = 0;
			else if (t>255)
				outData[i] = 255;
			else
				outData[i] = (uchar)t;
		}
	}
	namedWindow("origin");
	imshow("origin", image);
	namedWindow("lab11");
	imshow("lab11", outImage);
	waitKey();
	cvDestroyWindow("origin");
	cvDestroyWindow("lab11");
}
void lab12(Mat image){
	int x1, y1, x2, y2;
	double t;
	Mat outImage;
	cout << "input x1 y1 x2 y2" << endl;
	cin >> x1 >> y1 >> x2 >> y2;
	outImage.create(image.size(), image.type());
	int nr = image.rows;
	int nl = image.cols*image.channels();
	for (int k = 0; k < nr; k++){
		const uchar* inData = image.ptr<uchar>(k);
		uchar* outData = outImage.ptr<uchar>(k);
		for (int i = 0; i < nl; i++){
			if (inData[i] < x1){
				t = (double)y1*inData[i] / x1 + 0.5;
			}
			else if (inData[i]>x2){
				t = (double)(255 - y2)*(inData[i] - x2) / (255 - x2) + y2 + 0.5;
			}
			else{
				t = (double)(y2 - y1)*(inData[i] - x1) / (x2 - x1) + y1 + 0.5;
			}
			if (t < 0)
				outData[i] = 0;
			else if (t>255)
				outData[i] = 255;
			else
				outData[i] = (uchar)t;
		}
	}
	namedWindow("origin");
	imshow("origin", image);
	namedWindow("lab12");
	imshow("lab12", outImage);
	waitKey();
	cvDestroyWindow("origin");
	cvDestroyWindow("lab12");
}
void lab13(Mat image){
	int upper, lower;
	cout << "input lower upper" << endl;
	cin >> lower >> upper;
	int histSize = 256;
	unsigned int *hist = (unsigned int *)calloc(histSize, sizeof(unsigned int));
	int nr = image.rows;
	int nl = image.cols*image.channels();
	double range = (double)(upper - lower) / histSize;
	for (int k = 0; k < nr; k++){
		const uchar* inData = image.ptr<uchar>(k);
		for (int i = 0; i < nl; i++){
			if (inData[i] >= lower&&inData[i] < upper){
				hist[(int)((inData[i] - lower) / range + 0.5)]++;
			}
		}
	}
	cout << endl;
	unsigned max = 0, min = 0xFFFFFFFF;
	for (int i = 0; i < histSize; i++){
		max = max < hist[i] ? hist[i] : max;
		min = min > hist[i] ? hist[i] : min;
	}
	int hv = max / 0.9;
	Mat histogram(500, 500, CV_8U, Scalar(255));
	for (int i = 0; i < histSize; i++){
		rectangle(histogram, Point(i * 500 / histSize, 500), Point((i + 1) * 500 / histSize, 500 - hist[i] * 500 / hv), Scalar(128), -1, 8, 0);
	}
	namedWindow("origin");
	imshow("origin", image);
	namedWindow("lab13");
	imshow("lab13", histogram);
	waitKey();
	cvDestroyWindow("lab13");
	cvDestroyWindow("origin");
}
void lab14(Mat image){
	int L = 256;
	unsigned long *hist = (unsigned long*)calloc(L, sizeof(unsigned long));
	int nr = image.rows;
	int nc = image.cols;
	for (int k = 0; k < nr; k++){
		const uchar* inData = image.ptr<uchar>(k);
		for (int i = 0; i < nc; i++){
			hist[inData[i] / (256 / L)]++;
		}
	}
	unsigned int max = 0, min = 0xFFFFFFFF;
	for (int i = 0; i < L; i++){
		max = max < hist[i] ? hist[i] : max;
		min = min > hist[i] ? hist[i] : min;
	}
	int hv = max / 0.9;
	Mat oldHistogram(512, 512, CV_8U, Scalar(255));
	for (int i = 0; i < L; i++){
		rectangle(oldHistogram, Point(i * 512 / L, 512), Point((i + 1) * 512 / L, 512 - hist[i] * 512 / hv), Scalar(128), -1, 8, 0);
	}
	unsigned long sum = nr*nc;
	int temp = 0;
	int *f = (int *)calloc(L, sizeof(int));
	for (int i = 0; i < L; i++){
		temp += hist[i];
		f[i] = (int)((L - 1)*(double)temp / sum + 0.5) + (256 / L) / 2;
	}
	Mat newImage(image.size(), image.type(), Scalar(255));
	for (int k = 0; k < nr; k++){
		const uchar* inData = image.ptr<uchar>(k);
		uchar* outData = newImage.ptr<uchar>(k);
		for (int i = 0; i < nc; i++){
			outData[i] = f[inData[i] / (256 / L)];
		}
	}
	for (int k = 0; k < nr; k++){
		const uchar* inData = newImage.ptr<uchar>(k);
		for (int i = 0; i < nc; i++){
			hist[inData[i] / (256 / L)]++;
		}
	}
	max = 0, min = 0xFFFFFFFF;
	for (int i = 0; i < L; i++){
		max = max < hist[i] ? hist[i] : max;
		min = min > hist[i] ? hist[i] : min;
	}
	hv = max / 0.9;
	Mat newHistogram(512, 512, CV_8U, Scalar(255));
	for (int i = 0; i < L; i++){
		rectangle(newHistogram, Point(i * 512 / L, 512), Point((i + 1) * 512 / L, 512 - hist[i] * 512 / hv), Scalar(128), -1, 8, 0);
	}
	namedWindow("oldImage");
	imshow("oldImage", image);
	namedWindow("oldHistogram");
	imshow("oldHistogram", oldHistogram);
	namedWindow("newImage");
	imshow("newImage", newImage);
	namedWindow("newHistogram");
	imshow("newHistogram", newHistogram);
	waitKey();
	cvDestroyWindow("oldImage");
	cvDestroyWindow("oldHistogram");
	cvDestroyWindow("newImage");
	cvDestroyWindow("newHistogram");
}
int main(int argc, char** argv){
	String s;
	if (argc > 1) {
		s = argv[1];
	}
	else {
		s = "../source/pout.bmp";
	}
	Mat image = imread(s, CV_8U);
	if (image.empty()) {
		cout << "Can't open image " << s << endl;
		return -1;
	}
	int select = 5;
	do{
		cout << "1.灰度的线性变换\n2.灰度拉伸\n3.灰度直方图\n4.直方图均衡\n0.exit\ninput your select(1-4,0):" << endl;
		cin >> select;
		switch (select){
		case 1:lab11(image); break;
		case 2:lab12(image); break;
		case 3:lab13(image); break;
		case 4:lab14(image); break;
		case 0:break;
		default:break;
		}
	} while (select != 0);
	return 0;
}