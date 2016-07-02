#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

Mat lab31(Mat image){
	Mat RobertsImage = image.clone();
	int nr = image.rows;
	int nc = image.cols;
	for (int i = 0; i < nr-1; i++){
		const uchar* inData1 = image.ptr<uchar>(i);
		const uchar* inData2 = image.ptr<uchar>(i + 1);
		uchar* outDate = RobertsImage.ptr<uchar>(i);
		for (int j = 0; j < nc - 1; j++){
			outDate[j] = abs(inData1[j] - inData2[j + 1]) + abs(inData2[j] - inData1[j + 1]);
		}
	}
	for (int i = 0; i < nr - 1; i++){
		RobertsImage.at<uchar>(i, nc - 1) = RobertsImage.at<uchar>(i, nc - 2);
	}
	for (int j = 0; j < nc; j++){
		RobertsImage.at<uchar>(nr - 1, j) = RobertsImage.at<uchar>(nr - 2, j);
	}
	return RobertsImage;
}
Mat lab32(Mat image){
	Mat SobelImage = image.clone();
	int nr = image.rows;
	int nc = image.cols;
	for (int i = 1; i < nr - 1; i++){
		const uchar* inData0 = image.ptr<uchar>(i - 1);
		const uchar* inData1 = image.ptr<uchar>(i);
		const uchar* inData2 = image.ptr<uchar>(i + 1);
		uchar* outDate = SobelImage.ptr<uchar>(i);
		for (int j = 1; j < nc - 1; j++){
			outDate[j] = max(
				abs(inData2[j - 1] + 2 * inData2[j] + inData2[j + 1] - inData0[j - 1] - 2 * inData0[j] + inData0[j + 1]),
				abs(inData0[j + 1] + 2 * inData1[j + 1] + inData2[j + 1] - inData0[j - 1] - 2 * inData1[j - 1] - inData2[j - 1]));
		}
	}
	return SobelImage;
}
Mat lab33(Mat image){
	Mat PrewittImage = image.clone();
	int nr = image.rows;
	int nc = image.cols;
	for (int i = 1; i < nr - 1; i++){
		const uchar* inData0 = image.ptr<uchar>(i - 1);
		const uchar* inData1 = image.ptr<uchar>(i);
		const uchar* inData2 = image.ptr<uchar>(i + 1);
		uchar* outDate = PrewittImage.ptr<uchar>(i);
		for (int j = 1; j < nc - 1; j++){
			outDate[j] = max(
				abs(inData2[j - 1] + inData2[j] + inData2[j + 1] - inData0[j - 1] - inData0[j] - inData0[j + 1]),
				abs(inData0[j - 1] + inData1[j - 1] + inData2[j - 1] - inData0[j + 1] - inData1[j + 1] - inData2[j + 1]));
		}
	}
	return PrewittImage;
}
Mat lab14(Mat image) {
	int L = 256;
	unsigned long *hist = (unsigned long *)calloc(L, sizeof(unsigned long));
	int nr = image.rows;
	int nc = image.cols;
	for (int k = 0; k < nr; k++) {
		const uchar *inData = image.ptr<uchar>(k);
		for (int i = 0; i < nc; i++) {
			hist[inData[i]]++;
		}
	}
	unsigned long sum = nr * nc;
	int temp = 0;
	int *f = (int *)calloc(L, sizeof(int));
	for (int i = 0; i < L; i++) {
		temp += hist[i];
		f[i] = (int)((L - 1) * (double)temp / sum + 0.5) + (256 / L) / 2;
	}
	Mat newImage(image.size(), image.type(), Scalar(255));
	for (int k = 0; k < nr; k++) {
		const uchar *inData = image.ptr<uchar>(k);
		uchar * outData = newImage.ptr<uchar>(k);
		for (int i = 0; i < nc; i++) {
			outData[i] = f[inData[i] / (256 / L)];
		}
	}
	return newImage;
}
Mat lab341(Mat image){
	Mat LaplaceImage = image.clone();
	int nr = image.rows;
	int nc = image.cols;
	int max = -1024, min = 1023;
	int **matrix = new int*[nr];
	for (int i = 0; i < nr; i++){
		matrix[i] = new int[nc];
	}
	for (int i = 1; i < nr - 1; i++){
		const uchar* inData0 = image.ptr<uchar>(i - 1);
		const uchar* inData1 = image.ptr<uchar>(i);
		const uchar* inData2 = image.ptr<uchar>(i + 1);
		int* outDate = matrix[i];
		for (int j = 1; j < nc - 1; j++){
			outDate[j] = inData0[j] + inData1[j - 1] + inData1[j + 1] + inData2[j] - 4 * inData1[j];
			if (max < outDate[j]){
				max = outDate[j];
			}
			if (min > outDate[j]){
				min = outDate[j];
			}
		}
	}
	double size = (max - min) / 256.0;
	for (int i = 1; i < nr - 1; i++){
		int* inDate = matrix[i];
		uchar* outDate = LaplaceImage.ptr<uchar>(i);
		for (int j = 1; j < nc - 1; j++){
			outDate[j] = (inDate[j] - min) / size;
		}
	}
//	return lab14(LaplaceImage);
	return LaplaceImage;
}
Mat lab342(Mat image){
	Mat LaplaceImage = image.clone();
	int nr = image.rows;
	int nc = image.cols;
	int max = -2048, min = 2047;
	int **matrix = new int*[nr];
	for (int i = 0; i < nr; i++){
		matrix[i] = new int[nc];
	}
	for (int i = 1; i < nr - 1; i++){
		const uchar* inData0 = image.ptr<uchar>(i - 1);
		const uchar* inData1 = image.ptr<uchar>(i);
		const uchar* inData2 = image.ptr<uchar>(i + 1);
		int* outDate = matrix[i];
		for (int j = 1; j < nc - 1; j++){
			outDate[j] = inData0[j - 1] + inData0[j] + inData0[j + 1] + inData1[j - 1] + inData1[j + 1] + inData2[j - 1] + inData2[j] + inData2[j + 1] - 8 * inData1[j];
			if (max < outDate[j]){
				max = outDate[j];
			}
			if (min > outDate[j]){
				min = outDate[j];
			}
		}
	}
	double size = (max - min) / 256.0;
	for (int i = 1; i < nr - 1; i++){
		int* inDate = matrix[i];
		uchar* outDate = LaplaceImage.ptr<uchar>(i);
		for (int j = 1; j < nc - 1; j++){
			outDate[j] = (inDate[j] - min) / size;
		}
	}
	//	return lab14(LaplaceImage);
	return LaplaceImage;
}
int main(int argc, char** argv){
	String s;
	if (argc > 1) {
		s = argv[1];
	}
	else {
		s = "../source/blood1.bmp";
	}
	Mat image = imread(s, CV_8U);
	if (image.empty()) {
		cout << "Can't open image " << s << endl;
		return -1;
	}
	imshow("origin image", image);
	Mat RobertsImage = lab31(image);
	imshow("Roberts Image", RobertsImage);
	Mat SobelImage = lab33(image);
	imshow("Sobel Image", SobelImage);
	Mat PrewittImage = lab33(image);
	imshow("Prewitt Image", PrewittImage);
	Mat LaplaceImage1 = lab341(image);
	imshow("Laplace Image 1", LaplaceImage1);
	Mat LaplaceImage2 = lab342(image);
	imshow("Laplace Image 2", LaplaceImage2);
	waitKey();
	return 0;
}