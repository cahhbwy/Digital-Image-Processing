//《计算机图形与图像》---数字图像处理大作业
//图像配准
//手动选取图中点
//线性变换进行拼接

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace std;
using namespace cv;

Mat image1, image2;	//待配准图像
Mat selectImage;
vector<Point> points1, points2;
string result_name = "result.jpg";

//输入图像
bool inputImage() {
	cout << "输入两幅图像文件名" << endl;
	string source_name;
	cin >> source_name;
	image1 = imread(source_name);
	if (image1.empty()) {
		cout << "Can't read image '" << source_name << endl << "Please try again" << endl;
		return false;
	}
	cin >> source_name;
	image2 = imread(source_name);
	if (image2.empty()) {
		cout << "Can't read image '" << source_name << endl << "Please try again" << endl;
		return false;
	}
	return true;
}
//手动选点
void mouseHandler(int event, int x, int y, int flags, void* param) {
	switch (event) {
		case CV_EVENT_LBUTTONDOWN:
			if (x < image1.cols) {
				points1.push_back(Point(y + image1.rows, x + image1.cols));
			} else if (x > image1.cols + 10) {
				points2.push_back(Point(y, x - image1.cols - 10));
			}
			break;
		case CV_EVENT_RBUTTONDOWN:
			if (x < image1.cols) {
				if (!points1.empty())
					points1.pop_back();
			} else if (x > image1.cols + 10) {
				if (!points2.empty())
					points2.pop_back();
			}
			break;
		default:break;
	}
	Mat temp = selectImage.clone();
	for (unsigned int i = 0; i < points1.size(); ++i) {
		circle(temp, Point(points1[i].y - image1.cols, points1[i].x - image1.rows), 3, Scalar(255, 0, 0), 2);
	}
	for (unsigned int i = 0; i < points2.size(); ++i) {
		circle(temp, Point(points2[i].y + image1.cols + 10, points2[i].x), 3, Scalar(255, 0, 0), 2);
	}
	for (int i = 0; i < points1.size() && i < points2.size(); ++i) {
		line(temp, Point(points1[i].y - image1.cols, points1[i].x - image1.rows), Point(points2[i].y + image1.cols + 10, points2[i].x), Scalar(255, 0, 0));
	}
	imshow("select", temp);
}
int selectPoint() {
	cout << "please select points order by up to down and left to right" << endl;
	selectImage = Mat(max(image1.rows, image2.rows), image1.cols + image2.cols + 10, image1.type());
	Mat temp1(selectImage, Rect(0, 0, image1.cols, image1.rows));
	Mat temp2(selectImage, Rect(image1.cols + 10, 0, image2.cols, image2.rows));
	image1.copyTo(temp1);
	image2.copyTo(temp2);
	imshow("select", selectImage);
	setMouseCallback("select", mouseHandler, 0);
	waitKey();
	cvDestroyWindow("select");
	return min(points2.size(), points1.size());
}
//图像拼接
Mat cutBorder(Mat image) {
	int nr = image.rows;
	int nc = image.cols;
	int border[] = { 0, image.rows - 1, 0, image.cols - 1 };
	bool state;
	if (image.channels() == 3) {
		for (int i = 0; i < nr; ++i) {
			state = true;
			for (int j = 0; j < nc; ++j) {
				if (image.at<Vec3b>(i, j)[0] != 0 || image.at<Vec3b>(i, j)[1] != 0 || image.at<Vec3b>(i, j)[2] != 0)
					state = false;
			}
			if (state) {
				border[0] = i;
			} else {
				break;
			}
		}
		for (int i = nr - 1; i >= 0; --i) {
			state = true;
			for (int j = 0; j < nc; ++j) {
				if (image.at<Vec3b>(i, j)[0] != 0 || image.at<Vec3b>(i, j)[1] != 0 || image.at<Vec3b>(i, j)[2] != 0)
					state = false;
			}
			if (state) {
				border[1] = i;
			} else {
				break;
			}
		}
		for (int j = 0; j < nc; ++j) {
			state = true;
			for (int i = 0; i < nr; ++i) {
				if (image.at<Vec3b>(i, j)[0] != 0 || image.at<Vec3b>(i, j)[1] != 0 || image.at<Vec3b>(i, j)[2] != 0)
					state = false;
			}
			if (state) {
				border[2] = j;
			} else {
				break;
			}
		}
		for (int j = nc - 1; j >= 0; --j) {
			state = true;
			for (int i = 0; i < nr; ++i) {
				if (image.at<Vec3b>(i, j)[0] != 0 || image.at<Vec3b>(i, j)[1] != 0 || image.at<Vec3b>(i, j)[2] != 0)
					state = false;
			}
			if (state) {
				border[3] = j;
			} else {
				break;
			}
		}
	} else if (image.channels() == 1) {
		for (int i = 0; i < nr; ++i) {
			state = true;
			for (int j = 0; j < nc; ++j) {
				if (image.at<uchar>(i, j) != 205)
					state = false;
			}
			if (state) {
				border[0] = i;
			} else {
				break;
			}
		}
		for (int i = nr - 1; i >= 0; --i) {
			state = true;
			for (int j = 0; j < nc; ++j) {
				if (image.at<uchar>(i, j) != 205)
					state = false;
			}
			if (state) {
				border[1] = i;
			} else {
				break;
			}
		}
		for (int j = 0; j < nc; ++j) {
			state = true;
			for (int i = 0; i < nr; ++i) {
				if (image.at<uchar>(i, j) != 205)
					state = false;
			}
			if (state) {
				border[2] = j;
			} else {
				break;
			}
		}
		for (int j = nc - 1; j >= 0; --j) {
			state = true;
			for (int i = 0; i < nr; ++i) {
				if (image.at<uchar>(i, j) != 205)
					state = false;
			}
			if (state) {
				border[3] = j;
			} else {
				break;
			}
		}
	}
	return Mat(image, Rect(border[2], border[0], border[3] - border[2] + 1, border[1] - border[0] + 1));
}
Mat stitching1() {
	Mat outImage = Mat(image1.rows * 3, image1.cols * 3, image1.type());
	Mat temp(outImage, Rect(image1.cols, image1.rows, image1.cols, image1.rows));
	image1.copyTo(temp);
	//计算变换矩阵，相当于两次求解矛盾方程组
	double trans[2][3];
	//变换矩阵第一行系数计算
	double ATA[3][3], ATy[3];	//矛盾方程组法方程
	int n = min(points1.size(), points2.size());
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			ATA[i][j] = 0.0;
		}
		ATy[i] = 0.0;
	}
	for (int i = 0; i < n; ++i) {
		ATA[0][1] += points1[i].x;
		ATA[0][2] += points1[i].y;
		ATA[1][1] += points1[i].x * points1[i].x;
		ATA[1][2] += points1[i].x * points1[i].y;
		ATA[2][2] += points1[i].y * points1[i].y;
		ATy[0] += points2[i].x;
		ATy[1] += points1[i].x * points2[i].x;
		ATy[2] += points1[i].y * points2[i].x;
	}
	ATA[0][0] = n;
	ATA[1][0] = ATA[0][1];
	ATA[2][0] = ATA[0][2];
	ATA[2][1] = ATA[1][2];
	double AF[9];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			AF[i * 3 + j] = ATA[i][j];
		}
	}
	CvMat *My = cvCreateMat(3, 1, CV_64FC1);
	cvSetData(My, ATy, CV_AUTOSTEP);
	CvMat *MA = cvCreateMat(3, 3, CV_64FC1);
	cvSetData(MA, AF, CV_AUTOSTEP);
	CvMat *Mx = cvCreateMat(3, 1, CV_64FC1);
	cvSolve(MA, My, Mx, CV_LU);
	for (int i = 0; i < 3; ++i) {
		trans[0][i] = Mx->data.db[i];
	}
	//变换矩阵第二行系数计算
	for (int i = 0; i < 3; ++i) {
		ATy[i] = 0.0;
	}
	for (int i = 0; i < n; ++i) {
		ATy[0] += points2[i].y;
		ATy[1] += points1[i].x * points2[i].y;
		ATy[2] += points1[i].y * points2[i].y;
	}
	cvSetData(My, ATy, CV_AUTOSTEP);
	cvSetData(MA, AF, CV_AUTOSTEP);
	cvSolve(MA, My, Mx, CV_LU);
	for (int i = 0; i < 3; ++i) {
		trans[1][i] = Mx->data.db[i];
	}
	// 	trans[0][0] = sqrt(2) / 2;
	// 	trans[0][1] = sqrt(2) / 2;
	// 	trans[0][2] = -image1.rows;
	// 	trans[1][0] = -sqrt(2) / 2;
	// 	trans[1][1] = sqrt(2) / 2;
	// 	trans[1][2] = -image1.cols*2;

	cout << trans[0][0] << " " << trans[0][1] << " " << trans[0][2] << endl;
	cout << trans[1][0] << " " << trans[1][1] << " " << trans[1][2] << endl;
	//进行矩阵变换
	int nr = outImage.rows, nc = outImage.cols;
	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			double x = trans[0][0] + trans[0][1] * i + trans[0][2] * j;
			double y = trans[1][0] + trans[1][1] * i + trans[1][2] * j;
			if (x >= 0.5 && x < image2.rows - 1.5 && y >= 0.5 && y < image2.cols - 1.5) {
				if (outImage.channels() == 3) {
					double r = (image2.at<Vec3b>((int)x, (int)y)[0] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y)[0] * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<Vec3b>((int)x, (int)y + 1)[0] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y + 1)[0] * (x - (int)x))*(y - (int)y);
					double g = (image2.at<Vec3b>((int)x, (int)y)[1] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y)[1] * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<Vec3b>((int)x, (int)y + 1)[1] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y + 1)[1] * (x - (int)x))*(y - (int)y);
					double b = (image2.at<Vec3b>((int)x, (int)y)[2] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y)[2] * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<Vec3b>((int)x, (int)y + 1)[2] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y + 1)[2] * (x - (int)x))*(y - (int)y);
					if (i >= image1.rows&&i < image1.rows * 2 && j >= image1.cols&&j < image1.cols * 2) {
						outImage.at<Vec3b>(i, j)[0] = (outImage.at<Vec3b>(i, j)[0] + r) / 2;
						outImage.at<Vec3b>(i, j)[1] = (outImage.at<Vec3b>(i, j)[1] + g) / 2;
						outImage.at<Vec3b>(i, j)[2] = (outImage.at<Vec3b>(i, j)[2] + b) / 2;
					}
					outImage.at<Vec3b>(i, j)[0] = r;
					outImage.at<Vec3b>(i, j)[1] = g;
					outImage.at<Vec3b>(i, j)[2] = b;
				} else if (outImage.channels() == 1) {
					double gray = (image2.at<uchar>((int)x, (int)y) * ((int)x + 1 - x) + image2.at<uchar>((int)x + 1, (int)y) * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<uchar>((int)x, (int)y + 1) * ((int)x + 1 - x) + image2.at<uchar>((int)x + 1, (int)y + 1) * (x - (int)x))*(y - (int)y);
					if (i >= image1.rows&&i < image1.rows * 2 && j >= image1.cols&&j < image1.cols * 2) {
						outImage.at<uchar>(i, j) = (outImage.at<uchar>(i, j) + gray) / 2;
					}
					outImage.at<uchar>(i, j) = gray;
				}
			}
		}
	}
	return cutBorder(outImage);
}
Mat stitching2() {
	Mat outImage = Mat(image1.rows * 3, image1.cols * 3, image1.type());
	Mat temp(outImage, Rect(image1.cols, image1.rows, image1.cols, image1.rows));
	image1.copyTo(temp);
	//计算变换矩阵，相当于两次求解矛盾方程组
	//u=a0+a1*x+a2*y+a3*x*x+a4*y*y+a5*x*y
	//v=b0+b1*x+b2*y+b3*x*x+b4*y*y+b5*x*y
	double trans[2][6];
	//变换矩阵第一行系数计算
	double ATA[6][6], ATy[6];	//矛盾方程组法方程
	int n = min(points1.size(), points2.size());
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			ATA[i][j] = 0.0;
		}
		ATy[i] = 0.0;
	}
	double x, y;
	for (int i = 0; i < n; ++i) {
		x = points1[i].x;
		y = points1[i].y;
		ATA[0][1] += x;
		ATA[0][2] += y;
		ATA[0][3] += x*x;
		ATA[0][4] += y*y;
		ATA[0][5] += x*y;
		ATA[1][3] += x*x*x;
		ATA[1][4] += x*y*y;
		ATA[1][5] += x*x*y;
		ATA[2][4] += y*y*y;
		ATA[3][3] += x*x*x*x;
		ATA[3][4] += x*x*y*y;
		ATA[3][5] += x*x*x*y;
		ATA[4][4] += y*y*y*y;
		ATA[4][5] += x*y*y*y;
		ATy[0] += points2[i].x;
		ATy[1] += x * points2[i].x;
		ATy[2] += y * points2[i].x;
		ATy[3] += x * x * points2[i].x;
		ATy[4] += y * y * points2[i].x;
		ATy[5] += x * y * points2[i].x;
	}
	ATA[0][0] = n;
	ATA[1][0] = ATA[0][1];
	ATA[1][1] = ATA[0][3];
	ATA[1][2] = ATA[0][5];
	ATA[2][0] = ATA[0][2];
	ATA[2][1] = ATA[1][2];
	ATA[2][2] = ATA[0][4];
	ATA[2][3] = ATA[1][5];
	ATA[2][5] = ATA[1][4];
	ATA[3][0] = ATA[0][3];
	ATA[3][1] = ATA[1][3];
	ATA[3][2] = ATA[2][3];
	ATA[4][0] = ATA[0][4];
	ATA[4][1] = ATA[1][4];
	ATA[4][2] = ATA[2][4];
	ATA[4][3] = ATA[3][4];
	ATA[5][0] = ATA[0][5];
	ATA[5][1] = ATA[1][5];
	ATA[5][2] = ATA[2][5];
	ATA[5][3] = ATA[3][5];
	ATA[5][4] = ATA[4][5];
	ATA[5][5] = ATA[3][4];

	double AF[36];
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			AF[i * 6 + j] = ATA[i][j];
		}
	}
	CvMat *My = cvCreateMat(6, 1, CV_64FC1);
	cvSetData(My, ATy, CV_AUTOSTEP);
	CvMat *MA = cvCreateMat(6, 6, CV_64FC1);
	cvSetData(MA, AF, CV_AUTOSTEP);
	CvMat *Mx = cvCreateMat(6, 1, CV_64FC1);
	cvSolve(MA, My, Mx, CV_LU);
	for (int i = 0; i < 6; ++i) {
		trans[0][i] = Mx->data.db[i];
	}
	//变换矩阵第二行系数计算
	for (int i = 0; i < 6; ++i) {
		ATy[i] = 0.0;
	}
	for (int i = 0; i < n; ++i) {
		x = points1[i].x;
		y = points1[i].y;
		ATy[0] += points2[i].y;
		ATy[1] += x * points2[i].y;
		ATy[2] += y * points2[i].y;
		ATy[3] += x * x * points2[i].y;
		ATy[4] += y * y * points2[i].y;
		ATy[5] += x * y * points2[i].y;
	}
	cvSetData(My, ATy, CV_AUTOSTEP);
	cvSetData(MA, AF, CV_AUTOSTEP);
	cvSolve(MA, My, Mx, CV_LU);
	for (int i = 0; i < 6; ++i) {
		trans[1][i] = Mx->data.db[i];
	}
	cout << trans[0][0] << " " << trans[0][1] << " " << trans[0][2] << " " << trans[0][3] << " " << trans[0][4] << " " << trans[0][5] << " " << endl;
	cout << trans[1][0] << " " << trans[1][1] << " " << trans[1][2] << " " << trans[1][3] << " " << trans[1][4] << " " << trans[1][5] << " " << endl;
	//进行矩阵变换
	int nr = outImage.rows, nc = outImage.cols;
	for (int i = 0; i < nr; ++i) {
		for (int j = 0; j < nc; ++j) {
			double x = trans[0][0] + trans[0][1] * i + trans[0][2] * j + trans[0][3] * i*i + trans[0][4] * j*j + trans[0][5] * i*j;
			double y = trans[1][0] + trans[1][1] * i + trans[1][2] * j + trans[1][3] * i*i + trans[1][4] * j*j + trans[1][5] * i*j;
			if (x >= 0.5 && x < image2.rows - 1.5 && y >= 0.5 && y < image2.cols - 1.5) {
				if (outImage.channels() == 3) {
					double r = (image2.at<Vec3b>((int)x, (int)y)[0] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y)[0] * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<Vec3b>((int)x, (int)y + 1)[0] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y + 1)[0] * (x - (int)x))*(y - (int)y);
					double g = (image2.at<Vec3b>((int)x, (int)y)[1] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y)[1] * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<Vec3b>((int)x, (int)y + 1)[1] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y + 1)[1] * (x - (int)x))*(y - (int)y);
					double b = (image2.at<Vec3b>((int)x, (int)y)[2] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y)[2] * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<Vec3b>((int)x, (int)y + 1)[2] * ((int)x + 1 - x) + image2.at<Vec3b>((int)x + 1, (int)y + 1)[2] * (x - (int)x))*(y - (int)y);
					if (i >= image1.rows&&i < image1.rows * 2 && j >= image1.cols&&j < image1.cols * 2) {
						outImage.at<Vec3b>(i, j)[0] = (outImage.at<Vec3b>(i, j)[0] + r) / 2;
						outImage.at<Vec3b>(i, j)[1] = (outImage.at<Vec3b>(i, j)[1] + g) / 2;
						outImage.at<Vec3b>(i, j)[2] = (outImage.at<Vec3b>(i, j)[2] + b) / 2;
					} else {
						outImage.at<Vec3b>(i, j)[0] = r;
						outImage.at<Vec3b>(i, j)[1] = g;
						outImage.at<Vec3b>(i, j)[2] = b;
					}
				} else if (outImage.channels() == 1) {
					double gray = (image2.at<uchar>((int)x, (int)y) * ((int)x + 1 - x) + image2.at<uchar>((int)x + 1, (int)y) * (x - (int)x))*((int)y + 1 - y)
						+ (image2.at<uchar>((int)x, (int)y + 1) * ((int)x + 1 - x) + image2.at<uchar>((int)x + 1, (int)y + 1) * (x - (int)x))*(y - (int)y);
					if (i >= image1.rows&&i < image1.rows * 2 && j >= image1.cols&&j < image1.cols * 2) {
						outImage.at<uchar>(i, j) = (outImage.at<uchar>(i, j) + gray) / 2;
					} else {
						outImage.at<uchar>(i, j) = gray;
					}
				}
			}
		}
	}
	return cutBorder(outImage);
}
int main(int argc, char** argv) {
	if (!inputImage()) {
		return -1;
	}
	int select;
	cout << "选择变换方法：1.一阶变换 2.二阶变换" << endl;
	cin >> select;
	if (select == 1) {
		if (selectPoint() < 4) {
			cout << "too few points" << endl;
			return 2;
		}
		for (int i = 0; i < points1.size(); ++i) {
			cout << points1[i] << " " << points2[i] << endl;
		}
		Mat result = stitching1();
		imwrite("result.jpg", result);
	} else if (select == 2) {
		if (selectPoint() < 7) {
			cout << "too few points" << endl;
			return 2;
		}
		for (int i = 0; i < points1.size(); ++i) {
			cout << points1[i] << " " << points2[i] << endl;
		}
		Mat result = stitching2();
		imwrite("result.jpg", result);
	}
	return 0;
}