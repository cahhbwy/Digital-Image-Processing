#include <iostream>
#include <ctime>
#include <complex>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#define M_PI       3.14159265358979323846
using namespace std;
using namespace cv;
vector<vector<complex<double> > >fourier(Mat image) {
	int nr = image.rows;
	int nc = image.cols;
	complex<double >j(0, 1);
	complex<double >nj(0, -1);
	vector<vector<complex<double> > >outMatrix;
	outMatrix.resize(nr);
	for (int i = 0; i<nr; ++i) {
		outMatrix[i].resize(nc);
	}
	for (int u = 0; u<nr; ++u) {
		for (int v = 0; v<nc; ++v) {
			complex<double >t;
			for (int x = 0; x<nr; ++x) {
				const uchar *inData = image.ptr<uchar>(x);
				for (int y = 0; y<nc; ++y) {
					t = t +(double)inData[y] * exp(nj *(2 * M_PI *((double)u* x / nr +(double)v* y /nc)));
				}
			}
//                      t = t / (double)(nr*nc);
			outMatrix[u][v] = t;
		}
	}
	return outMatrix;
}
vector<vector<complex<double> > > ifourier(vector<vector<complex<double> > >inMatrix) {
	int nr = inMatrix.size();
	int nc = inMatrix[0].size();
	complex<double >j(0, 1);
	complex<double >nj(0, -1);
	vector<vector<complex<double> > >outMatrix;
	outMatrix.resize(nr);
	for (int i = 0; i<nr; ++i) {
		outMatrix[i].resize(nc);
	}
	for (int x = 0; x<nr; ++x) {
		for (int y = 0; y<nc; ++y) {
			complex<double >t;
			for (int u = 0; u<nr; ++u) {
				for (int v = 0; v<nc; ++v) {
					t = t + complex<double>(inMatrix[u][v].real(),inMatrix[u][v].imag()) *exp(j *(2 * M_PI *((double)u * x / nr +(double)v * y / nc)));
				}
			}
			//t = t / (double)(nr*nc);
			outMatrix[x][y] = t;
		}
	}
	return outMatrix;
}
Mat lab14(Mat image) {
	int L = 256;
	unsigned long *hist =(unsigned long *)calloc(L, sizeof(unsigned long));
	int nr = image.rows;
	int nc = image.cols;
	for (int k = 0; k<nr; k++) {
		const uchar *inData = image.ptr<uchar>(k);
		for (int i = 0; i<nc; i++) {
			hist[inData[i]]++;
		}
	}
	unsigned long sum = nr * nc;
	int temp = 0;
	int *f = (int *)calloc(L, sizeof(int));
	for (int i = 0; i<L; i++) {
		temp += hist[i];
		f[i] = (int)((L - 1) * (double)temp / sum + 0.5) + (256 / L) / 2;
	}
	Mat newImage(image.size(), image.type(), Scalar(255));
	for (int k = 0; k<nr; k++) {
		const uchar *inData = image.ptr<uchar>(k);
		uchar * outData = newImage.ptr<uchar>(k);
		for (int i = 0; i<nc; i++) {
			outData[i] = f[inData[i] / (256 / L)];
		}
	}
	return newImage;
}
vector<complex<double> >fft(vector<complex<double> >x) {
	int n = x.size();
	if (n == 1)
		return x;
	complex<double >j(0, 1);
	complex<double >nj(0, -1);
	complex<double >wn = exp(nj * (2 * M_PI / n));
	complex<double >w(1, 0);
	vector<complex<double >>x0, x1;
	for (int i = 0; i<n;) {
		x0.push_back(x[i++]);
		x1.push_back(x[i++]);
	}
	vector<complex<double >>y0 = fft(x0);
	vector<complex<double >>y1 = fft(x1);
	vector<complex<double >>y;
	y.resize(n);
	for (int k = 0; k<n / 2; ++k) {
		y[k] = y0[k] + w * y1[k];
		y[k + n / 2] = y0[k] - w * y1[k];
		w = w * wn;
	}
	return y;
}
vector<vector<complex<double> > >fastFourier(Mat image) {
	int nr = image.rows;
	int nc = image.cols;
	vector<vector<complex<double> > >outMatrix;
	outMatrix.resize(nr);
	vector<complex<double> >tempCol, tempRow, tempCol2;
	tempCol.resize(nr);
	tempRow.resize(nc);
	int t;
	for (int i = 0; i < nr; ++i) {
		const uchar *inData = image.ptr<uchar>(i);
		for (int j = 0; j<nc; ++j) {
			t = (int)inData[j];
			tempRow[j] = complex<double>(t);
		}
		outMatrix[i] = fft(tempRow);
	}
	for (int j = 0; j<nc; ++j) {
		for (int i = 0; i<nr; ++i) {
			tempCol[i] = outMatrix[i][j];
		}
		tempCol2 = fft(tempCol);
		for (int i = 0; i<nr; i++) {
			outMatrix[i][j] = tempCol2[i];
		}
	}
	return outMatrix;
}
vector<vector<complex<double> > >ifft(vector<vector<complex<double> > >inMatrix) {
	int nr = inMatrix.size();
	int nc = inMatrix[0].size();
	vector<vector<complex<double> > >outMatrix;
	outMatrix.resize(nr);
	vector<complex<double >>tempCol, tempRow, tempCol2;
	tempCol.resize(nr);
	tempRow.resize(nc);
	for (int i = 0; i<nr; ++i) {
		for (int j = 0; j<nc; ++j) {
			tempRow[j] =complex<double>(inMatrix[i][j].real(), -inMatrix[i][j].imag());
		}
		outMatrix[i] = fft(tempRow);
	}
	for (int j = 0; j<nc; ++j) {
		for (int i = 0; i<nr; ++i) {
			tempCol[i] = outMatrix[i][j];
		}
		tempCol2 = fft(tempCol);
		for (int i = 0; i<nr; i++) {
			outMatrix[i][j] =complex<double >(tempCol2[i].real(), -tempCol2[i].imag());
		}
	}
	return outMatrix;
}
vector<vector<complex<double> > > filter(vector<vector<complex<double> > >inMatrix, double r0,double r1, bool pass) {
	vector<vector<complex<double> > >outMatrix;
	int nr = inMatrix.size();
	int nc = inMatrix[0].size();
	outMatrix.resize(nr);
	for (int i = 0; i<nr; ++i) {
		outMatrix[i].resize(nc);
	}
	double dis;
	for (int i = 0; i<nr / 2; ++i) {
		for (int j = 0; j<nc / 2; ++j) {
			dis = sqrt(i * i + j * j);
			if (dis >= r0 && dis<r1) {
				outMatrix[i][j] =pass ? inMatrix[i][j] : complex<double>(0, 0);
			} else {
				outMatrix[i][j] =pass ? complex<double>(0,0) :inMatrix[i][j];
			}
		}
	}
	for (int i = nr / 2; i<nr; ++i) {
		for (int j = 0; j<nc / 2; ++j) {
			dis = sqrt((i - nr + 1) * (i - nr + 1) + j * j);
			if (dis >= r0 && dis<r1) {
				outMatrix[i][j] =pass ? inMatrix[i][j] : complex<double>(0, 0);
			} else {
				outMatrix[i][j] =pass ? complex<double>(0,0) : inMatrix[i][j];
			}
		}
	}
	for (int i = 0; i<nr / 2; ++i) {
		for (int j = nc / 2; j<nc; ++j) {
			dis = sqrt(i * i + (j - nc + 1) * (j - nc + 1));
			if (dis >= r0 && dis<r1) {
				outMatrix[i][j] =pass ? inMatrix[i][j] : complex<double>(0, 0);
			} else {
				outMatrix[i][j] =pass ? complex<double>(0,0) :inMatrix[i][j];
			}
		}
	}
	for (int i = nr / 2; i<nr; ++i) {
		for (int j = nc / 2; j<nc; ++j) {
			dis =sqrt((i - nr + 1) * (i - nr + 1) +(j - nc + 1) * (j - nc + 1));
			if (dis >= r0 && dis<r1) {
				outMatrix[i][j] =pass ? inMatrix[i][j] : complex <double >(0, 0);
			} else {
				outMatrix[i][j] =pass ? complex<double >(0,0) :inMatrix[i][j];
			}
		}
	}
	return outMatrix;
}
vector<vector<complex<double> > > ifftMargin(vector<vector<complex<double > > > inMatrix) {
	int nr = inMatrix.size();
	int nc = inMatrix[0].size();
	vector<vector<complex<double> > >outMatrix;
	outMatrix.resize(nr);
	vector<complex<double> >tempCol, tempRow, tempCol2;
	tempCol.resize(nr);
	tempRow.resize(nc);
	for (int i = 0; i<nr; ++i) {
		for (int j = 0; j<nc; ++j) {
			tempRow[j] = complex<double>(abs(inMatrix[i][j]));
		}
		outMatrix[i] = fft(tempRow);
	}
	for (int j = 0; j<nc; ++j) {
		for (int i = 0; i<nr; ++i) {
			tempCol[i] = outMatrix[i][j];
		}
		tempCol2 = fft(tempCol);
		for (int i = 0; i<nc; i++) {
			outMatrix[i][j] = tempCol2[i];
		}
	}
	return outMatrix;
}
vector<vector<complex<double> > >ifftPhase(vector<vector<complex<double> > >inMatrix) {
	int nr = inMatrix.size();
	int nc = inMatrix[0].size();
	vector<vector<complex<double> > >outMatrix;
	outMatrix.resize(nr);
	vector<complex<double> >tempCol, tempRow, tempCol2;
	tempCol.resize(nr);
	tempRow.resize(nc);
	for (int i = 0; i<nr; ++i) {
		for (int j = 0; j<nc; ++j) {
			tempRow[j] =complex<double>(inMatrix[i][j].imag(),inMatrix[i][j].real()) /(abs(inMatrix[i][j]) + 0.00001);
		}
		outMatrix[i] = fft(tempRow);
	}
	for (int j = 0; j<nc; ++j) {
		for (int i = 0; i<nr; ++i) {
			tempCol[i] = outMatrix[i][j];
		}
		tempCol2 = fft(tempCol);
		for (int i = 0; i<nc; i++) {
			outMatrix[i][j] = tempCol2[i];
		}
	}
	return outMatrix;
}
Mat showMargin(vector<vector<complex<double> > >matrix) {
	int nr = matrix.size();
	int nc = matrix[0].size();
	Mat outImage(nr, nc, CV_8U);
	double max = 0.0;
	for (int i = 0; i<nr; ++i) {
		for (int j = 0; j<nc; ++j) {
			max =max >log(abs(matrix[i][j]) +1.0) ? max : log(abs(matrix[i][j]) + 1.0);
		}
	}
	max = 255.0 / max;
	for (int i = 0; i<nr; ++i) {
		uchar * outData = outImage.ptr<uchar>(i);
		for (int j = 0; j<nc; ++j) {
			outData[j] =(uchar)(log(abs(matrix[i][j]) + 1.0) * max);
		}
	}
	int cx = outImage.cols / 2;
	int cy = outImage.rows / 2;
	Mat q0(outImage, Rect(0, 0, cx, cy));	
	Mat q1(outImage, Rect(cx, 0, cx, cy));	
	Mat q2(outImage, Rect(0, cy, cx, cy));	
	Mat q3(outImage, Rect(cx, cy, cx, cy));	
	Mat tmp;
	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);
	q1.copyTo(tmp);
	q2.copyTo(q1);
	tmp.copyTo(q2);
	return outImage;
}
Mat showPhase(vector<vector<complex<double > > >matrix) {
	int nr = matrix.size();
	int nc = matrix[0].size();
	Mat outImage(nr, nc, CV_8U);
	double max = 0.0;
	for (int i = 0; i<nr; ++i) {
		for (int j = 0; j<nc; ++j) {
			max = max >atan(matrix[i][j].imag() /matrix[i][j].real())? max : atan(matrix[i][j].imag() /matrix[i][j].real());
		}
	}
	max = 255.0 / max;
	for (int i = 0; i<nr; ++i) {
		uchar * outData = outImage.ptr<uchar>(i);
		for (int j = 0; j<nc; ++j) {
			outData[j] = (uchar) (atan (matrix[i][j].imag() / matrix[i][j].real()) * max);
		}
	}
	return outImage;
}
Mat showInverse(vector<vector<complex<double> > > matrix) {
	int nr = matrix.size();
	int nc = matrix[0].size();
	Mat outImage(nr, nc, CV_8U);
	double max = 0.0, min = 1e10;
	for (int i = 0; i<nr; ++i) {
		for (int j = 0; j<nc; ++j) {
			max =max > abs(matrix[i][j]) ? max : abs(matrix[i][j]);
			min =min < abs(matrix[i][j]) ? min : abs(matrix[i][j]);
		}
	}
	double scale = (max - min) / 255.0;
	for (int i = 0; i<nr; ++i) {
		uchar * outData = outImage.ptr<uchar>(i);
		for (int j = 0; j<nc; ++j) {
			outData[j] =(uchar) ((abs(matrix[i][j]) - min) / scale);
		}
	}
	return outImage;
}
int main(int argc, char **argv) {
	String s;
	if (argc > 1) {
		s = argv[1];
	} else {
		s = "../source/rect1.bmp";
	}
	Mat image = imread(s, CV_8U);
	if (image.empty()) {
		cout << "Can't open image " << s << endl;
		return -1;
	}
	imshow("image", image);
// 	vector<vector<complex<double> > > dftMatrix = fourier(image);
// 	Mat fourierImage = showMargin(dftMatrix);
// 	imshow("fourier image", fourierImage);
// 	vector<vector<complex<double> > > idftMatrix = ifourier(dftMatrix);
// 	Mat ifourierImage = showInverse(idftMatrix);
// 	imshow("inverse fourier image", ifourierImage);

	vector<vector<complex<double> > >fftMatrix = fastFourier(image);

	Mat fastfourierImageMargin = showMargin(fftMatrix);
//	fastfourierImage = lab14(fastfourierImage);
	imshow("fast fourier Margin image", fastfourierImageMargin);

	Mat fastfourierImagePhase = showPhase(fftMatrix);
//	fastfourierImage = lab14(fastfourierImage);
	imshow("fast fourier Phase image", fastfourierImagePhase);
	
	vector<vector<complex<double> > >ifftMatrixMargin =ifftMargin(fftMatrix);
	Mat ifastfourierImageMargin = showInverse(ifftMatrixMargin);
//	ifastfourierImageMargin = lab14(ifastfourierImageMargin);
	imshow("inverse fast fourier Margin image",ifastfourierImageMargin);
	
	vector<vector<complex<double> > >ifftMatrixPhase =ifftPhase(fftMatrix);
	Mat ifastfourierImagePhase = showInverse(ifftMatrixPhase);
//	ifastfourierImagePhase = lab14(ifastfourierImagePhase);
	imshow("inverse fast fourier Phase image", ifastfourierImagePhase);

	fftMatrix = filter(fftMatrix,0,25,true);
	
	vector<vector<complex<double> > > ifftMatrix = ifft(fftMatrix);
	Mat ifastfourierImage = showInverse(ifftMatrix);
	imshow("inverse fast fourier image", ifastfourierImage);

	
	waitKey();
	return 0;
}
