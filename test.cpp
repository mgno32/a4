#include "CImg.h"
#include <iostream>
#include <string>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
using namespace cimg_library;
using namespace std;
typedef unsigned char uchar;
int LINENUM = 4;
float max_grad = 400;
const double PI = 3.1415926;
const double DIFF = 200;
int Threshold = 650;

void RGBtoGray(CImg<uchar> &image, CImg<uchar> &gray){
	cimg_forXY(image,x,y) {
		int R = (int)image(x,y,0,0);
		int G = (int)image(x,y,0,1);
		int B = (int)image(x,y,0,2);
		int grayWeight = (int)(0.299*R + 0.587*G + 0.114*B + 0.5);
		if (grayWeight > 255)
			grayWeight = 255;
		gray(x,y,0,0) = uchar(grayWeight);
	}
} 

float dis(float x, float y){
	return sqrt(x*x + y*y);
}

class pix{
	public:
		int x;
		int y;
		int value;
		pix(int x,int y,int value){
			this->x = x;
			this->y = y;
			this->value = value;
		}
}; 

class dot{ 
	public:
	double x;
	double y;
	dot(double a, double b){
		x = a;
		y = b;
	}
};

int P2X(double theta, double dist, int x){
	double ang = theta * PI / 180.0;
	double m = -cos(ang) / sin(ang);
	double b = dist / sin(ang);
	return m * x + b;
}

int P2Y(double theta, double dist, int y){
	double ang = theta * PI / 180.0;
	double b = dist / sin(ang);
	return (y - b) * (-tan(ang));
}


double get_distance(double dx, double dy){
	return sqrt(dx * dx + dy * dy);
}

vector<pix> findLines(CImg<float> &hough, CImg<uchar> image){
	vector<pix> maxP;
	const int ymin = 0;
	const int ymax = image.height() - 1;
	const int xmin = 0;
	const int xmax = image.width() - 1;
	cimg_forXY(hough,angle,polar){
		// 判断是否为峰值
		float v = hough(angle, polar); // 投票数
		if(v > Threshold){
			// 图像范围
			int x0 = P2Y(angle, polar, ymin);
			int x1 = P2Y(angle, polar, ymax);
			int y0 = P2X(angle, polar, xmin);
			int y1 = P2X(angle, polar, xmax);
			//cout << x0 << "!"<< x1 << "|" << xmax << endl;
			//cout << y0 << "?"<< y1 << "}" << ymax << endl;
			if ((x0 < 0 || x0 > xmax) && \
			(x1 < 0 || x1 > xmax) && \
			(y0 < 0 || y0 > ymax) && \
			(y1 < 0 || y1 > ymax)) continue;
			bool flag = false;
			for (pix &p : maxP){
				if (get_distance(p.x - angle, p.y - polar) < DIFF){
					flag = true;
					if (p.value < v){
						p = pix(angle, polar, v);
					}
				}
			}
			if (!flag){
				maxP.push_back(pix(angle, polar, v));
			}
 		}
	}

	for (pix p: maxP){
		cout << "+=" << p.x << "," << p.y<<endl;
	}
	
 	return maxP;
}

void printOutLineAndDots(vector<pix> hough_line,CImg<uchar> image){
	double k,b;
	vector<dot> kbLine;
	vector<dot> intersec;
	for(size_t i = 0; i < hough_line.size(); i++){
	    double angle  = hough_line[i].x * PI / 180;	
		k = -cos(angle)/sin(angle);
		b = hough_line[i].y / sin(angle);
		kbLine.push_back(dot(k,b));
		cout<<"y = "<<k<<" x + "<<b<<endl;
	}
	for(size_t i = 0; i < kbLine.size(); i++){
		for(size_t j = i + 1; j < kbLine.size(); j++){
			double m0 = kbLine[i].x;
			double b0 = kbLine[i].y;
			double m1 = kbLine[j].x;
			double b1 = kbLine[j].y;
			double x1 = (b1-b0) / (m0-m1);
			double y1 = (m0 * b1 - m1 * b0)/(m0 - m1);
			cout << m0 << ", " << b0 << endl;
			cout << m1 << ",**" << b1 << endl;
			if( x1 >= 0 && y1 >= 0 && x1 < image.width() && y1 < image.height()){
				intersec.push_back(dot(x1,y1)); // 添加交点
				cout << "DOT:(" << x1 << ", " << y1 << ")" << endl; 
			}
		}	
	}
	for(size_t i = 0; i < intersec.size();i++){
		int dot_x = intersec[i].x;
		int dot_y = intersec[i].y;
		uchar color[] = { 100,100,100};
		image.draw_circle(dot_x,dot_y,50,color);

	}
	image.display();
}

int main() {
	CImg<uchar> image("1.jpg");	
	CImg<uchar> grayImage(image.width(),image.height());
	CImg<uchar> blurImage(image.width(),image.height());
	CImg<float> gbImage(image.width(), image.height());
	
	RGBtoGray(image,grayImage);
	gbImage = grayImage.get_blur_median(3);
	gbImage.blur(3);
	gbImage.display();
	
	CImg<float> gradnum(image.width(), image.height());
	CImg_3x3(I,float);
	CImg<float> hough_vote(360,dis(image.width(),image.height()),1,1,0);
	cimg_for3x3(gbImage,x,y,0,0,I,float){
		double ix = Inc - Ipc;
		double iy = Icp - Icn;
		double graditude = ix * ix + iy * iy;
		if (graditude > max_grad){
			gradnum(x,y) = graditude;
			cimg_forX(hough_vote,angle){
				double r = (double)angle * PI / 180;
				int polar = (int)(x * cos(r) + y * sin(r));
				if(polar >= 0 && polar < hough_vote.height()){
					hough_vote(angle,polar) += 1;
				}
			}
		}
	}
	gradnum.display();
	hough_vote.display();
	vector<pix> hough_line = findLines(hough_vote, image);
	cout << hough_line.size() << endl;
	printOutLineAndDots(hough_line,image);
	return 0;
}
