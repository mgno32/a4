#include "CImg.h"
#include <iostream>
#include <string>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include "interpolation.h"
using namespace cimg_library;
using namespace std;
typedef unsigned char uchar;
const int LINENUM = 4;
float max_grad = 400;
const double PI = 3.1415926;
double DIFF = 200;
int Threshold = 750;
const float PAPER_WIDTH = 1190.0;
const float PAPER_HEIGHT = 1684.0;
class dot;
vector<dot> intersec;
	
void RGBtoGray(CImg<uchar> &image, CImg<uchar> &gray){
	cimg_forXY(image,x, y) {
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

 	return maxP;
}

void printOutHough(vector<pix> houghLine, CImg<float> hough){
	dot houghMaxDot(0,0);
	uchar color[] = { 100,100,100};
		for (size_t i = 0; i < houghLine.size(); i++){
		houghMaxDot.x = houghLine[i].x;
		houghMaxDot.y = houghLine[i].y;
	
		hough.draw_circle(houghMaxDot.x,houghMaxDot.y,50,color);
	}
	hough.display();
}

void printOutLineAndDots(vector<pix>  hough_line,CImg<uchar> image){
	double k,b;
	vector<dot> kbLine;
	
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

	
			if( x1 >= 0 && y1 >= 0 && x1 < image.width() && y1 < image.height()){
				intersec.push_back(dot(x1,y1)); // 添加交点
			
			}
		}	
	}
	for(size_t i = 0; i < intersec.size();i++){
		int dot_x = intersec[i].x;
		int dot_y = intersec[i].y;
		uchar color[] = { 100,100,100};
		image.draw_circle(dot_x,dot_y,50,color);
 		cout<<"--------------------"<<endl;
		cout<<dot_x<<"  "<<dot_y<<endl;
	}
	image.display();
	image.save("DotDetect.jpg");
}

int get_thr(char* filename){
	int thr[6] = {500,750,650,800,300,650};
	int numOfTest = filename[0] - '1';
	return thr[numOfTest];
}
 
double get_diff(char* filename){
	int dif[6] = {50,200,200,200,200,100};
	int numOfTest = filename[0] - '1';
	return dif[numOfTest];
}

float get_grad(char* filename){
	float grad[6] = {400,400,400,200,400,400};
	int numOfTest = filename[0] - '1';
	return grad[numOfTest];
}

void swap(float & a, float& b){
	float tmp = a;
	a = b;
	b = tmp;
}

CImg<float> warp(vector<dot>& hough, const CImg<float> &src_img) {
	float u0 = 0, v0 = 0;
	float u1 = PAPER_WIDTH, v1 = 0;
	float u2 = 0, v2 = PAPER_HEIGHT;
	float u3 = PAPER_WIDTH, v3 = PAPER_HEIGHT;

	float x0 = hough[0].x, y0 = hough[0].y;
	float x1 = hough[1].x, y1 = hough[1].y;
	float x2 = hough[2].x, y2 = hough[2].y;
	float x3 = hough[3].x, y3 = hough[3].y;
	float c1, c2, c3, c4, c5, c6, c7, c8;

	c1 = -(u0*v0*v1*x2 - u0*v0*v2*x1 - u0*v0*v1*x3 + u0*v0*v3*x1 - u1*v0*v1*x2 + u1*v1*v2*x0 + u0*v0*v2*x3 - u0*v0*v3*x2 + u1*v0*v1*x3 - u1*v1*v3*x0 + u2*v0*v2*x1 - u2*v1*v2*x0
		- u1*v1*v2*x3 + u1*v1*v3*x2 - u2*v0*v2*x3 + u2*v2*v3*x0 - u3*v0*v3*x1 + u3*v1*v3*x0 + u2*v1*v2*x3 - u2*v2*v3*x1 + u3*v0*v3*x2 - u3*v2*v3*x0 - u3*v1*v3*x2 + u3*v2*v3*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c2 = (u0*u1*v0*x2 - u0*u2*v0*x1 - u0*u1*v0*x3 - u0*u1*v1*x2 + u0*u3*v0*x1 + u1*u2*v1*x0 + u0*u1*v1*x3 + u0*u2*v0*x3 + u0*u2*v2*x1 - u0*u3*v0*x2 - u1*u2*v2*x0 - u1*u3*v1*x0
		- u0*u2*v2*x3 - u0*u3*v3*x1 - u1*u2*v1*x3 + u1*u3*v1*x2 + u1*u3*v3*x0 + u2*u3*v2*x0 + u0*u3*v3*x2 + u1*u2*v2*x3 - u2*u3*v2*x1 - u2*u3*v3*x0 - u1*u3*v3*x2 + u2*u3*v3*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c3 = (u0*v1*x2 - u0*v2*x1 - u1*v0*x2 + u1*v2*x0 + u2*v0*x1 - u2*v1*x0 - u0*v1*x3 + u0*v3*x1 + u1*v0*x3 - u1*v3*x0 - u3*v0*x1 + u3*v1*x0
		+ u0*v2*x3 - u0*v3*x2 - u2*v0*x3 + u2*v3*x0 + u3*v0*x2 - u3*v2*x0 - u1*v2*x3 + u1*v3*x2 + u2*v1*x3 - u2*v3*x1 - u3*v1*x2 + u3*v2*x1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c4 = (u0*u1*v0*v2*x3 - u0*u1*v0*v3*x2 - u0*u2*v0*v1*x3 + u0*u2*v0*v3*x1 + u0*u3*v0*v1*x2 - u0*u3*v0*v2*x1 - u0*u1*v1*v2*x3 + u0*u1*v1*v3*x2 + u1*u2*v0*v1*x3 - u1*u2*v1*v3*x0 - u1*u3*v0*v1*x2 + u1*u3*v1*v2*x0
		+ u0*u2*v1*v2*x3 - u0*u2*v2*v3*x1 - u1*u2*v0*v2*x3 + u1*u2*v2*v3*x0 + u2*u3*v0*v2*x1 - u2*u3*v1*v2*x0 - u0*u3*v1*v3*x2 + u0*u3*v2*v3*x1 + u1*u3*v0*v3*x2 - u1*u3*v2*v3*x0 - u2*u3*v0*v3*x1 + u2*u3*v1*v3*x0)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c5 = -(u0*v0*v1*y2 - u0*v0*v2*y1 - u0*v0*v1*y3 + u0*v0*v3*y1 - u1*v0*v1*y2 + u1*v1*v2*y0 + u0*v0*v2*y3 - u0*v0*v3*y2 + u1*v0*v1*y3 - u1*v1*v3*y0 + u2*v0*v2*y1 - u2*v1*v2*y0
		- u1*v1*v2*y3 + u1*v1*v3*y2 - u2*v0*v2*y3 + u2*v2*v3*y0 - u3*v0*v3*y1 + u3*v1*v3*y0 + u2*v1*v2*y3 - u2*v2*v3*y1 + u3*v0*v3*y2 - u3*v2*v3*y0 - u3*v1*v3*y2 + u3*v2*v3*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c6 = (u0*u1*v0*y2 - u0*u2*v0*y1 - u0*u1*v0*y3 - u0*u1*v1*y2 + u0*u3*v0*y1 + u1*u2*v1*y0 + u0*u1*v1*y3 + u0*u2*v0*y3 + u0*u2*v2*y1 - u0*u3*v0*y2 - u1*u2*v2*y0 - u1*u3*v1*y0
		- u0*u2*v2*y3 - u0*u3*v3*y1 - u1*u2*v1*y3 + u1*u3*v1*y2 + u1*u3*v3*y0 + u2*u3*v2*y0 + u0*u3*v3*y2 + u1*u2*v2*y3 - u2*u3*v2*y1 - u2*u3*v3*y0 - u1*u3*v3*y2 + u2*u3*v3*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c7 = (u0*v1*y2 - u0*v2*y1 - u1*v0*y2 + u1*v2*y0 + u2*v0*y1 - u2*v1*y0 - u0*v1*y3 + u0*v3*y1 + u1*v0*y3 - u1*v3*y0 - u3*v0*y1 + u3*v1*y0
		+ u0*v2*y3 - u0*v3*y2 - u2*v0*y3 + u2*v3*y0 + u3*v0*y2 - u3*v2*y0 - u1*v2*y3 + u1*v3*y2 + u2*v1*y3 - u2*v3*y1 - u3*v1*y2 + u3*v2*y1)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	c8 = (u0*u1*v0*v2*y3 - u0*u1*v0*v3*y2 - u0*u2*v0*v1*y3 + u0*u2*v0*v3*y1 + u0*u3*v0*v1*y2 - u0*u3*v0*v2*y1 - u0*u1*v1*v2*y3 + u0*u1*v1*v3*y2 + u1*u2*v0*v1*y3 - u1*u2*v1*v3*y0 - u1*u3*v0*v1*y2 + u1*u3*v1*v2*y0
		+ u0*u2*v1*v2*y3 - u0*u2*v2*v3*y1 - u1*u2*v0*v2*y3 + u1*u2*v2*v3*y0 + u2*u3*v0*v2*y1 - u2*u3*v1*v2*y0 - u0*u3*v1*v3*y2 + u0*u3*v2*v3*y1 + u1*u3*v0*v3*y2 - u1*u3*v2*v3*y0 - u2*u3*v0*v3*y1 + u2*u3*v1*v3*y0)
		/ (u0*u1*v0*v2 - u0*u2*v0*v1 - u0*u1*v0*v3 - u0*u1*v1*v2 + u0*u3*v0*v1 + u1*u2*v0*v1 + u0*u1*v1*v3 + u0*u2*v0*v3 + u0*u2*v1*v2 - u0*u3*v0*v2 - u1*u2*v0*v2 - u1*u3*v0*v1
		- u0*u2*v2*v3 - u0*u3*v1*v3 - u1*u2*v1*v3 + u1*u3*v0*v3 + u1*u3*v1*v2 + u2*u3*v0*v2 + u0*u3*v2*v3 + u1*u2*v2*v3 - u2*u3*v0*v3 - u2*u3*v1*v2 - u1*u3*v2*v3 + u2*u3*v1*v3);

	CImg<float> res(PAPER_WIDTH, PAPER_HEIGHT, 1, src_img.spectrum(), 0);

	for (int cur_x = 0; cur_x < PAPER_WIDTH; cur_x++) {
		for (int cur_y = 0; cur_y < PAPER_HEIGHT; cur_y++) {
			double src_x = c1 * cur_x + c2 * cur_y + c3 * cur_x * cur_y + c4;
			double src_y = c5 * cur_x + c6 * cur_y + c7 * cur_x * cur_y + c8;

			double src_col = src_x, src_row = src_img.height() - src_y - 1;
			int cur_col = cur_x, cur_row = PAPER_HEIGHT - cur_y - 1;

			for (int channel = 0; channel < src_img.spectrum(); channel++) {
				res(cur_col, cur_row, channel) = bilinear_interpolation(src_img, src_col, src_row, channel);
			}

		}
	}
	
	return res;
}

int main(int agrc,char **argv) {	
	char *filename = argv[1];
	Threshold = get_thr(filename);
	DIFF = get_diff(filename);
	max_grad = get_grad(filename);
	
	CImg<uchar> image(filename);	
	CImg<uchar> grayImage(image.width(),image.height());
	CImg<uchar> blurImage(image.width(),image.height());
	CImg<float> gbImage(image.width(), image.height());
	
	RGBtoGray(image,grayImage);
	gbImage = grayImage.get_blur_median(3);
	gbImage.blur(3);
	//gbImage.vanvliet(0.5,0);
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
    
	vector<pix> hough_line = findLines(hough_vote, image);
	printOutHough(hough_line,hough_vote);
	printOutLineAndDots(hough_line,image);
	CImg<float> backupImage(filename);
	CImg<float> res = warp(intersec,backupImage);
	res.save("warp.jpg");
	cout<<filename<<endl;
	return 0;
}
