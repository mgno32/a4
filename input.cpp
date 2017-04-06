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
const int LINENUM = 4;
const float max_grad = 400;
const double PI = 3.1415926;
const double DIFF = 200;
const int Threshold = 300;
int main(int agrc,char **argv) {
	char *filename = argv[1];
	//cout<<filename<<" "<<&filename<<endl;
	CImg<uchar> image(filename);
	image.display();
	return 0;
}
