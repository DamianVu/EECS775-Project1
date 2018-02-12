
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
using namespace std;

#include "ImageReader.h"
#include "ImageWriter.h"

float cubicInterpolation(float* p, float x) {
	float p0 = p[3] - p[2] - p[0] + p[1];
	float p1 = p[0] - p[1] - p0;
	float p2 = p[2] - p[0];
	float p3 = p[1];
	float x2 = x * x; // Optimization
	return (p0 * x * x2 + p1 * x2 + p2 * x + p3);
}

// Saturate function since we wish to keep all calculations in float
inline unsigned char saturate( float x ) {
	if (x > 255.0f) {
		return 255;
	} else if (x < 0.0f) {
		return 0;
	} else {
		return (unsigned char)x;
	}
}

unsigned char bicubicInterpolation(float** p, float x, float y) {
	// Since bicubic interpolation is basically cubic interpolation in two dimensions, if we do 4 cubics then one cubic of the 4, we get what we want.
	float temp[4];
	temp[0] = cubicInterpolation(p[0], y);
	temp[1] = cubicInterpolation(p[1], y);
	temp[2] = cubicInterpolation(p[2], y);
	temp[3] = cubicInterpolation(p[3], y);

	return saturate(cubicInterpolation(temp, x));
}

int roundNum(float input) { return (input - floor(input < .5)) ? floor(input) : ceil(input); }

cryph::Packed3DArray<unsigned char>* createNewImage(ImageReader* ir, int method, int newXres, int newYres, int newChannels) {

	cryph::Packed3DArray<unsigned char>* orig = ir->getInternalPacked3DArrayImage();
	cryph::Packed3DArray<unsigned char>* out = new cryph::Packed3DArray<unsigned char>(newYres, newXres, newChannels);

	int oldXres = ir->getWidth();
	int oldYres = ir->getHeight();
	int oldChannels = ir->getNumChannels();

	float xRatio = oldXres / (float)newXres;
	float yRatio = oldYres / (float)newYres;

	if (method == 1) {
		// Nearest
		for (int i = 0; i < newYres; i++)
			for (int j = 0; j < newXres; j++) {
				int oldX = roundNum((float)j * xRatio);
				int oldY = roundNum((float)i * yRatio);

				for (int k = 0; k < newChannels; k++)
					out->setDataElement(i, j, k, orig->getDataElement(oldY, oldX, k));
			}
	} else if (method == 2) {
		// Bilinear
		for (int i = 0; i < newYres; i++)
			for (int j = 0; j < newXres; j++) {
				int xi = floor(j * xRatio);
				int yi = floor(i * yRatio);
				float xf = (j * xRatio) - xi;
				float yf = (i * yRatio) - yi;

				// Flags for if we're at the last pixel. We will essentially duplicate the previous pixel
				bool xDuplicateLastPixelFlag = (xi == oldXres - 1);
				bool yDuplicateLastPixelFlag = (yi == oldYres - 1);

				for (int k = 0; k < newChannels; k++) {
					int xk = (xDuplicateLastPixelFlag) ? xi : xi + 1;
					int yk = (yDuplicateLastPixelFlag) ? yi : yi + 1;
					unsigned char value = ((1.f - yf) * (((1.f-xf)*orig->getDataElement(yi, xi, k)) + (xf * orig->getDataElement(yi, xk, k)))) + (yf * (((1.f-xf)*orig->getDataElement(yk, xi, k)) + (xf*orig->getDataElement(yk, xk, k))));
					out->setDataElement(i, j, k, value);
				}
			}
	} else if (method == 3) {
		// Bicubic
		for (int i = 0; i < newYres; i++)
			for (int j = 0; j < newXres; j++) {
				// Find the location in the old image first.
				float x = j * xRatio;
				float y = i * yRatio;
				int xi = floor(x);
				int yi = floor(y);
				float xf = x - xi;
				float yf = y - yi;

				bool xDuplicateFirstPixelFlag = (xi == 0);
				bool yDuplicateFirstPixelFlag = (yi == 0);
				bool xDuplicateThirdPixelFlag = (xi >= oldXres - 2);
				bool yDuplicateThirdPixelFlag = (yi >= oldYres - 2);
				bool xDuplicateLastPixelFlag = (xi >= oldXres - 3);
				bool yDuplicateLastPixelFlag = (yi >= oldYres - 3);

				int xInitial = (xDuplicateFirstPixelFlag) ? xi : xi - 1;
				int yInitial = (yDuplicateFirstPixelFlag) ? yi : yi - 1;
				int xThird = (xDuplicateThirdPixelFlag) ? xi : xi + 1;
				int yThird = (yDuplicateThirdPixelFlag) ? yi : yi + 1;
				int xFinal = (xDuplicateLastPixelFlag) ? xi : xi + 2;
				int yFinal = (yDuplicateLastPixelFlag) ? yi : yi + 2;

				// Let's get all our p values.
				for (int k = 0; k < newChannels; k++) {
					float** currentPvalues = new float*[4];
					for (int l = 0; l < 4; l++) 
						currentPvalues[l] = new float[4];

					currentPvalues[0][0] = orig->getDataElement(yInitial, xInitial, k);
					currentPvalues[0][1] = orig->getDataElement(yInitial, xi, k);
					currentPvalues[0][2] = orig->getDataElement(yInitial, xThird, k);
					currentPvalues[0][3] = orig->getDataElement(yInitial, xFinal, k);

					currentPvalues[1][0] = orig->getDataElement(yi, xInitial, k);
					currentPvalues[1][1] = orig->getDataElement(yi, xi, k);
					currentPvalues[1][2] = orig->getDataElement(yi, xThird, k);
					currentPvalues[1][3] = orig->getDataElement(yi, xFinal, k);

					currentPvalues[2][0] = orig->getDataElement(yThird, xInitial, k);
					currentPvalues[2][1] = orig->getDataElement(yThird, xi, k);
					currentPvalues[2][2] = orig->getDataElement(yThird, xThird, k);
					currentPvalues[2][3] = orig->getDataElement(yThird, xFinal, k);

					currentPvalues[3][0] = orig->getDataElement(yFinal, xInitial, k);
					currentPvalues[3][1] = orig->getDataElement(yFinal, xi, k);
					currentPvalues[3][2] = orig->getDataElement(yFinal, xThird, k);
					currentPvalues[3][3] = orig->getDataElement(yFinal, xFinal, k);


					out->setDataElement(i, j, k, bicubicInterpolation(currentPvalues, xf, yf));

					for (int l = 0; l < 4; l++) 
						delete[] currentPvalues[l];
					delete[] currentPvalues;
				}
			}
	}
	delete orig;
	return out;
}

const unsigned char* createAndWriteTheOutputFile(ImageReader* ir, int method, int newXres, int newYres, int newChannels) {

	cryph::Packed3DArray<unsigned char>* out = createNewImage(ir, method, newXres, newYres, newChannels);
	return out->getData();

}

int main(int argc, char* argv[]) {
	if (argc < 4 || (strcmp(argv[1],"c") != 0 && strcmp(argv[1],"n") != 0 && strcmp(argv[1],"l") != 0)) {
		cerr << "Usage: " << argv[0] << " interpolationType(n,l,c) inputImageFile outputImageFile\n";
	} else {
		//cout << "Args: " << argv[1] << argv[2] << argv[3];
		int method;
		if (strcmp(argv[1],"n") == 0) method = 1;
		else if (strcmp(argv[1],"l") == 0) method = 2;
		else if (strcmp(argv[1],"c") == 0) method = 3;

		ImageReader* ir = ImageReader::create(argv[2]);
		if (ir == nullptr) exit(1);
		int xResIn = ir->getWidth();
		int yResIn = ir->getHeight();
		int nChannelsIn = ir-> getNumChannels();
		cout << "Dimensions of input image: " << xResIn << " x " << yResIn << endl;
		string inWidth, inHeight;
		cout << "Enter desired resampled width: ";
		cin >> inWidth;
		cout << "Enter desired resampled height: ";
		cin >> inHeight;

		int xResOut = stoi(inWidth);
		int yResOut = stoi(inHeight);
		int nChannelsOut = nChannelsIn;
		ImageWriter* iw = ImageWriter::create(argv[3], xResOut, yResOut, nChannelsOut);
		if (iw == nullptr) exit(1);

		const unsigned char* out = createAndWriteTheOutputFile(ir, method, xResOut, yResOut, nChannelsOut);

		float xRat = xResIn / (float)xResOut;
		float yRat = yResIn / (float)yResOut;
		delete ir;
		iw->writeImage(out);
		iw->closeImageFile();
		delete iw;
	}

	return 0;
}