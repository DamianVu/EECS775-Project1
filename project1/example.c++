// example.c++ -- Illustrate use of ImageReader and ImageWriter
// This program will read the image whose file name is in argv[1],
// shift colors R->G, G->B, B->R, and then write the file to argv[2].

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#include "ImageReader.h"
#include "ImageWriter.h"

cryph::Packed3DArray<unsigned char>* shiftColors(ImageReader* ir)
{
	cryph::Packed3DArray<unsigned char>* rgb = ir->getInternalPacked3DArrayImage();
	int xres = ir->getWidth();
	int yres = ir->getHeight();
	for (int r=0 ; r<yres ; r++)
		for (int c=0 ; c<xres ; c++)
		{
			unsigned char R = rgb->getDataElement(r, c, 0);
			unsigned char G = rgb->getDataElement(r, c, 1);
			unsigned char B = rgb->getDataElement(r, c, 2);
			// (This simple demo does not modify the alpha channel, if present.)
			// Circular shift the colors:
			rgb->setDataElement(r, c, 0, G);
			rgb->setDataElement(r, c, 1, B);
			rgb->setDataElement(r, c, 2, R);
		}
	return rgb;
}

const unsigned char* createAndWriteTheOutputFile(ImageReader* ir, ImageWriter* iw)
{
	// This simple demo program assumes the size of input and output images are
	// the same, and it simply modifies the input data by shifting the colors.
	// In general, you will (i) dynamically allocate a Packed3DArray of size
	// (xresOut, yresOut, nChannelsOut) and (ii) fill it from the Packed3DArray
	// retrieved from the ImageReader.
	// Recall you can query the xres, yres, and numChannels of both images, either
	// from the ImageReader/ImageWriter or equivalently from the Packed3DArray
	// returned by the >getInternalPacked3DArrayImage() method of class ImageReader.
	cryph::Packed3DArray<unsigned char>* pa = shiftColors(ir);
	return pa->getData();
}

int main(int argc, char* argv[])
{
	if (argc < 3)
		cerr << "Usage: " << argv[0] << " inputImageFile outputImageFile\n";
	else
	{
		ImageReader* ir = ImageReader::create(argv[1]);
		if (ir == nullptr)
			exit(1);
		int xresIn = ir->getWidth();
		int yresIn = ir->getHeight();
		int nChannelsIn = ir->getNumChannels();

		// Create the output image:
		int xresOut = xresIn; // you will generally want a different resolution
		int yresOut = yresIn; // ditto
		int nChannelsOut = nChannelsIn; // probably WILL be the same
		ImageWriter* iw = ImageWriter::create(argv[2], xresOut, yresOut, nChannelsOut);
		if (iw == nullptr)
			exit(1);

		const unsigned char* fb = createAndWriteTheOutputFile(ir, iw);

		delete ir;
		iw->writeImage(fb);
		iw->closeImageFile();
		delete iw;
	}
	return 0;
}
