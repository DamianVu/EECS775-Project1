CPP = g++ -std=c++11 -fPIC

INCLUDES = -I../Packed3DArray

CFLAGS = -O -c $(INCLUDES)

OBJS = ImageReader.o BMPImageReader.o BMPLoader.o JPEGImageReader.o TGAImageReader.o PNGImageReader.o

../lib/libCOGLImageReader.so : libCOGLImageReader.so
	cp libCOGLImageReader.so ../lib/

libCOGLImageReader.so: $(OBJS)
	g++ -shared -o libCOGLImageReader.so $(OBJS) -lpng -ljpeg

ImageReader.o: ImageReader.c++
	$(CPP) $(CFLAGS) ImageReader.c++

BMPImageReader.o: BMPImageReader.c++
	$(CPP) $(CFLAGS) BMPImageReader.c++

BMPLoader.o: BMPLoader.c++
	$(CPP) $(CFLAGS) BMPLoader.c++

JPEGImageReader.o: JPEGImageReader.c++
	$(CPP) $(CFLAGS) JPEGImageReader.c++

PNGImageReader.o: PNGImageReader.c++
	$(CPP) $(CFLAGS) PNGImageReader.c++

TGAImageReader.o: TGAImageReader.c++
	$(CPP) $(CFLAGS) TGAImageReader.c++
