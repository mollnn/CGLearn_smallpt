del bin\smallpt.exe
g++ -o bin\smallpt.exe src\smallpt.cpp
del img\image.ppm
bin\smallpt.exe img\image.ppm
bin\ppmviewer.exe img\image.ppm
