del smallpt.exe
g++ -o smallpt.exe smallpt.cpp
del image.ppm
smallpt.exe
ppmgl.exe image.ppm