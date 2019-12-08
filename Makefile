rasterizer_cpp:
	g++ -o rasterizer *.cpp
clean:
	rm -rf *.ppm *.png ./rasterizer
