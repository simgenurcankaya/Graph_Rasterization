#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Matrix4.h"

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;
	int projectionType;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Model* > models;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

	void vertex_processing();
	Matrix4 modelingTransform(Model model);
	Matrix4 rotation(Matrix4, Rotation rot);
	Matrix4 translation(Matrix4, Translation tra);
	Matrix4 scaling(Matrix4, Scaling sca);
	Matrix4 viewingTransform(Camera c);
	int culling(int modelID, Camera cam, Triangle tri);
	void midPointF(int , int, int modelId,Camera);
	void draw(int x, int y, Vec3 a, Vec3 b);
	void rasterization(Triangle tri, int modelId,Camera);
	bool clipping(Vec3 a, Vec3 b,Camera);
	Color colorFinder(int x,int y , Vec3 a, Vec3 b);


};

#endif
