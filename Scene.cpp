#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;

#define PI 3.14159265
/*
	Transformations, clipping, culling, rasterization are done here.
	You can define helper functions inside Scene class implementation.
*/

vector<vector<Vec3>> verticesOfVertices;
vector<vector<Vec3>> viewportVertices;
vector<vector<Vec3>> CVVVertices;

Vec3 vectorA;
Vec3 vectorB;
bool isClipped;

Vec3 constructVector(Vec3 u)
{
	double temp_min = ABS(u.x);
	Vec3 *v = new Vec3(0, -u.z, u.y, -1);

	if (ABS(u.y) < temp_min)
	{
		temp_min = ABS(u.y);
		v->x = -u.z;
		v->y = 0;
		v->z = u.x;
	}
	else if (ABS(u.z) < temp_min)
	{
		v->x = -u.y;
		v->y = u.x;
		v->z = 0;
	}
	return normalizeVec3(*v);
}

double slopeCalculator(Vec3 a, Vec3 b)
{

	if (a.x != b.x)
		return (a.y - b.y) / (a.x - b.x);
	else
	{
		return __DBL_MAX__;
	}
}

void Scene::draw(int x, int y, Vec3 a, Vec3 b)
{
	double alphaX = (x - a.x) / (b.x - a.x);
	double alphaY = (y - a.y) / (b.y - a.y);

	Color *color_a = colorsOfVertices[a.colorId - 1];
	Color *color_b = colorsOfVertices[b.colorId - 1];
	double cX_r = (1 - alphaX) * (color_a->r) + alphaX * color_b->r;
	double cX_g = (1 - alphaX) * (color_a->g) + alphaX * color_b->g;
	double cX_b = (1 - alphaX) * (color_a->b) + alphaX * color_b->b;

	if (cX_r > 255)
		cX_r = 255;
	else if (cX_r < 0)
		cX_r = 0;
	if (cX_g > 255)
		cX_g = 255;
	else if (cX_g < 0)
		cX_g = 0;
	if (cX_b > 255)
		cX_b = 255;
	else if (cX_b < 0)
		cX_b = 0;

	Color c = Color(cX_r, cX_g, cX_b);

	image[x][y].r = c.r;
	image[x][y].g = c.g;
	image[x][y].b = c.b;
}

double f01(double x, double y, Triangle tri, int i)
{ // i -> modelId

	return (x * (viewportVertices[i][tri.vertexIds[0]].y - viewportVertices[i][tri.vertexIds[1]].y) +
			y * (viewportVertices[i][tri.vertexIds[1]].x - viewportVertices[i][tri.vertexIds[0]].x) +
			viewportVertices[i][tri.vertexIds[0]].x * viewportVertices[i][tri.vertexIds[1]].y -
			viewportVertices[i][tri.vertexIds[0]].y * viewportVertices[i][tri.vertexIds[1]].x);
}

double f12(double x, double y, Triangle tri, int i)
{ // i -> modelId

	return (x * (viewportVertices[i][tri.vertexIds[1]].y - viewportVertices[i][tri.vertexIds[2]].y) +
			y * (viewportVertices[i][tri.vertexIds[2]].x - viewportVertices[i][tri.vertexIds[1]].x) +
			viewportVertices[i][tri.vertexIds[1]].x * viewportVertices[i][tri.vertexIds[2]].y -
			viewportVertices[i][tri.vertexIds[1]].y * viewportVertices[i][tri.vertexIds[2]].x);
}

double f20(double x, double y, Triangle tri, int i)
{ // i -> modelId

	return (x * (viewportVertices[i][tri.vertexIds[2]].y - viewportVertices[i][tri.vertexIds[0]].y) +
			y * (viewportVertices[i][tri.vertexIds[0]].x - viewportVertices[i][tri.vertexIds[2]].x) +
			viewportVertices[i][tri.vertexIds[2]].x * viewportVertices[i][tri.vertexIds[0]].y -
			viewportVertices[i][tri.vertexIds[2]].y * viewportVertices[i][tri.vertexIds[0]].x);
}

bool Scene::clipping(Vec3 a, Vec3 b, Camera c)
{
	isClipped = false;

	int abit = 0;
	int bbit = 0;
	if (a.x < 0)
		abit += 1; //else abit[3] = 0;
	if (a.x > c.horRes)
		abit += 2; // else abit[2] = 0;
	if (a.y > c.verRes)
		abit += 8; //else abit[0] = 0;
	if (a.y < 0)
		abit += 4; // else abit[1] = 0;

	if (b.x < 0)
		bbit += 1; // else bbit[3] = 0;
	if (b.x > c.horRes)
		bbit += 2; //else bbit[2] = 0;
	if (b.y > c.verRes)
		bbit += 8; //else bbit[0] = 0;
	if (b.y < 0)
		bbit += 4; // else bbit[1] = 0;

	if ((abit | bbit) == 0)
		return false;
	else if (abit & bbit)
	{ // trivial reject
		vectorA = Vec3(0, 0, 0, -1);
		vectorB = Vec3(0, 0, 0, -1);
		isClipped = true;
		return true;
	}
	else
	{
		int xxx;
		double x, y;

		if (abit != 0)
			xxx = abit;
		else
			xxx = bbit;

		//intersection

		if (xxx & 8)
		{
			x = a.x + (b.x - a.x) * (c.verRes - 1 - a.y) / (b.y - a.y);
			y = c.verRes - 1;
		}
		else if (xxx & 4)
		{
			x = a.x + (b.x - a.x) * (-a.y) / (b.y - a.y);
			y = 0;
		}
		else if (xxx & 2)
		{
			y = a.y + (b.y - a.y) * (c.horRes - 1 - a.x) / (b.x - a.x);
			x = c.horRes - 1;
		}
		else if (xxx & 1)
		{
			y = a.y + (b.y - a.y) * (-a.x) / (b.x - a.x);
			x = 0;
		}
		if (xxx == abit)
		{
			vectorA = Vec3(x, y, 1, a.colorId);
			vectorB = Vec3(b.x, b.y, b.z, b.colorId);
			return true;
		}

		else
		{
			vectorB = Vec3(x, y, 1, b.colorId);
			vectorA = Vec3(a.x, a.y, a.z, a.colorId);
			return true;
		}
	}
}

void Scene::midPointF(int i, int j, int modelId, Camera c)
{
	Triangle tri = models[i]->triangles[j];

	double m;
	for (int i = 0; i < 3; i++)
	{

		Vec3 a = viewportVertices[modelId - 1][tri.vertexIds[i]];
		Vec3 b = viewportVertices[modelId - 1][tri.vertexIds[(i + 1) % 3]];

		//clipping varsa
		if (clipping(a, b, c) && (isClipped == false))
		{
			a = vectorA;
			b = vectorB;
		}
		if (isClipped)
			continue;
		m = slopeCalculator(a, b);
		if (m == 0)
			continue;
		if ((b.y - a.y) >= 0 && (b.x - a.x) >= 0)
		{ // 1.çeyrek
			if (m > 0 && m <= 1)
			{ //1.bölge
				int x = a.x;
				int y = a.y;
				// y=a.y;
				double M = 1.0 * (a.y - b.y) + 0.5 * (b.x - a.x);
				for (x = a.x; x < int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x + 1)
				{
					draw(x, y, a, b);
					if (M < 0)
					{ // NE
						y += 1;
						M += ((a.y - b.y) + (b.x - a.x));
					}
					else
					{ // E
						M += (a.y - b.y);
					}
				}
			}
			else if (m > 1)
			{ //2.bölge
				int x = a.x;
				int y = a.y;
				double M = 0.5 * (a.y - b.y) + 1.0 * (b.x - a.x);
				for (y = a.y; y < int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y + 1)
				{
					draw(x, y, a, b);
					if (M <= 0)
					{ // N
						M += (b.x - a.x);
					}
					else
					{ // NE
						x += 1;
						M += ((a.y - b.y) + (b.x - a.x));
					}
				}
			}
		}
		else if ((b.y - a.y) > 0 && (b.x - a.x) < 0) // 2.çeyrek
		{
			if (m < -1)
			{ //2.bölge
				int x = a.x;
				int y = a.y;
				double M = -0.5 * (a.y - b.y) + 1.0 * (b.x - a.x);
				for (y = a.y; y < int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y + 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // N
						M += ((b.x - a.x));
					}
					else
					{ // NW
						x -= 1;
						M += (-(a.y - b.y) + (b.x - a.x));
					}
				}
			}
			else if (m < 0 && m >= -1)
			{ // 1.bölge
				int x = a.x;
				int y = a.y;
				double M = -1.0 * (a.y - b.y) + 0.5 * (b.x - a.x);
				for (x = a.x; x > int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x - 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // NW
						y += 1;
						M += (-(a.y - b.y) + (b.x - a.x));
					}
					else
					{ // W
						M += (-(a.y - b.y));
					}
				}
			}
		}
		else if ((b.y - a.y) <= 0 && (b.x - a.x) <= 0) // 3.çeyrek
		{
			if (m > 1)
			{ // 2.bölge
				int x = a.x;
				int y = a.y;
				double M = -0.5 * (a.y - b.y) - 1.0 * (b.x - a.x);
				for (y = a.y; y > int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y - 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // SW
						x -= 1;
						M += (-(a.y - b.y) - (b.x - a.x));
					}
					else
					{ // S
						M += (-(b.x - a.x));
					}
				}
			}
			else if (m > 0 && m <= 1)
			{ //1.bölge
				int x = a.x;
				int y = a.y;
				double M = -1.0 * (a.y - b.y) - 0.5 * (b.x - a.x);
				for (x = a.x; x > int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x - 1)
				{
					draw(x, y, a, b);
					if (M > 0)
					{ // W
						M += (-(a.y - b.y));
					}
					else
					{ // SW
						y -= 1;
						M += (-(a.y - b.y) - (b.x - a.x));
					}
				}
			}
		}
		else if ((b.y - a.y) < 0 && (b.x - a.x) > 0) //4.çeyrek
		{
			if (m < 0 && m >= -1)
			{ //1.bölge
				int x = a.x;
				int y = a.y;
				double M = 1.0 * (a.y - b.y) - 0.5 * (b.x - a.x);
				for (x = a.x; x < int(b.x) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; x = x + 1)
				{
					draw(x, y, a, b);
					if (M <= 0)
					{ // E
						M += ((a.y - b.y));
					}
					else
					{ // SE
						y -= 1;
						M += ((a.y - b.y) - (b.x - a.x));
					}
				}
			}
			else if (m < -1)
			{ //2.bölge
				int x = a.x;
				int y = a.y;
				double M = 0.5 * (a.y - b.y) - 1.0 * (b.x - a.x);
				for (y = a.y; y > int(b.y) && x < c.horRes && x >= 0 && y < c.verRes && y >= 0; y = y - 1)
				{
					draw(x, y, a, b);
					if (M < 0)
					{ // SE
						x += 1;
						M += ((a.y - b.y) - (b.x - a.x));
					}
					else
					{ // S
						M += (-(b.x - a.x));
					}
				}
			}
		}
	}
}

void Scene::rasterization(Triangle tri, int modelId, Camera c)
{
	double xmin = std::min(viewportVertices[modelId - 1][tri.vertexIds[0]].x,
						   std::min(viewportVertices[modelId - 1][tri.vertexIds[1]].x, viewportVertices[modelId - 1][tri.vertexIds[2]].x));
	double ymin = std::min(viewportVertices[modelId - 1][tri.vertexIds[0]].y,
						   std::min(viewportVertices[modelId - 1][tri.vertexIds[1]].y, viewportVertices[modelId - 1][tri.vertexIds[2]].y));
	double xmax = std::max(viewportVertices[modelId - 1][tri.vertexIds[0]].x,
						   std::max(viewportVertices[modelId - 1][tri.vertexIds[1]].x, viewportVertices[modelId - 1][tri.vertexIds[2]].x));
	double ymax = std::max(viewportVertices[modelId - 1][tri.vertexIds[0]].y,
						   std::max(viewportVertices[modelId - 1][tri.vertexIds[1]].y, viewportVertices[modelId - 1][tri.vertexIds[2]].y));

	ymin = std::max((int)ymin, 0);
	xmin = std::max((int)xmin, 0);
	ymax = std::min((int)ymax, c.verRes - 1);
	xmax = std::min((int)xmax, c.horRes - 1);
	if (ymax < 0)
		ymax = 0;
	if (xmax < 0)
		xmax = 0;
	if (ymin < 0)
		ymin = 0;
	if (xmin < 0)
		xmin = 0;

	for (int j = ymin; j <= ymax; j = j + 1)
	{
		for (int i = xmin; i <= xmax; i = i + 1)
		{

			double alpha = f12(i, j, tri, modelId - 1) / f12(viewportVertices[modelId - 1][tri.vertexIds[0]].x,
															 viewportVertices[modelId - 1][tri.vertexIds[0]].y, tri, modelId - 1);
			double beta = f20(i, j, tri, modelId - 1) / f20(viewportVertices[modelId - 1][tri.vertexIds[1]].x,
															viewportVertices[modelId - 1][tri.vertexIds[1]].y, tri, modelId - 1);
			double gama = f01(i, j, tri, modelId - 1) / f01(viewportVertices[modelId - 1][tri.vertexIds[2]].x,
															viewportVertices[modelId - 1][tri.vertexIds[2]].y, tri, modelId - 1);

			if (alpha >= 0 && beta >= 0 && gama >= 0)
			{
				Color c0 = *colorsOfVertices[tri.vertexIds[0] - 1];
				c0.r *= alpha;
				c0.g *= alpha;
				c0.b *= alpha;
				Color c1 = *colorsOfVertices[tri.vertexIds[1] - 1];
				c1.r *= beta;
				c1.g *= beta;
				c1.b *= beta;
				Color c2 = *colorsOfVertices[tri.vertexIds[2] - 1];
				c2.r *= gama;
				c2.g *= gama;
				c2.b *= gama;
				Color c;
				c.r = c0.r + c1.r + c2.r;
				c.g = c0.g + c1.g + c2.g;
				c.b = c0.b + c1.b + c2.b;

				if (c.r > 255)
					c.r = 255;
				if (c.g > 255)
					c.g = 255;
				if (c.b > 255)
					c.b = 255;

				image[i][j].r = (c.r);
				image[i][j].g = (c.g);
				image[i][j].b = (c.b);
			}
		}
	}
}

Matrix4 Scene::rotation(Matrix4 res, Rotation rotation)
{
	Vec3 u = Vec3(rotation.ux, rotation.uy, rotation.uz, -1);
	Vec3 v = constructVector(u);
	Vec3 w = crossProductVec3(u, v);
	w.colorId = -1;

	double M[4][4] = {{u.x, u.y, u.z, 0},
					  {v.x, v.y, v.z, 0},
					  {w.x, w.y, w.z, 0},
					  {0, 0, 0, 1}};

	double inverse_M[4][4] = {{u.x, v.x, w.x, 0},
							  {u.y, v.y, w.y, 0},
							  {u.z, v.z, w.z, 0},
							  {0, 0, 0, 1}};

	double theta = rotation.angle * (PI / 180.0);
	double R_x[4][4] = {{1, 0, 0, 0},
						{0, cos(theta), -sin(theta), 0},
						{0, sin(theta), cos(theta), 0},
						{0, 0, 0, 1}};

	Matrix4 r_x_m;
	r_x_m = multiplyMatrixWithMatrix(Matrix4(R_x), Matrix4(M));
	res = multiplyMatrixWithMatrix(Matrix4(inverse_M), r_x_m);
	return res;
}

Matrix4 Scene::scaling(Matrix4 res, Scaling scaling)
{
	double temp[4][4] = {{scaling.sx, 0, 0, 0},
						 {0, scaling.sy, 0, 0},
						 {0, 0, scaling.sz, 0},
						 {0, 0, 0, 1}};
	res = Matrix4(temp);
	return res;
}

Matrix4 Scene::translation(Matrix4 res, Translation translation)
{
	double temp[4][4] = {{1, 0, 0, translation.tx},
						 {0, 1, 0, translation.ty},
						 {0, 0, 1, translation.tz},
						 {0, 0, 0, 1}};
	res = Matrix4(temp);
	return res;
}

Matrix4 Scene::viewingTransform(Camera c)
{
	double cam[4][4] = {{c.u.x, c.u.y, c.u.z, -(c.u.x * c.pos.x + c.u.y * c.pos.y + c.u.z * c.pos.z)},
						{c.v.x, c.v.y, c.v.z, -(c.v.x * c.pos.x + c.v.y * c.pos.y + c.v.z * c.pos.z)},
						{c.w.x, c.w.y, c.w.z, -(c.w.x * c.pos.x + c.w.y * c.pos.y + c.w.z * c.pos.z)},
						{0, 0, 0, 1}};

	double otrh[4][4] = {{(2) / (c.right - c.left), 0, 0, -(c.top + c.bottom) / (c.top - c.bottom)},
						 {0, 2 / (c.top - c.bottom), 0, -(c.top + c.bottom) / (c.top - c.bottom)},
						 {0, 0, -2 / (c.far - c.near), -(c.far + c.near) / (c.far - c.near)},
						 {0, 0, 0, 1}};

	double per[4][4] = {{(2 * c.near) / (c.right - c.left), 0, (c.right + c.left) / (c.right - c.left), 0},
						{0, 2 * c.near / (c.top - c.bottom), (c.top + c.bottom) / (c.top - c.bottom), 0},
						{0, 0, -(c.far + c.near) / (c.far - c.near), -(2 * c.far * c.near) / (c.far - c.near)},
						{0, 0, -1, 0}};

	Matrix4 perMtx = Matrix4(per);
	Matrix4 camMtx = Matrix4(cam);
	Matrix4 othMtx = Matrix4(otrh);
	Matrix4 ret;
	if (projectionType == 1) //perspective
		ret = multiplyMatrixWithMatrix(perMtx, camMtx);
	else
		ret = multiplyMatrixWithMatrix(othMtx, camMtx);

	return ret;
}

//0 not culled 1 culled
int Scene::culling(int modelId, Camera cam, Triangle tri)
{
	Vec3 a = CVVVertices[modelId - 1][tri.vertexIds[0]];
	Vec3 b = CVVVertices[modelId - 1][tri.vertexIds[1]];
	Vec3 c = CVVVertices[modelId - 1][tri.vertexIds[2]];
	Vec3 center = multiplyVec3WithScalar(addVec3(addVec3(a, b), c), 1.0 / 3.0);

	Vec3 v = subtractVec3(center, cam.pos);
	Vec3 n = crossProductVec3(subtractVec3(b, a), subtractVec3(c, a));

	if (dotProductVec3(n, center) > 0) // backface ise
		return 1;
	else // front ise
		return 0;
}

void Scene::forwardRenderingPipeline(Camera *cam)
{
	vertex_processing();
	for (int x = 0; x < verticesOfVertices.size(); x++)
	{
		vector<Vec3> vect(verticesOfVertices[x]);
		viewportVertices.push_back(vect);
		CVVVertices.push_back(vect);
	}

	Matrix4 viewingMatrix;
	viewingMatrix = viewingTransform(*cam);
	double viewportMatrix[4][4] = {{(cam->horRes) / 2.0, 0, 0, (cam->horRes - 1) / 2.0},
								   {0, (cam->verRes) / 2.0, 0, (cam->verRes - 1) / 2.0},
								   {0, 0, 0.5, 0.5},
								   {0, 0, 0, 0}};
	Vec4 temp, temp2;

	for (int i = 0; i < verticesOfVertices.size(); i++)
	{
		for (int j = 0; j < verticesOfVertices[i].size(); j++)
		{
			temp = Vec4(verticesOfVertices[i][j].x, verticesOfVertices[i][j].y, verticesOfVertices[i][j].z, 1, verticesOfVertices[i][j].colorId);
			temp2 = multiplyMatrixWithVec4(viewingMatrix, temp);
			temp2.x /= temp2.t;
			temp2.y /= temp2.t;
			temp2.z /= temp2.t;
			temp2.t = 1;

			Vec3 v2 = Vec3(temp2.x, temp2.y, temp2.z, verticesOfVertices[i][j].colorId);
			CVVVertices[i][j] = v2;

			temp = multiplyMatrixWithVec4(viewportMatrix, temp2);
			Vec3 v3 = Vec3(temp.x, temp.y, temp.z, verticesOfVertices[i][j].colorId);
			viewportVertices[i][j] = v3;
		}
	}
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			if (cullingEnabled && culling(i + 1, *cam, models[i]->triangles[j]) == 0)
				continue;
			if (models[i]->type == 0) // wireframe
			{
				midPointF(i, j, models[i]->modelId, *cam);
			}
			else
			{
				rasterization(models[i]->triangles[j], models[i]->modelId, *cam);
			}
		}
	}
}

void Scene::vertex_processing()
{
	for (int i = 0; i < models.size(); i++)
	{
		Matrix4 modelingMatrix;
		Vec4 hom_point;
		Vec4 result_point;
		modelingMatrix = modelingTransform(*models.at(i));
		vector<Vec3> newVertices;
		Vec3 dummy_vector = Vec3(0, 0, 0, -1);
		newVertices.push_back(dummy_vector);
		for (int j = 0; j < vertices.size(); j++)
		{
			hom_point = Vec4(vertices[j]->x, vertices[j]->y, vertices[j]->z, 1, -1);
			result_point = multiplyMatrixWithVec4(modelingMatrix, hom_point);
			Vec3 v = Vec3(result_point.x, result_point.y, result_point.z, vertices[j]->colorId);
			newVertices.push_back(v);
		}
		verticesOfVertices.push_back(newVertices);
	}
}

Matrix4 Scene::modelingTransform(Model model)
{
	Matrix4 temp, *temp2, temp3;
	double temp_matrix[4][4] = {{1, 0, 0, 0},
								{0, 1, 0, 0},
								{0, 0, 1, 0},
								{0, 0, 0, 1}};
	temp2 = new Matrix4(temp_matrix);

	for (int i = 0; i < model.numberOfTransformations; i++)
	{
		char type = model.transformationTypes[i];
		if (type == 'r')
		{
			temp = rotation(temp, *rotations[model.transformationIds[i] - 1]);
			temp3 = multiplyMatrixWithMatrix(temp, *temp2);
			temp2 = new Matrix4(temp3);
		}
		else if (type == 's')
		{
			temp = scaling(temp, *scalings[model.transformationIds[i] - 1]);
			temp3 = multiplyMatrixWithMatrix(temp, *temp2);
			temp2 = new Matrix4(temp3);
		}
		else if (type == 't')
		{
			temp = translation(temp, *translations[model.transformationIds[i] - 1]);
			temp3 = multiplyMatrixWithMatrix(temp, *temp2);
			temp2 = new Matrix4(temp3);
		}
		else
		{
			printf("Error on modeling transform \n");
		}
	}
	return *temp2;
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL)
		pElement->QueryBoolText(&cullingEnabled);

	// read projection type
	pElement = pRoot->FirstChildElement("ProjectionType");
	if (pElement != NULL)
		pElement->QueryIntText(&projectionType);

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read models
	pElement = pRoot->FirstChildElement("Models");

	XMLElement *pModel = pElement->FirstChildElement("Model");
	XMLElement *modelElement;
	while (pModel != NULL)
	{
		Model *model = new Model();

		pModel->QueryIntAttribute("id", &model->modelId);
		pModel->QueryIntAttribute("type", &model->type);

		// read model transformations
		XMLElement *pTransformations = pModel->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		pTransformations->QueryIntAttribute("count", &model->numberOfTransformations);

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			model->transformationTypes.push_back(transformationType);
			model->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		// read model triangles
		XMLElement *pTriangles = pModel->FirstChildElement("Triangles");
		XMLElement *pTriangle = pTriangles->FirstChildElement("Triangle");

		pTriangles->QueryIntAttribute("count", &model->numberOfTriangles);

		while (pTriangle != NULL)
		{
			int v1, v2, v3;

			str = pTriangle->GetText();
			sscanf(str, "%d %d %d", &v1, &v2, &v3);

			model->triangles.push_back(Triangle(v1, v2, v3));

			pTriangle = pTriangle->NextSiblingElement("Triangle");
		}

		models.push_back(model);

		pModel = pModel->NextSiblingElement("Model");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	// if image is filled before, just change color rgb values with the background color
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}