#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

struct Intersection {
	vec4 position;
	float distance;
	int triangleIndex;
};

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

//Variables to store the camera parameters
float focalLength = SCREEN_HEIGHT;
vec4 cameraPos(0, 0, -3, 1);
mat4 cameraRotMatrix;

//Define light by its position and colour
vec4 lightPos(0, -0.5, -0.7, 1.0);
vec3 lightColor = 14.f * vec3(1, 1, 1);

vector<Triangle> triangles;

SDL_Event event;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);

bool ClosestIntersection(vec4 start,
	vec4 dir,
	const vector<Triangle>& triangles,
	Intersection& closestIntersection);

void Rotate(mat3 rotation);
mat3 RotMatrixX(float angle);
mat3 RotMatrixY(float angle);

vec3 DirectLight(const Intersection& i);

double Find3Distance(vec3 vectorOne, vec3 vectorTwo);

void NormaliseVec(vec3 vec);

int main(int argc, char* argv[])
{

	screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
	LoadTestModel(triangles);

	//Fill the cameraRotMatrix with the identity matrix
	vec4 x0(1,0,0,0);
	vec4 x1(0,1,0,0);
	vec4 x2(0,0,1,0);
	vec4 x3(0,0,0,1);
	cameraRotMatrix = mat4(x0,x1,x2,x3);

	while (Update())
	{
		Draw(screen);
		SDL_Renderframe(screen);
	}

	SDL_SaveImage(screen, "screenshot.bmp");

	KillSDL(screen);
	return 0;
}

//Function that takes a ray and finds the closest intersecting geometry
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection) {

	bool closestIntersectionFound = false;

	//Max size of float
	closestIntersection.distance = std::numeric_limits<float>::max();

	//For each triangle find intersection point x
	for (int i = 0; i < triangles.size(); i++) {
		vec4 v0 = triangles[i].v0;
		vec4 v1 = triangles[i].v1;
		vec4 v2 = triangles[i].v2;

		vec3 e1 = vec3(v1 - v0);
		vec3 e2 = vec3(v2 - v0);
		vec3 b = vec3(start - v0);

		mat3 A(-vec3(dir), e1, e2);
		if (glm::determinant(A) != 0) {

			vec3 x = glm::inverse(A) * b;

			float t = x.x;
			float u = x.y;
			float v = x.z;

			if (t >= 0) {

				if ((0 <= u) && (0 <= v) && (u + v) <= 1) {

					if (t < closestIntersection.distance) {
						closestIntersectionFound = true;
						//Set values for the intersection
						closestIntersection.position = start + dir * t;
						closestIntersection.distance = t;
						closestIntersection.triangleIndex = i;
					}
				}
			}

		}
	}
	return closestIntersectionFound;
}

/*Place your drawing here*/
void Draw(screen* screen) {
	//Clear buffer
	memset(screen->buffer, 0, screen->height*screen->width * sizeof(uint32_t));

	//Iterate through every pixel in the image
	for (int y = 0; y < SCREEN_HEIGHT; y++) {
		for (int x = 0; x < SCREEN_WIDTH; x++) {

			//Compute the corresponding ray direction
			vec4 rayDirection(x - SCREEN_WIDTH * 0.5, y - SCREEN_HEIGHT * 0.5, focalLength, 1);

			//
			rayDirection = rayDirection*cameraRotMatrix;

			//Call the function ClosestIntersection to get the closest intersection in that direction
			Intersection closestIntersectionItem;
			if (ClosestIntersection(cameraPos, rayDirection, triangles, closestIntersectionItem)) {
				//The color of the pixel should be set to the color of that triangle
				//PutPixelSDL(screen, x, y, triangles[closestIntersectionItem.triangleIndex].color);
				//The color of the pixel is set to the percentage of light that hits it
				PutPixelSDL(screen, x, y, DirectLight(closestIntersectionItem));

			}
			else {
				//It should be black
				vec3 colour(0.0, 0.0, 0.0);
				PutPixelSDL(screen, x, y, colour);
			}
		}
	}
}

mat3 RotMatrixX(float angle){
	vec3 x0(1, 0, 0);
	vec3 x1(0, cosf(angle), sinf(angle));
	vec3 x2(0, -sinf(angle), cosf(angle));

	return mat3(x0,x1,x2);
}

mat3 RotMatrixY(float angle){
	vec3 y1(cosf(angle), 0, -sinf(angle));
	vec3 y2(0, 1, 0);
	vec3 y3(sinf(angle), 0, cosf(angle));

	return mat3(y1,y2,y3);
}

void Rotate(mat3 rotation){
	mat3 cameraRotExtract(cameraRotMatrix[0][0], cameraRotMatrix[0][1], cameraRotMatrix[0][2],
													 cameraRotMatrix[1][0], cameraRotMatrix[1][1], cameraRotMatrix[1][2],
												 	 cameraRotMatrix[2][0], cameraRotMatrix[2][1], cameraRotMatrix[2][2]);

	cameraRotExtract = rotation * cameraRotExtract;

	cameraRotMatrix[0][0] = cameraRotExtract[0][0];
	cameraRotMatrix[0][1] = cameraRotExtract[0][1];
	cameraRotMatrix[0][2] = cameraRotExtract[0][2];
	cameraRotMatrix[1][0] = cameraRotExtract[1][0];
	cameraRotMatrix[1][1] = cameraRotExtract[1][1];
	cameraRotMatrix[1][2] = cameraRotExtract[1][2];
	cameraRotMatrix[2][0] = cameraRotExtract[2][0];
	cameraRotMatrix[2][1] = cameraRotExtract[2][1];
	cameraRotMatrix[2][2] = cameraRotExtract[2][2];
}

/*Place updates of parameters here*/
bool Update()
{
	static int t = SDL_GetTicks();
	/* Compute frame time */
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;

	float angle = 0.01;
	int step = 1;

	SDL_Event e;
	while (SDL_PollEvent(&e))
	{
		if (e.type == SDL_QUIT)
		{
			return false;
		}
		else
			if (e.type == SDL_KEYDOWN)
			{
				int key_code = e.key.keysym.sym;
				switch (key_code)
				{
				case SDLK_UP:
					/* Move camera forward */
					//cameraPos.z += step;
					/*Rotate camera upwards on X axis */
					Rotate(RotMatrixX(-angle));
					break;
				case SDLK_DOWN:
					/* Move camera backwards */
					//cameraPos.z -= step;
					/*Rotate camera downwards on X axis */
					Rotate(RotMatrixX(angle));
					break;
				case SDLK_LEFT:
					/* Move camera left */
					//cameraPos.x += step;
					/*Rotate camera left on Y axis */
					Rotate(RotMatrixY(angle));
					break;
				case SDLK_RIGHT:
					/* Move camera right */
					//cameraPos.x -= step;
					/*Rotate camera right on Y axis */
					Rotate(RotMatrixY(-angle));
					break;

				//Camera Movement
				case SDLK_w:
					/*Move light forward*/
					lightPos.z += step;
					break;
				case SDLK_s:
					/*Move light backwards*/
					lightPos.z -= step;
					break;
				case SDLK_q:
					/*Move light up*/
					lightPos.y += step;
					break;
				case SDLK_e:
					/*Move light down*/
					lightPos.y -= step;
					break;
				case SDLK_a:
					/*Move light left*/
					lightPos.x -= step;
					break;
				case SDLK_d:
					/*Move light right*/
					lightPos.x += step;
					break;

				case SDLK_ESCAPE:
					/* Move camera quit */
					return false;
				}
			}
	}
	return true;
}

//Function that takes an intersection and returns the color of the triangle
vec3 DirectLight(const Intersection& i) {

	vec3 updatedColor;

	vec3 n = triangles[i.triangleIndex].normal;
	vec3 r;
	r.x = lightPos.x - i.position.x;
	r.y = lightPos.y - i.position.y;
	r.z = lightPos.z - i.position.z;

	//Normalise vectors n and r
	NormaliseVec(n);
	NormaliseVec(r);

	double dotProduct = (n.x * r.x) + (n.y * r.y) + (n.z * r.z);

	if (dotProduct > 0) {

		//Find the distance between the lightPos and the intersection
		double d = Find3Distance(lightPos, vec3(i.position));

		//Calculate the number (0 < n < 1) with which the power P of the light (lightColor) will be multiplied with
		double percentage = dotProduct / (4 * M_PI * d);
		
		//Calculate the new value of light for the intersection
		updatedColor.x = lightColor.x * percentage;
		updatedColor.y = lightColor.y * percentage;
		updatedColor.z = lightColor.z * percentage;
	}
	else {
		updatedColor.x = 0;
		updatedColor.y = 0;
		updatedColor.z = 0;
	}

	return updatedColor;
}

double Find3Distance(vec3 vectorOne, vec3 vectorTwo) {

	double sqrDifferenceX = (vectorOne.x - vectorTwo.x)*(vectorOne.x - vectorTwo.x);
	double sqrDifferenceY = (vectorOne.y - vectorTwo.y)*(vectorOne.y - vectorTwo.y);
	double sqrDifferenceZ = (vectorOne.z - vectorTwo.z)*(vectorOne.z - vectorTwo.z);
	double root = sqrt(sqrDifferenceX + sqrDifferenceY + sqrDifferenceZ);

	return root;
}

void NormaliseVec(vec3 vec) {
	double length = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);

	vec.x = vec.x / length;
	vec.y = vec.y / length;
	vec.z = vec.z / length;
}