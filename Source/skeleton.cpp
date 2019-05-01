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

//#define SCREEN_WIDTH 320
//#define SCREEN_HEIGHT 256
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 500
#define FULLSCREEN_MODE false

//Variables to store the camera parameters
float focalLength = SCREEN_HEIGHT;
vec4 cameraPos(0, 0, -3, 1);
mat4 cameraRotMatrix;

//Define light by its position and colour
vec4 lightPos(0, -0.5, -0.7, 1.0);
vec3 lightColor = 14.f * vec3(1, 1, 1);

vec3 indirectLight = 0.5f*vec3(1, 1, 1);

//Sphere parameters
vec4 sphereCenter(0, 0, 0.02, 1);
const float sphereRadius = 0.2;

//Indeces of refraction
float air = 1;
float glass = 1.5;
float etai = air; 
float etat = glass;
float eta = etai / etat;

vector<Triangle> triangles;

//Count of recursion
int recursionCount = 0;

SDL_Event event;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);

bool ClosestIntersection(vec4 start,
	vec4 dir,
	const vector<Triangle>& triangles,
	Intersection& closestIntersection);

bool sphereClosestIntersection(vec4 origin, vec4 dir, Intersection& closestSphereIntersection);

void Rotate(mat3 rotation);
mat3 RotMatrixX(float angle);
mat3 RotMatrixY(float angle);
vec3 DirectLight(const Intersection& i);
vec3 CombineReflectionRefraction(const Intersection& intersection);
vec3 ReflectionColor(const Intersection& intersection, vec4 i, vec4 n, double dotProduct);
vec3 RefractionColor(const Intersection& intersection, vec4 i, vec4 n, double dotProduct);
vec4 RefractionDirection(const vec4 i, const vec4 n, double dotProduct);
void FresnelR(double costheta1, double costheta2, double fresnelOrthogonal, double fresnelParallel, double fresnel);
double Find3Distance(vec3 vectorOne, vec3 vectorTwo);
vec4 NormaliseVec(vec4 vec);
double angle(vec4 u, vec4 v);
vec3 RayTracerColor(float x, float y);

int main(int argc, char* argv[])
{

	screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
	LoadTestModel(triangles);

	//Fill the cameraRotMatrix with the identity matrix
	vec4 x0(1, 0, 0, 0);
	vec4 x1(0, 1, 0, 0);
	vec4 x2(0, 0, 1, 0);
	vec4 x3(0, 0, 0, 1);
	cameraRotMatrix = mat4(x0, x1, x2, x3);

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
		float aDet = glm::determinant(A);

		if (aDet != 0) {

			//t found using Cramers rule
			mat3 A1(b, e1, e2);
			float x1Det = glm::determinant(A1);
			float t = x1Det / aDet;

			if (t >= 0) {

				mat3 A2(-vec3(dir), b, e2);
				mat3 A3(-vec3(dir), e1, b);

				float x2Det = glm::determinant(A2);
				float x3Det = glm::determinant(A3);

				float u = x2Det / aDet;
				float v = x3Det / aDet;

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
	Intersection intersectionSphere;
	//Check if there is a sphere in front of it
	if (sphereClosestIntersection(start, dir, intersectionSphere)) {
		//Set values for the intersection
		if (intersectionSphere.distance < closestIntersection.distance) {
			closestIntersectionFound = true;
			closestIntersection = intersectionSphere;
			closestIntersection.triangleIndex = -1;
		}
	}
	return closestIntersectionFound;
}

//Function that takes a ray and finds if there is an intersection with the sphere
bool sphereClosestIntersection(vec4 origin, vec4 dir, Intersection& closestSphereIntersection) {

	bool closestSphereIntersectionFound = false;

	//Max size of float
	closestSphereIntersection.distance = std::numeric_limits<float>::max();

	//Find intersection points with sphere if any

	//First calculate a, b, c for f(t) = at^2 + bt + c
	vec3 dirVec3 = dir;
	vec3 originVec3 = origin;
	vec3 center = sphereCenter;

	//float a = glm::dot(dirVec3, dirVec3);
	//float b = 2 * glm::dot(originVec3, dirVec3);
	//float c = glm::dot(originVec3, originVec3) - (sphereRadius * sphereRadius);
	float a = 1.f;
	float b = 2 * glm::dot(originVec3 - center, dirVec3);
	float dist = (originVec3 - center).length();
	float c = dist*dist - sphereRadius * sphereRadius;

	float t;

	//Calculate Delta to see how many intersection points we have

	float Delta = (b*b) - (4 * a * c);

	if (Delta < 0) {
		return closestSphereIntersectionFound;
	}
	else if (Delta == 0) {
		closestSphereIntersectionFound = true;

		t = (-b) / (2 * a);
		//Set values for the intersection
		closestSphereIntersection.position = origin + dir * t;
		closestSphereIntersection.distance = t;
	}
	else { //Delta > 0 and we have 2 values for t. We are going to keep the smallest one because that will be the closest intersection

		float tOne = -(b + sqrt(Delta)) / (2 * a);

		float tTwo = -(b - sqrt(Delta)) / (2 * a);

		if (tOne < tTwo && tOne > 0) {
			closestSphereIntersectionFound = true;
			closestSphereIntersection.position = origin + dir * tOne;
			closestSphereIntersection.distance = tOne;
		}
		else if (tTwo > 0) {
			closestSphereIntersectionFound = true;
			closestSphereIntersection.position = origin + dir * tTwo;
			closestSphereIntersection.distance = tTwo;
		}
	}
	return closestSphereIntersectionFound;
}

void Draw(screen* screen) {
	//Clear buffer
	memset(screen->buffer, 0, screen->height*screen->width * sizeof(uint32_t));

	//Iterate through every pixel in the image
	for (int y = 0; y < SCREEN_HEIGHT; y++) {
		for (int x = 0; x < SCREEN_WIDTH; x++) {
			vec3 color(0);
			int count = 2;

			float step = 1/(float)(count*2+1);
			for (int i = -count; i <= count; i ++) {
				for (int j = -count; j <= count; j++) {
					color += RayTracerColor(x + i*step, y + j*step);
				}
			}
			color = color / (float)pow((count*2+1),2);
			//vec3 color = RayTracerColor(x, y);
			PutPixelSDL(screen, x, y, color);
		}
	}
}

vec3 RayTracerColor(float x, float y) {
	//Compute the corresponding ray direction
	vec4 rayDirection(x - SCREEN_WIDTH * 0.5, y - SCREEN_HEIGHT * 0.5, focalLength, 1);

	rayDirection = rayDirection * cameraRotMatrix;

	rayDirection = NormaliseVec(rayDirection);

	//Call the function ClosestIntersection to get the closest intersection in that direction
	Intersection closestIntersectionItem;
	if (ClosestIntersection(cameraPos, rayDirection, triangles, closestIntersectionItem)) {
		if (closestIntersectionItem.triangleIndex == -1) {
			vec3 sphereColor = CombineReflectionRefraction(closestIntersectionItem);
			return sphereColor;
		}
		else {
			vec3 color = triangles[closestIntersectionItem.triangleIndex].color * (DirectLight(closestIntersectionItem) + indirectLight);
			return color;
		}
	}
	else {
		//It should be black
		vec3 color(0.0, 0.0, 0.0);
		return color;
	}
}

mat3 RotMatrixX(float angle) {
	vec3 x0(1, 0, 0);
	vec3 x1(0, cosf(angle), sinf(angle));
	vec3 x2(0, -sinf(angle), cosf(angle));

	return mat3(x0, x1, x2);
}

mat3 RotMatrixY(float angle) {
	vec3 y1(cosf(angle), 0, -sinf(angle));
	vec3 y2(0, 1, 0);
	vec3 y3(sinf(angle), 0, cosf(angle));

	return mat3(y1, y2, y3);
}

void Rotate(mat3 rotation) {
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

bool Update()
{
	static int t = SDL_GetTicks();
	/* Compute frame time */
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;

	//std::cout << "Render time: " << dt << " ms." << std::endl;

	float angle = 0.01;
	float step = 0.2f;

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

					//Light Movement
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
				case SDLK_i:
					sphereCenter.z += step/10;
					break;
				case SDLK_k:
					sphereCenter.z -= step/10;
					break;
				/*case SDLK_l:
					cameraPos.x += step;
					break;
				case SDLK_j:
					cameraPos.x -= step;
					break;*/

				case SDLK_ESCAPE:
					/* Move camera quit */
					return false;
				}
			}
	}
	return true;
}

vec3 DirectLight(const Intersection& i) {

	vec3 updatedColor = vec3(0, 0, 0);

	//n is a unit vector describing the normal pointing out from the surface
	vec4 n = triangles[i.triangleIndex].normal;

	//r is a unit vector describing the direction from the surface point to the light source
	vec4 r = lightPos - i.position;

	//Normalise vectors n and r
	n = NormaliseVec(n);
	r = NormaliseVec(r);

	/*Check if the intersection given receives direct illumination*/

	//Calculate the distance between the intersection and the lightPos
	double d = Find3Distance(lightPos, vec3(i.position));

	//Calculate the closestIntersection of the intersection in the direction of the light
	Intersection closestIntersectionItem;

	if (ClosestIntersection(i.position + 0.0001f*n, r, triangles, closestIntersectionItem)) {
		if (closestIntersectionItem.distance>0.f && closestIntersectionItem.distance < d) {
			return updatedColor;
		}
	}

	double dotProduct = (n.x * r.x) + (n.y * r.y) + (n.z * r.z);

	if (dotProduct > 0) {

		//Calculate the number (0 < n < 1) with which the power P of the light (lightColor) will be multiplied with
		double percentage = dotProduct / (4 * M_PI * d * d);

		//Calculate the new value of light for the intersection
		updatedColor.x = lightColor.x * percentage;
		updatedColor.y = lightColor.y * percentage;
		updatedColor.z = lightColor.z * percentage;
	}
	return updatedColor;
}

vec3 CombineReflectionRefraction(const Intersection& intersection) {

	//Find vector from camera to intersection
	vec4 i = cameraPos - intersection.position;
	//Find norman of Intersection
	vec4 n = intersection.position - sphereCenter;

	//Normalise vectors ir and n ???maybe after magnitude???
	i = NormaliseVec(i);
	n = NormaliseVec(n);

	////Find angle of incidence
	double dotProduct = (i.x * n.x) + (i.y * n.y) + (i.z * n.z);

	//Figure out what is etai and what is etat
	if (dotProduct > 0) { //We go from air inside glass
		etai = air;
		etat = glass;
	}
	else { //We go from glass to air
		etai = glass;
		etat = air;
	}

	vec4 t = RefractionDirection(i, n, dotProduct);

	//Calculate angles for Fresnel
	double costheta1 = angle(i, n);
	double costheta2 = angle(n, -t);

	//Madness because Fresnel wasn't working
	double fresnelOrthogonal = costheta1;
	double fresnelParallel = costheta2;
	double fresnel = costheta2;
	
	//Call FresnelR to calculate the percentage that gets reflected
	FresnelR(costheta1, costheta2, fresnelOrthogonal, fresnelParallel, fresnel);

	vec3 reflectionColor = ReflectionColor(intersection, i, n, dotProduct);
	vec3 refractionColor = RefractionColor(intersection, i, n, dotProduct);

	//Return the combined color of reflection and refraction
	return ((float)fresnel*reflectionColor) + ((1 - (float)fresnel)*refractionColor);
	//return reflectionColor;
	//return refractionColor;
}

vec3 ReflectionColor(const Intersection& intersection, vec4 i, vec4 n, double dotProduct) {

	vec3 intersectionColor(1, 1, 1);

	float w = (float)(2 * dotProduct);

	//Find reflective ray
	vec4 r = i - (w * n);

	r = NormaliseVec(r);

	//Calculate the closestIntersection in the direction of reflection
	Intersection closestIntersectionItem;

	//Check if there is a SphereClosestIntersection with the reflection Ray
	if (sphereClosestIntersection(intersection.position + 0.001f*r, r, closestIntersectionItem)) {
		if (recursionCount < 5) {
			//If there is the closest intersection is with the shpere then we are inside the sphere are we want to have RECURSION
			recursionCount++;
			CombineReflectionRefraction(closestIntersectionItem);
		}
		else {
			//Zero contribution to the color
			return vec3(0, 0, 0);
		}
	}
	else if (ClosestIntersection(intersection.position + 0.001f*r, r, triangles, closestIntersectionItem)) {
		if (closestIntersectionItem.distance>0.f) {
			intersectionColor = triangles[closestIntersectionItem.triangleIndex].color;
			return intersectionColor;
		}
	}
	return intersectionColor;
}

vec3 RefractionColor(const Intersection& intersection, vec4 i, vec4 n, double dotProduct) {
	vec3 refractionColor(0, 0, 0);

	vec4 t = RefractionDirection(i, n, dotProduct);
	t = NormaliseVec(t);

	//Calculate the closestIntersection in the direction of refraction
	Intersection closestIntersectionItem;

	//If an intersection exists
	if (ClosestIntersection(intersection.position - 0.001f*t, -t, triangles, closestIntersectionItem)) {

		// --------------------------- 2 OPTIONS --------------------------------

		// --- 1: The ray just entered the sphere and we want to make a recursion, in which case the closest intersection is a sphere
		if (closestIntersectionItem.triangleIndex == -1) {
			CombineReflectionRefraction(closestIntersectionItem);
		}

		// --- 2: The rounds of recursion are complete and we want to exit the sphere
		else if (recursionCount > 5) {
			//Make it 100% refractive and exit the sphere to get a color back (a.k.a. just make another refraction round)

			//Find the point in the 3D space where the intersection with the other end of the sphere happened
			//Compute the normal
			vec4 normalTwo = closestIntersectionItem.position - sphereCenter;
			normalTwo = NormaliseVec(normalTwo);

			double dotProductTwo = (t.x * normalTwo.x) + (t.y * normalTwo.y) + (t.z * normalTwo.z);

			//Find the refraction direction of the exiting ray
			vec4 secondT = RefractionDirection(-t, normalTwo, dotProductTwo);

			//Send a ray towards this direction and return the color
			Intersection secondClosestIntersectionItem;

			if (ClosestIntersection(closestIntersectionItem.position + 0.01f*secondT, secondT, triangles, secondClosestIntersectionItem)) {
				if (closestIntersectionItem.triangleIndex == -1) {
					return vec3(0, 0, 0);
				}
				refractionColor = triangles[secondClosestIntersectionItem.triangleIndex].color;
				return refractionColor;
			}
		}
		//If the closest intersection is a triangle
		else if (closestIntersectionItem.distance>0.f) {
			refractionColor = triangles[closestIntersectionItem.triangleIndex].color;
			return refractionColor;
		}
	}
	//If no intersection is found
	return vec3(1, 1 , 1);//refractionColor;
}

vec4 RefractionDirection(const vec4 i, const vec4 n, double dotProduct) {

	vec4 nRefr = n;

	float c2 = sqrt(1 - (eta*eta)*(1 - (dotProduct*dotProduct)));

	//Calculate transmission ray
	vec4 t = (eta*i) + ((eta*(float)dotProduct) - c2)*n; // Not sure about this float casting

	return t;
}

void FresnelR(double costheta1, double costheta2, double fresnelOrthogonal, double fresnelParallel, double fresnel) {
	costheta1 = etat * costheta1;
	costheta2 = etai * costheta2;

	fresnelParallel = ((costheta1 - costheta2) / (costheta1 + costheta2)) * ((costheta1 - costheta2) / (costheta1 + costheta2));
	fresnelOrthogonal = ((costheta2 - costheta1) / (costheta2 + costheta1)) * ((costheta2 - costheta1) / (costheta2 + costheta1));

	fresnel = 0.5 * fresnelParallel * fresnelOrthogonal;
	
	return;
}

double Find3Distance(vec3 vectorOne, vec3 vectorTwo) {

	double sqrDifferenceX = (vectorOne.x - vectorTwo.x)*(vectorOne.x - vectorTwo.x);
	double sqrDifferenceY = (vectorOne.y - vectorTwo.y)*(vectorOne.y - vectorTwo.y);
	double sqrDifferenceZ = (vectorOne.z - vectorTwo.z)*(vectorOne.z - vectorTwo.z);
	double root = sqrt(sqrDifferenceX + sqrDifferenceY + sqrDifferenceZ);

	return root;
}

vec4 NormaliseVec(vec4 vec) {
	vec4 normalisedPoint;
	double length = sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);

	normalisedPoint.x = vec.x / length;
	normalisedPoint.y = vec.y / length;
	normalisedPoint.z = vec.z / length;
	normalisedPoint.w = vec.w;

	return normalisedPoint;
}

double angle(vec4 u, vec4 v) {

	//Calculate dot product
	double dotProduct = (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
	//Calculate magnitudes
	double uMagnitude = sqrt((u.x*u.x) + (u.y*u.y) + (u.z*u.z));
	double vMagnitude = sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z));
	//Calculate angle
	double cosangle = dotProduct / (uMagnitude * vMagnitude);

	return cosangle;
}