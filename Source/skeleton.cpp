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

struct Intersection
{
  vec4 position;
  float distance;
  int triangleIndex;
};
float focalLength = 1;
vec4 cameraPos(0, 0, 0, 1);

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen);
bool ClosestIntersection ( vec4 start,
                           vec4 dir,
                           const vector<Triangle>& triangles,
                           Intersection& closestIntersection );

int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  vector<Triangle> triangles;
  LoadTestModel( triangles );

  while( NoQuitMessageSDL() )
    {
      Update();
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

//Finds the closest intersection between a ray and a triangle
bool ClosestIntersection ( vec4 start,
                           vec4 dir,
                           const vector<Triangle>& triangles,
                           Intersection& closestIntersection )
{
  //Find size of triangles vector
  int tri_size = triangles.size();
  //For each triangle finx intersection point x
  for (int i = 0; i < tri_size; i++) {
    vec4 v0 = triangles[i].v0;
    vec4 v1 = triangles[i].v1;
    vec4 v2 = triangle[i].v2;

    vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z)
    vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z)
    vec3 b = vec3(s.x-v0.x,s.y-v0.y,s.z-v0.z);

    mat3 A( -d, e1, e2 );
    vec3 x = glm::inverse( A ) * b;

    float t = x.x;
    float u = x.y;
    float v = x.z;

    //Check intersection conditions satisfied
    if (0 < u && 0 < v && u + v < 1 && 0 <= t) {
      //Set values for the intersection
      closestIntersection.position = i;
      closestIntersection.distance = i;
      closestIntersection.triangleIndex = i;
      return true;
    }
    else {
      return false;
    }
  }
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 colour(1.0,0.0,0.0);
  for(int i=0; i<1000; i++)
    {
      uint32_t x = rand() % screen->width;
      uint32_t y = rand() % screen->height;
      PutPixelSDL(screen, x, y, colour);
    }
}

/*Place updates of parameters here*/
void Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
}
