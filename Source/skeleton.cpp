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

vector<Triangle> triangles;
struct Intersection
{
  vec4 position;
  float distance;
  int triangleIndex;
};
float focalLength = 256;
vec4 cameraPos(0, 0, -3, 1);

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
  //Max size of float
  float m = std::numeric_limits<float>::max();
  float smallest_d = m;
  bool intersection_found = false;
  //Find size of triangles vector
  int tri_size = triangles.size();
  //For each triangle find an intersection point x
  for (int i = 0; i < tri_size; i++) {
    vec4 v0 = triangles[i].v0;
    vec4 v1 = triangles[i].v1;
    vec4 v2 = triangles[i].v2;

    vec3 e1 = vec3(v1 - v0);
    vec3 e2 = vec3(v2 - v0);
    vec3 b = vec3(start - v0);

    mat3 A( -vec3(dir), e1, e2 );
    if (glm::determinant(A) != 0) {
      vec3 x = glm::inverse( A ) * b;

      float t = x.x;
      float u = x.y;
      float v = x.z;

      if (t < smallest_d) {
        //Check intersection conditions satisfied
        if (0 <= u && 0 <= v && u + v <= 1 && 0 <= t) {
          //Set values for the intersection
          intersection_found = true;
          closestIntersection.position = start + t * dir;
          closestIntersection.distance = t;
          closestIntersection.triangleIndex = i;
        }
      }
    }
  }
  return intersection_found;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  //Iterate over each pixel in the image
  for (int x = 0; x < SCREEN_WIDTH; x++) {
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
      //Compute ray direction & start point
      vec4 s = cameraPos;
      vec4 d = vec4(x - (SCREEN_WIDTH/2), y - (SCREEN_HEIGHT/2), focalLength, 1);
      //Get closest intersection
      Intersection closest = {
        .position = vec4(0, 0, 0, 0),
        .distance = 0,
        .triangleIndex = 0
      };
      if (ClosestIntersection(s, d, triangles, closest)) {
        //Set colour of pixel to colour of intersected triangle
        vec3 col = triangles[closest.triangleIndex].color;
        PutPixelSDL(screen, x, y, col);
      }
    }
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
