#ifndef CSE168_TEX_GEN
#define CSE168_TEX_GEN

#include "Vector3.h"
#include "Texture.h"
#include "Perlin.h"
/**
 * Procedurally makes a stone texture using cellular texturing
 */
class TextureGenerator {
public:
   static Texture* generateCellTexture(int width, int height, int numCells);
   static Texture* generateStoneTexture(int width, int height, int numCells);

};



#endif