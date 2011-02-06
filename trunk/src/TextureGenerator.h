#ifndef CSE168_TEX_GEN
#define CSE168_TEX_GEN

#include "Vector3.h"
#include "Texture.h"
#include "Perlin.h"
/**
 * Procedurally makes a stone texture using cellular texturing
 */
const class TextureGenerator {
public:
   Texture* const generateCellTexture(int width, int height, int numCells);
   Texture* const generateStoneTexture(int width, int height, int numCells);

};



#endif