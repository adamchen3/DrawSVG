#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
  //Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    //Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];
    MipLevel& upmip = tex.mipmap[i - 1];

    //for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      //float_to_uint8( &mip.texels[i], &c.r );
    //}

    for (int i = 0; i < mip.width; i++) {
      for (int j = 0; j < mip.height; j++) {

        int r{ 0 }, g{ 0 }, b{ 0 }, a{ 0 };
        for (int upi = i * 2; upi < i * 2 + 2; upi++) {
          for (int upj = j * 2; upj < j * 2 + 2; upj++) {
            int upindex = upi * upmip.width + upj;
            r += upmip.texels[upindex * 4];
            g += upmip.texels[upindex * 4 + 1];
            b += upmip.texels[upindex * 4 + 2];
            a += upmip.texels[upindex * 4 + 3];
          }
        }
        r /= 4; g /= 4; b /= 4; a /= 4;

        int index = i * mip.width + j;
        mip.texels[index * 4] = r;
        mip.texels[index * 4 + 1] = g;
        mip.texels[index * 4 + 2] = b;
        mip.texels[index * 4 + 3] = a;
      }
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
  MipLevel mipmap = tex.mipmap[level];
  int i = (int)floor(u * mipmap.width);
  int j = (int)floor(v * mipmap.height);
  int index = i * mipmap.width + j;
  Color c;
  c.r = mipmap.texels[index * 4] / 255.f;
  c.g = mipmap.texels[index * 4 + 1] / 255.f;
  c.b = mipmap.texels[index * 4 + 2] / 255.f;
  c.a = mipmap.texels[index * 4 + 3] / 255.f;
  return c;
  
  // return magenta for invalid level
  return Color(1,0,1,1);

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering
  MipLevel mipmap = tex.mipmap[level];
  if (level == tex.mipmap.size() - 1) {
    Color c;
    c.r = tex.mipmap[level].texels[0] / 255.f;
    c.g = tex.mipmap[level].texels[1] / 255.f;
    c.b = tex.mipmap[level].texels[2] / 255.f;
    c.a = tex.mipmap[level].texels[3] / 255.f;
    return c;
  }

  int i = (int)floor(u * mipmap.width - 0.5f);
  int j = (int)floor(v * mipmap.height - 0.5f);
  float s = u * mipmap.width - (i + 0.5f);
  float t = v * mipmap.height - (j + 0.5f);

  int a = i, b = j;
  if (a < 0) {
    a = mipmap.width - 1;
  }
  else if (a >= mipmap.width) {
    a = 0;
  }
  if (b < 0) {
    b = mipmap.height - 1;
  }
  else if (b >= mipmap.height) {
    b = 0;
  }
  int index00 = a * mipmap.width + b;

  a = i + 1; b = j;
  if (a < 0) {
    a = mipmap.width - 1;
  }
  else if (a >= mipmap.width) {
    a = 0;
  }
  if (b < 0) {
    b = mipmap.height - 1;
  }
  else if (b >= mipmap.height) {
    b = 0;
  }
  int index01 = a * mipmap.width + b;

  a = i; b = j + 1;
  if (a < 0) {
    a = mipmap.width - 1;
  }
  else if (a >= mipmap.width) {
    a = 0;
  }
  if (b < 0) {
    b = mipmap.height - 1;
  }
  else if (b >= mipmap.height) {
    b = 0;
  }
  int index10 = a * mipmap.width + b;

  a = i + 1; b = j + 1;
  if (a < 0) {
    a = mipmap.width - 1;
  }
  else if (a >= mipmap.width) {
    a = 0;
  }
  if (b < 0) {
    b = mipmap.height - 1;
  }
  else if (b >= mipmap.height) {
    b = 0;
  }
  int index11 = a * mipmap.width + b;

  Color c00;
  c00.r = mipmap.texels[index00 * 4] / 255.f;
  c00.g = mipmap.texels[index00 * 4 + 1] / 255.f;
  c00.b = mipmap.texels[index00 * 4 + 2] / 255.f;
  c00.a = mipmap.texels[index00 * 4 + 3] / 255.f;
  Color c10;
  c10.r = mipmap.texels[index10 * 4] / 255.f;
  c10.g = mipmap.texels[index10 * 4 + 1] / 255.f;
  c10.b = mipmap.texels[index10 * 4 + 2] / 255.f;
  c10.a = mipmap.texels[index10 * 4 + 3] / 255.f;
  Color c01;
  c01.r = mipmap.texels[index01 * 4] / 255.f;
  c01.g = mipmap.texels[index01 * 4 + 1] / 255.f;
  c01.b = mipmap.texels[index01 * 4 + 2] / 255.f;
  c01.a = mipmap.texels[index01 * 4 + 3] / 255.f;
  Color c11;
  c11.r = mipmap.texels[index11 * 4] / 255.f;
  c11.g = mipmap.texels[index11 * 4 + 1] / 255.f;
  c11.b = mipmap.texels[index11 * 4 + 2] / 255.f;
  c11.a = mipmap.texels[index11 * 4 + 3] / 255.f;

  // PS:别删除下面注释最后的那个.，不然编译不过。傻逼编译器
  // 系数s和t互换一下，效果也差太远了吧，手动捂脸.
  Color c0;
  c0 = (1 - t) * c00 + t * c10;
  Color c1;
  c1 = (1 - t) * c01 + t * c11;

  Color c;
  c = (1 - s) * c0 + s * c1;
  return c;

  // return magenta for invalid level
  return Color(1,0,1,1);

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering
  float scale = 1 / u_scale;
  float s = max(0.f, log2f(scale));
  int level = (int)floor(s);
  Color c0 = sample_bilinear(tex, u, v, level);
  Color c1;
  if (level >= tex.mipmap.size() - 1) {
    c1 = sample_bilinear(tex, u, v, tex.mipmap.size() - 1);
  }
  else {
    c1 = sample_bilinear(tex, u, v, level + 1);
  }
  Color c;
  s -= level;
  c = (1 - s) * c0 + s * c1;
  return c;

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CMU462
