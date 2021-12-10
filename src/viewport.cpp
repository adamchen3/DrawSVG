#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan;
  double matrix[] = { 1.f / (2 * vspan), 0.f, (vspan - centerX) / (2 * vspan), 0.f, 1.f / (2 * vspan), (vspan - centerY) / (2 * vspan), 0.f, 0.f, 1.f };
  set_svg_2_norm(Matrix3x3(matrix));

}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
