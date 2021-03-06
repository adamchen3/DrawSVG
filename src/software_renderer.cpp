#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  //transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    transformation = svg_2_screen;
    draw_element(svg.elements[i]);
  }

  transformation = svg_2_screen;

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

  this->buffer_w = this->target_w * this->sample_rate;
  this->buffer_h = this->target_h * this->sample_rate;
  this->sample_buffer.resize(4 * this->buffer_w * this->buffer_h);
  memset(sample_buffer.data(), 255, 4 * buffer_w * buffer_h);
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

  this->buffer_w = this->target_w * this->sample_rate;
  this->buffer_h = this->target_h * this->sample_rate;
  this->sample_buffer.resize(4 * this->buffer_w * this->buffer_h);
  memset(sample_buffer.data(), 255, 4 * buffer_w * buffer_h);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack
  transformation = transformation * element->transform;
  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  Matrix3x3 groupTransformation = transformation;
  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    transformation = groupTransformation;
    draw_element(group.elements[i]);
  }

}

// Rasterization //

void SoftwareRendererImp::rasterize_sample_point(float x, float y, Color color) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= buffer_w) return;
  if (sy < 0 || sy >= buffer_h) return;

  // fill sample - NOT doing alpha blending!
  //sample_buffer[4 * (sx + sy * buffer_w)] = (uint8_t)(color.r * 255);
  //sample_buffer[4 * (sx + sy * buffer_w) + 1] = (uint8_t)(color.g * 255);
  //sample_buffer[4 * (sx + sy * buffer_w) + 2] = (uint8_t)(color.b * 255);
  //sample_buffer[4 * (sx + sy * buffer_w) + 3] = (uint8_t)(color.a * 255);

  // doing alpha blending using premultiplied alpha method
  float r = sample_buffer[4 * (sx + sy * buffer_w)] / 255.f;
  float g = sample_buffer[4 * (sx + sy * buffer_w) + 1] / 255.f;
  float b = sample_buffer[4 * (sx + sy * buffer_w) + 2] / 255.f;
  float a = sample_buffer[4 * (sx + sy * buffer_w) + 3] / 255.f;

  r *= a; g *= a; b *= a;
  color.r *= color.a; color.g *= color.a; color.b *= color.a;
  float ac = color.a + (1 - color.a) * a;
  float cr = color.r + (1 - color.a) * r;
  float cg = color.g + (1 - color.a) * g;
  float cb = color.b + (1 - color.a) * b;
  cr /= ac;
  cg /= ac;
  cb /= ac;
  sample_buffer[4 * (sx + sy * buffer_w)] = (uint8_t)(cr * 255);
  sample_buffer[4 * (sx + sy * buffer_w) + 1] = (uint8_t)(cg * 255);
  sample_buffer[4 * (sx + sy * buffer_w) + 2] = (uint8_t)(cb * 255);
  sample_buffer[4 * (sx + sy * buffer_w) + 3] = (uint8_t)(ac * 255);

}

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {
  // render in sample_buffer
  x *= sample_rate;
  y *= sample_rate;
  rasterize_sample_point(x, y, color);
  return;

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);

}

float point_to_line_dis(Vector2D p0, Vector2D p1, Vector2D p) {
  return (float) std::abs((p.x - p0.x) * (-p1.y + p0.y) + (p.y - p0.y) * (p1.x - p0.x)) / (p1 - p0).norm();
}

float sign(float x0, float y0, float x1, float y1, float x2, float y2)
{
  return (x0 - x2) * (y1 - y2) - (x1 - x2) * (y0 - y2);
}

bool PointInTriangle(float x0, float y0, float x1, float y1, float x2, float y2, float x, float y)
{
  float d1, d2, d3;
  bool has_neg, has_pos;

  d1 = sign(x, y, x0, y0, x1, y1);
  d2 = sign(x, y, x1, y1, x2, y2);
  d3 = sign(x, y, x2, y2, x0, y0);

  has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  return !(has_neg && has_pos);
}

bool coverage(float x0, float y0, float x1, float y1, float x2, float y2, float x, float y) {
  return PointInTriangle(x0, y0, x1, y1, x2, y2, x, y);

  // the below method has bug
  //float w1 = (x0 * (y2 - y0) + (y - y0) * (x2 - x0) - x * (y2 - y0)) / ((y1 - y0) * (x2 - x0) - (x1 - x0) * (y2 - y0));
  //float w2 = (y - y0 - w1 * (y1 - y0)) / (y2 - y0);
  //return w1 >= 0 && w2 >= 0 && (w1 + w2) <= 1;
}

void SoftwareRendererImp::draw_triangle(Vector2D p0, Vector2D p1, Vector2D t0, Vector2D t1, Vector2D t2, Color color, int width = 2) {
  float minx = std::min({ t0.x, t1.x, t2.x });
  float maxx = std::max({ t0.x, t1.x, t2.x });
  float miny = std::min({ t0.y, t1.y, t2.y });
  float maxy = std::max({ t0.y, t1.y, t2.y });
  int xstart = (int)floor(minx);
  int xend = (int)ceil(maxx);
  int ystart = (int)floor(miny);
  int yend = (int)ceil(maxy);
  for (int x = xstart; x < xend; x++) {
    for (int y = ystart; y < yend; y++) {
      float dis = point_to_line_dis(p0, p1, { x + 0.5f, y + 0.5f });
      float factor = dis * 2 / width;
      if (coverage(t0.x, t0.y, t1.x, t1.y, t2.x, t2.y, x + 0.5f, y + 0.5f)) {
        rasterize_point(x, y, color * (1 - factor));
      }
    }
  }
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {
  // Task 2: 
  // Implement line rasterization

  // if render in sample_buffer, coordinate must multiple by sample_rate
  x0 *= sample_rate;
  y0 *= sample_rate;
  x1 *= sample_rate;
  y1 *= sample_rate;

  // Draw line as a triangle
  // ?????????????????????
  //if (x0 == x1) {
  //  draw_triangle({ x0, y0 }, { x1, y1 }, { x0 - 1, y0 }, { x0 + 1, y0 }, { x1 + 1, y1 }, color);
  //  draw_triangle({ x0, y0 }, { x1, y1 }, { x1 + 1, y1 }, { x1 - 1, y1 }, { x0 - 1, y0 }, color);
  //}
  //else if (y0 == y1) {
  //  draw_triangle({ x0, y0 }, { x1, y1 }, { x0, y0 - 1 }, { x1, y1 - 1 }, { x1, y1 + 1 }, color);
  //  draw_triangle({ x0, y0 }, { x1, y1 }, { x1, y1 + 1 }, { x0, y0 + 1 }, { x0, y0 - 1 }, color);
  //}
  //else {
  //  Vector2D v0{ -1, 0 };
  //  Vector2D v1{x1 - x0, y1 - y0};
  //  float offsetY = (float)std::abs(dot(v0, v1) / v1.norm());
  //  float slope = y1 - y0 / x1 - x0;
  //  Vector2D v2{ 1, -1 / slope };
  //  Vector2D v3{ 1, 0 };
  //  float offsetX = (float)std::abs(dot(v2, v3) / v2.norm());
  //  if (slope > 0) {
  //    offsetY = -offsetY;
  //  }
  //  draw_triangle({ x0, y0 }, { x1, y1 }, { x0 - offsetX, y0 - offsetY }, { x0 + offsetX, y0 + offsetY  }, { x1 + offsetX, y1 + offsetY }, color);
  //  draw_triangle({ x0, y0 }, { x1, y1 }, { x1 + offsetX, y1 + offsetY }, { x1 - offsetX, y1 - offsetY }, { x0 - offsetX , y0 - offsetY }, color);
  //}


  int  x_start, x_end, y_start;
  bool steep = false; // ????????????y???????????????
  float dy = abs(y0 - y1);
  float dx = abs(x0 - x1);
  if (dy > dx) {
    std::swap(x0, y0);
  std:swap(x1, y1);
    steep = true;
  }

  if (x0 > x1) {
    std::swap(x0, x1);
    std::swap(y0, y1);
  }

  /** Bresenham start **/
 /* dy = y1 - y0;
  dx = x1 - x0;
  x_start = (int)round(x0);
  x_end = (int)round(x1);
  y_start = (int)round(y0);

  float error = 0;
  for (int i = x_start; i <= x_end; i++) {
      if (steep) {
          rasterize_point(y_start, i, color);
      }
      else {
          rasterize_point(i, y_start, color);
      }
      error += 2 * abs(dy);
      if (abs(error) >= dx) {
          y_start += (dy > 0 ? 1 : -1);
          error -= 2 * dx;
      }
  }
  return;*/
  /** Bresenham end **/

  /** Xiaolin Wu's line algorithm */
  dy = y1 - y0;
  dx = x1 - x0;
  x_start = (int)round(x0);
  x_end = (int)round(x1);
  float slope = dy / dx;
  float y_f = y0 + slope * (x_start - x0);

  for (int i = x_start; i <= x_end; i++) {
    float y_start = (int)floor(y_f);
    float y_l = y_f + 0.5f * slope;
    float y_c = y_start + 0.5f;
    float dis = abs(y_c - y_l);
    float dis2;
    int index;
    if (abs(y_c + 1 - y_l) < abs(y_c - 1 - y_l)) {
      index = 1;
      dis2 = abs(y_c + 1 - y_l);
    }
    else {
      index = -1;
      dis2 = abs(y_c - 1 - y_l);
    }
    if (steep) {
      rasterize_sample_point(y_start, i, color * (1 - dis));
      rasterize_sample_point(y_start + index, i, color * (1 - dis2));
    }
    else {
      rasterize_sample_point(i, y_start, color * (1 - dis));
      rasterize_sample_point(i, y_start + index, color * (1 - dis2));
    }
    y_f += slope;
  }
}

bool line_segment_interset(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3) {
  float ax = x1 - x0;
  float ay = y1 - y0;
  float bx = x3 - x2;
  float by = y3 - y2;
  float cx = x0 - x2;
  float cy = y0 - y2;
  float denominator = ay * bx - ax * by;
  float numerator = by * cx - bx * cy;
  if (denominator == 0) {
    return true;
  }
  else if (denominator < 0) {
    if (numerator <= 0 && numerator >= denominator) {
      return true;
    }
  }
  else if (denominator > 0) {
    if (numerator >= 0 && numerator <= denominator) {
      return true;
    }
  }
  numerator = ax * cy - ay * cx;
  if (denominator < 0) {
    if (numerator <= 0 && numerator >= denominator) {
      return true;
    }
  }
  else if (denominator > 0) {
    if (numerator >= 0 && numerator <= denominator) {
      return true;
    }
  }
  return false;
}

void SoftwareRendererImp::rasterize_triangle(float x0, float y0,
  float x1, float y1,
  float x2, float y2,
  Color color) {
  // Task 3: 
  // Implement triangle rasterization
  
  x0 *= sample_rate;
  y0 *= sample_rate;
  x1 *= sample_rate;
  y1 *= sample_rate;
  x2 *= sample_rate;
  y2 *= sample_rate;

  // draw triangle outline
  //rasterize_line(x0, y0, x1, y1, color);
  //rasterize_line(x1, y1, x2, y2, color);
  //rasterize_line(x2, y2, x0, y0, color);
  //return;

  float minx = std::min({ x0, x1, x2 });
  float maxx = std::max({ x0, x1, x2 });
  float miny = std::min({ y0, y1, y2 });
  float maxy = std::max({ y0, y1, y2 });
  int xstart = (int)floor(minx);
  int xend = (int)ceil(maxx);
  int ystart = (int)floor(miny);
  int yend = (int)ceil(maxy);

  // using all bounding box
  //for (int x = xstart; x < xend; x++) {
  //  for (int y = ystart; y < yend; y++) {
  //    if (coverage(x0, y0, x1, y1, x2, y2, x + 0.5f, y + 0.5f)) {
  //      rasterize_sample_point(x, y, color);
  //    }
  //  }
  //}
  //return;

  // using blockwise method
  const int blocksize = 16;
  for (int offsetX = 0; xstart + offsetX < xend; offsetX += blocksize) {
    for (int offsetY = 0; ystart + offsetY < yend; offsetY += blocksize) {
      int x_lt = xstart + offsetX;
      int y_lt = ystart + offsetY;
      int x_rt = std::min(x_lt + blocksize, xend);
      int y_rt = ystart + offsetY;
      int x_lb = xstart + offsetX;
      int y_lb = std::min(y_lt + blocksize, yend);
      int x_rb = x_rt;
      int y_rb = y_lb;

      if (line_segment_interset(x_lt, y_lt, x_rt, y_rt, x0, y0, x1, y1)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lt, y_lt, x_rt, y_rt, x1, y1, x2, y2)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lt, y_lt, x_rt, y_rt, x2, y2, x0, y0)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lt, y_lt, x_lb, y_lb, x0, y0, x1, y1)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lt, y_lt, x_lb, y_lb, x1, y1, x2, y2)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lt, y_lt, x_lb, y_lb, x2, y2, x0, y0)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_rt, y_rt, x_rb, y_rb, x0, y0, x1, y1)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_rt, y_rt, x_rb, y_rb, x1, y1, x2, y2)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_rt, y_rt, x_rb, y_rb, x2, y2, x0, y0)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lb, y_lb, x_rb, y_rb, x0, y0, x1, y1)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lb, y_lb, x_rb, y_rb, x1, y1, x2, y2)) {
        goto check_and_draw_tri;
      }

      if (line_segment_interset(x_lb, y_lb, x_rb, y_rb, x2, y2, x0, y0)) {
        goto check_and_draw_tri;
      }

      if (coverage(x0, y0, x1, y1, x2, y2, x_lt, y_lt)) {
        for (int x = x_lt; x < x_rt; x++) {
          for (int y = y_lt; y < y_lb; y++) {
            rasterize_sample_point(x, y, color);
          }
        }
      }

      continue;

check_and_draw_tri:
      for (int x = x_lt; x < x_rt; x++) {
        for (int y = y_lt; y < y_lb; y++) {
          if (coverage(x0, y0, x1, y1, x2, y2, x + 0.5f, y + 0.5f)) {
            rasterize_sample_point(x, y, color);
          }
        }
      }
    }
  }
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization
  int xstart = (int)floor(x0 + 0.5f);
  int ystart = (int)floor(y0 + 0.5f);
  int xend = (int)round(x1 + 0.5f);
  int yend = (int)round(y1 + 0.5f);
  int width = xend - xstart;
  int height = yend - ystart;
  float scale = max(width, height) * 1.f / max(tex.width, tex.height);

  for (int x = xstart; x < xend; x++) {
    for (int y = ystart; y < yend; y++) {
      float u = ((x - xstart) * 1.f + .5f) / width;
      float v = ((y - ystart) * 1.f + .5f) / height;
      //rasterize_point(x, y, sampler->sample_nearest(tex, u, v, 0));
      //rasterize_point(x, y, sampler->sample_bilinear(tex, u, v, 0));
      rasterize_point(x, y, sampler->sample_trilinear(tex, u, v, scale, scale));
    }
  }
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  memset(render_target, 255, 4 * target_w * target_h);
  if (sample_rate == 1) {
    memcpy(render_target, sample_buffer.data(), 4 * target_w * target_h);
  }
  else {
    for (int x = 0; x < target_w; x++) {
      for (int y = 0; y < target_h; y++) {
        unsigned int r = 0, g = 0, b = 0, a = 0;
        for (int i = x * sample_rate; i < (x + 1) * sample_rate; i++){
          for (int j = y * sample_rate; j < (y + 1) * sample_rate; j++) {
            r += sample_buffer[4 * (i + j * buffer_w) + 0];
            g += sample_buffer[4 * (i + j * buffer_w) + 1];
            b += sample_buffer[4 * (i + j * buffer_w) + 2];
            a += sample_buffer[4 * (i + j * buffer_w) + 3];
          }
        }
        render_target[4 * (x + y * target_w) + 0] = (uint8_t)(r / (sample_rate * sample_rate));
        render_target[4 * (x + y * target_w) + 1] = (uint8_t)(g / (sample_rate * sample_rate));
        render_target[4 * (x + y * target_w) + 2] = (uint8_t)(b / (sample_rate * sample_rate));
        render_target[4 * (x + y * target_w) + 3] = (uint8_t)(a / (sample_rate * sample_rate));
      }
    }
  }
  memset(sample_buffer.data(), 255, 4 * buffer_w * buffer_h);
  return;
}


} // namespace CMU462
