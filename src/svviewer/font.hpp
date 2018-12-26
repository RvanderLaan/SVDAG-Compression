
//=============================================================
//  This file is part of the SymVox (Symmetry Voxelization) software
//  Copyright (C) 2016 by CRS4 Visual Computing Group, Pula, Italy
//
//  For more information, visit the CRS4 Visual Computing Group 
//  web pages at http://vic.crs4.it
//
//  This file may be used under the terms of the GNU General Public
//  License as published by the Free Software Foundation and appearing
//  in the file LICENSE included in the packaging of this file.
//
//  CRS4 reserves all rights not expressly granted herein.
//  
//  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
//  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
//  FOR A PARTICULAR PURPOSE.
//=============================================================

#pragma once

#include <utility>

namespace vic {
  
  namespace gl {
    
    class font_glyph {
    protected:
      float u_, v_;
      float width_;
      float left_margin_;
    public:
      font_glyph(float u=0.0f, float v=0.0f, float w=1.0f, float l=0.0f)
          : u_(u), v_(v), width_(w), left_margin_(l) {
      }
      
      inline float u() const {return u_; }
      inline float v() const {return v_; }
      inline float width() const {return width_; }
      inline float left_margin() const { return left_margin_; }
    };
    
    class font {
    protected:
      unsigned char *texture_data_;
      int texture_dx_;
      int texture_dy_;
      int character_dx_;
      int character_dy_;
      float character_fixed_width_;
      float character_fixed_height_;
      bool  is_proportional_;
      float scale_x_;
      float scale_y_;
      float spacing_;
      font_glyph glyphs_[256];
    protected:
      mutable unsigned int texture_id_;
      mutable int begin_depth_;
    protected:
      void construct_from_one_channel(const unsigned char* tiled,
                                      int sz,
                                      bool proportional,
                                      bool shadowed,
                                      bool embossed);
                                   
      void reset_texture_data(const unsigned char* idata,bool shadowed,bool embossed);
      void reset_glyphs();
    public:

      /**
       *  Create a default font from an embedded image
       *  spec.
       */
      font();
      
      /**
       *  Create a font from a tiled 8bit image. The image
       *  must be a square of sz x sz pixels, with sz=2^n. It must be
       *  subdivided into 16x16 tiles which contain all the
       *  characters of the ASCII mapping.
       */
      font(const unsigned char* tiled,
           int sz,
           bool proportional = true);

      virtual ~font();

      float scale_x() const { return scale_x_; }
      float scale_y() const { return scale_y_; }
      bool is_proportional() const { return is_proportional_; }
      float spacing() const { return spacing_; }
      
      void set_spacing(float spacing);
      void set_scaling(float sx, float sy);
      void set_proportional(bool b);

      float string_width(const char* s) const;

      void string_bbox(const char* s, float& x0, float&y0, float&x1, float& y1) const;
      
      void begin() const;
      void end() const;

      void cleanup() const;
      
      void write(const char* s) const;
      std::pair<int, int> render_size(const char* s) const;

    };
  }
}
