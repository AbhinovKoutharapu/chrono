// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Abhinov Koutharapu
// =============================================================================

#ifndef CH_DELAUNAY_3D_H
#define CH_DELAUNAY_3D_H

#include <vector>
#include "chrono/core/ChVector3.h"

namespace chrono {

/// A triangular face of a 3D tetrahedron.
struct ChDelaunayFace {
    ChVector3d p1;
    ChVector3d p2;
    ChVector3d p3;

    /// Returns true if both faces share the exact same three vertices.
    bool operator==(const ChDelaunayFace& other) const;
};

/// A 3D tetrahedron defined by four vertices.
struct ChDelaunayTetrahedron {
    ChVector3d a;
    ChVector3d b;
    ChVector3d c;
    ChVector3d d;

    /// Returns true if the given point lies inside this tetrahedron's circumsphere.
    bool CircumsphereContains(const ChVector3d& pt) const;

    /// Returns true if this tetrahedron shares one or more vertices with another tetrahedron.
    bool SharesVertexWith(const ChDelaunayTetrahedron& other) const;
};

/// Creates the 3D Delaunay tetrahedral meshes.
class ChDelaunay3D {
  public:
    /// Builds a 3D Delaunay mesh from a point cloud using the Boyer-Watson algorithm.
    static std::vector<ChDelaunayTetrahedron> CreateTetrahedralization(const std::vector<ChVector3d>& points);
};

} 

#endif  