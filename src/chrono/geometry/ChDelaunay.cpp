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

#include "chrono/geometry/ChDelaunay3D.h"
#include <algorithm>
#include <cmath>

namespace chrono {

bool ChDelaunayFace::operator==(const ChDelaunayFace& other) const {
    // Check all combinations to see if the faces share the exact same 3 vertices
    return (p1 == other.p1 || p1 == other.p2 || p1 == other.p3) &&
           (p2 == other.p1 || p2 == other.p2 || p2 == other.p3) && (p3 == other.p1 || p3 == other.p2 || p3 == other.p3);
}

bool ChDelaunayTetrahedron::CircumsphereContains(const ChVector3d& pt) const {
    // Treat a as the center point of the tetrahedron
    ChVector3d r1 = b - a;
    ChVector3d r2 = c - a;
    ChVector3d r3 = d - a;

    // Calculate the circumcenter relative to 'a'
    double det = r1.x() * (r2.y() * r3.z() - r2.z() * r3.y()) - r1.y() * (r2.x() * r3.z() - r2.z() * r3.x()) +
                 r1.z() * (r2.x() * r3.y() - r2.y() * r3.x());

    double val1 = r1.Length2();
    double val2 = r2.Length2();
    double val3 = r3.Length2();

    ChVector3d cross_r1_r2 = r1.Cross(r2);  
    ChVector3d cross_r2_r3 = r2.Cross(r3);
    ChVector3d cross_r3_r1 = r3.Cross(r1);

    if (std::abs(det) < 1e-12) {
        return false; 
    }

    ChVector3d center_offset = (cross_r2_r3 * val1 + cross_r3_r1 * val2 + cross_r1_r2 * val3) / (2.0 * det);
    ChVector3d circumcenter = a + center_offset;

    double radius_squared = center_offset.Length2();
    double distance_squared = (pt - circumcenter).Length2();

    return distance_squared <= radius_squared;
}

bool ChDelaunayTetrahedron::SharesVertexWith(const ChDelaunayTetrahedron& other) const {
    const ChVector3d* superTetrahedraVertices[] = {&a, &b, &c, &d};
    const ChVector3d* otherTetrahedronVertices[] = {&other.a, &other.b, &other.c, &other.d};

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (*superTetrahedraVertices[i] == *otherTetrahedronVertices[j])
                // Returns true if the super tetrahedron vertice is the same as another tetrahedron
                return true;
        }
    }
    return false;
}

std::vector<ChDelaunayTetrahedron> ChDelaunay3D::CreateTetrahedralization(const std::vector<ChVector3d>& points) {
    if (points.size() < 4) {
        return {};
    }
    std::vector<ChDelaunayTetrahedron> tetrahedra;

    // Calculate bounding box
    double minX = points[0].x(), minY = points[0].y(), minZ = points[0].z();
    double maxX = minX, maxY = minY, maxZ = minZ;

    // find min x, y, and z values
    for (const auto& pt : points) {
        minX = std::min(minX, pt.x());
        minY = std::min(minY, pt.y());
        minZ = std::min(minZ, pt.z());
        maxX = std::max(maxX, pt.x());
        maxY = std::max(maxY, pt.y());
        maxZ = std::max(maxZ, pt.z());
    }

    double deltaMax = std::max({maxX - minX, maxY - minY, maxZ - minZ});
    double midX = (minX + maxX) / 2.0;
    double midY = (minY + maxY) / 2.0;
    double midZ = (minZ + maxZ) / 2.0;

    // Setup the encompassing tetrahedron 
    ChVector3d p1(midX, midY + 5 * deltaMax, midZ);
    ChVector3d p2(midX - 5 * deltaMax, midY - 5 * deltaMax, midZ - 5 * deltaMax);
    ChVector3d p3(midX + 5 * deltaMax, midY - 5 * deltaMax, midZ - 5 * deltaMax);
    ChVector3d p4(midX, midY - 5 * deltaMax, midZ + 5 * deltaMax);

    ChDelaunayTetrahedron superTetrahedra = {p1, p2, p3, p4};
    tetrahedra.push_back(superTetrahedra);

 
    for (const auto& pt : points) {
        std::vector<ChDelaunayTetrahedron> badTetrahedra;
        std::vector<ChDelaunayFace> cavityFaces;

        // Find bad tetrahedra
        for (const auto& tet : tetrahedra) {
            if (tet.CircumsphereContains(pt)) {
                badTetrahedra.push_back(tet);
            }
        }

        // Find the boundary faces of the 3D cavity
        for (size_t i = 0; i < badTetrahedra.size(); ++i) {
            ChDelaunayFace faces[4] = {{badTetrahedra[i].a, badTetrahedra[i].b, badTetrahedra[i].c},
                                       {badTetrahedra[i].a, badTetrahedra[i].b, badTetrahedra[i].d},
                                       {badTetrahedra[i].a, badTetrahedra[i].c, badTetrahedra[i].d},
                                       {badTetrahedra[i].b, badTetrahedra[i].c, badTetrahedra[i].d}};

            // Check if a face is shared by two bad tetrahedra
            for (int j = 0; j < 4; ++j) {
                bool isShared = false;
                for (size_t k = 0; k < badTetrahedra.size(); ++k) {
                    if (i == k)
                        continue;
                    ChDelaunayFace otherFaces[4] = {{badTetrahedra[k].a, badTetrahedra[k].b, badTetrahedra[k].c},
                                                    {badTetrahedra[k].a, badTetrahedra[k].b, badTetrahedra[k].d},
                                                    {badTetrahedra[k].a, badTetrahedra[k].c, badTetrahedra[k].d},
                                                    {badTetrahedra[k].b, badTetrahedra[k].c, badTetrahedra[k].d}};
                    if (faces[j] == otherFaces[0] || faces[j] == otherFaces[1] || faces[j] == otherFaces[2] ||
                        faces[j] == otherFaces[3]) {
                        isShared = true;
                        break;
                    }
                }
                if (!isShared) {
                    cavityFaces.push_back(faces[j]);
                }
            }
        }

        // Remove bad tetrahedra from the main list
        tetrahedra.erase(std::remove_if(tetrahedra.begin(), tetrahedra.end(),
                                        [&badTetrahedra](const ChDelaunayTetrahedron& t) {
                                            for (const auto& bt : badTetrahedra) {
                                                if (t.a == bt.a && t.b == bt.b && t.c == bt.c && t.d == bt.d)
                                                    return true;
                                            }
                                            return false;
                                        }),
                         tetrahedra.end());

        // Connect the new point to the cavity's boundary faces to form new tetrahedra
        for (const auto& face : cavityFaces) {
            tetrahedra.push_back({face.p1, face.p2, face.p3, pt});
        }
    }

    // Remove all temporary tetrahedra that touch the vertices of the encompassing Tetrahedron
    tetrahedra.erase(std::remove_if(tetrahedra.begin(), tetrahedra.end(),
        [&superTetrahedra](const ChDelaunayTetrahedron& t) {
            return t.SharesVertexWith(superTetrahedra);
        }),  
    tetrahedra.end());

    return tetrahedra;
}


}  

