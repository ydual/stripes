// -----------------------------------------------------------------------------
// libDDG -- Mesh.h
// -----------------------------------------------------------------------------
//
// Mesh represents a polygonal surface mesh using the halfedge data structure.
// It is essentially a large collection of disjoint vertices, edges, and faces
// that are ``glued together'' by halfedges which encode connectivity (see
// the documentation for an illustration).  By construction, the halfedge data
// structure cannot represent nonorientable surfaces or meshes with nonmanifold
// edges.
//
// Mesh elements are referenced using iterators -- common usage of these
// iterators is to either traverse an entire vector of mesh elements:
//
//    // visit all vertices
//    for( VertexIter i = vertices.begin(); i != vertices.end(); i++ )
//    {
//       //...
//    }
//
// or to perform a local traversal over the neighborhood of some mesh element:
//
//    // visit both halfedges of edge e
//    HalfEdgeIter he = e->he;
//    do
//    {
//       // ...
//
//       he = he->flip;
//    }
//    while( he != e->he );
//
// (See Types.h for an explicit definition of iterator types.)
//
// Meshes with boundary are handled by creating an additional face for each
// boundary loop (the method Face::isBoundary() determines whether a given
// face is a boundary loop).  Isolated vertices (i.e., vertiecs not contained
// in any edge or face) reference a dummy halfedge and can be checked via
// the method Vertex::isIsolated().
//

#ifndef DDG_MESH_H
#define DDG_MESH_H

#include <vector>
#include <string>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/upsample.h>
#include <igl/barycentric_coordinates.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/loop.h>

#include "HalfEdge.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "Complex.h"
#include <eigen3/Eigen/Dense>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Eigen/SparseCholesky>

#define WRITE_TEXTURE_COOR True

namespace DDG
{
   class Mesh
   {
      public:
         Mesh( void );
         // constructs an empty mesh

         Mesh( const Mesh& mesh );
         // constructs a copy of mesh

         const Mesh& operator=( const Mesh& mesh );
         // copies mesh

         int read( const std::string& filename, const std::string& filename2 );
         // reads a mesh from a Wavefront OBJ file; return value is nonzero
         // only if there was an error

         int write( const std::string& filename );
         // writes a mesh to a Wavefront OBJ file; return value is nonzero
         // only if there was an error

         bool reload( void );
         // reloads a mesh from disk using the most recent input filename
         
         void normalize( void );
         // centers around the origin and rescales to have unit radius

         int eulerCharacteristic( void ) const;
         // returns V-E+F

         std::vector<HalfEdge> halfedges;
         std::vector<Vertex>   vertices;
         std::vector<Edge>     edges;
         std::vector<Face>     faces;
         // storage for mesh elements

         unsigned int fieldDegree;
         double fieldOffsetAngle = 0;
         void computeSmoothestSection( void );
         void computeCurvatureAlignedSection( void );
         void computeTrivialSection( void );
         void extractSingularities( void );
         void alignTrivialSection( void );
         void computeParameterization( int coordinate );
         void glueParameterization( void );
         void parameterize( void );
         void buildEnergy( SparseMatrix<Real>& A, int coordinate );
         void buildDirichletEnergy( SparseMatrix<Real>& A );
         void energyGradient(SparseMatrix<Real>& A);
         void checkenergyGradient(int coordinate);
         void parameterizatioGradient(SparseMatrix<Real>& A, Eigen::MatrixXd& groundState, double l, Eigen::SparseMatrix<double>& H);
         double computeSmallestEigenValueMe(SparseMatrix<Real>& A, Eigen::MatrixXd& groundState, int index);
         void generateConformalMesh();
         double SampleTexturePointDerivative(int num, int triangle, Eigen::Vector2d& bary);
         void assignTextureCoordinatesDerivative(int coordinate);

         double energy( const SparseMatrix<Real>& A,
                        const DenseMatrix<Real>& x,
                        double eps );

         double lambda;
         // target line frequency

         int nCoordinateFunctions;
         // specifies whether to compute just the stripe
         // pattern (n=1) or both the stripe pattern and
         // an orthogonal coordinate (n=2), so that we
         // can construct a 2D parameterization

         bool target;
         Eigen::VectorXd v_tar;
         Eigen::VectorXd textureCoor;
         Eigen::VectorXd gs;
         Eigen::VectorXd ns;
         Eigen::MatrixXd dtexturedthetas;
         std::vector<Eigen::MatrixXd> dcolorsdx;
         Eigen::VectorXd colors;

      protected:
         std::string inputFilename;
         std::string inputFilename2;

         void indexElements( void );
         void initializeSmoothStructure( void );
         void buildMassMatrix( void );
         void buildFieldEnergy( void );
         void assignTextureCoordinates( int coordinate );
         void buildDualLaplacian( void );

         SparseMatrix<Complex> massMatrix;
         SparseMatrix<Real> realMassMatrix;
         SparseMatrix<Complex> energyMatrix;
         SparseMatrix<Real> dualLaplacian;

         int nComputedCoordinateFunctions = 1;
         // keeps track of how many coordinate functions were
         // computed the last time we ran Mesh::parameterize()
   };
}

#endif

