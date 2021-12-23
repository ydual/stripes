#include <map>
#include <fstream>
#include <queue>
#include "Mesh.h"
#include "MeshIO.h"
#include "Utility.h"


using namespace std;
using namespace DDGConstants;
using namespace Spectra;

namespace DDG
{
   Mesh :: Mesh( void )
   : fieldDegree( 2 ),  // line field, cross field, etc.---one can change this to compute n-direction fields for n other than 2, but currently parameterization only works for n=2
     lambda( 130. ), // initial global line frequency
     nCoordinateFunctions( 1 )
   {}

   Mesh :: Mesh( const Mesh& mesh )
   {
      *this = mesh;
   }

   class  HalfEdgeIterCompare { public: bool operator()( const  HalfEdgeIter& i, const  HalfEdgeIter& j ) const { return &*i < &*j; } };
   class HalfEdgeCIterCompare { public: bool operator()( const HalfEdgeCIter& i, const HalfEdgeCIter& j ) const { return &*i < &*j; } };
   class    VertexIterCompare { public: bool operator()( const    VertexIter& i, const    VertexIter& j ) const { return &*i < &*j; } };
   class   VertexCIterCompare { public: bool operator()( const   VertexCIter& i, const   VertexCIter& j ) const { return &*i < &*j; } };
   class      FaceIterCompare { public: bool operator()( const      FaceIter& i, const      FaceIter& j ) const { return &*i < &*j; } };
   class     FaceCIterCompare { public: bool operator()( const     FaceCIter& i, const     FaceCIter& j ) const { return &*i < &*j; } };
   class      EdgeIterCompare { public: bool operator()( const      EdgeIter& i, const      EdgeIter& j ) const { return &*i < &*j; } };
   class     EdgeCIterCompare { public: bool operator()( const     EdgeCIter& i, const     EdgeCIter& j ) const { return &*i < &*j; } };

   const Mesh& Mesh :: operator=( const Mesh& mesh )
   {
      map< HalfEdgeCIter, HalfEdgeIter, HalfEdgeCIterCompare > halfedgeOldToNew;
      map<   VertexCIter,   VertexIter,   VertexCIterCompare >   vertexOldToNew;
      map<     EdgeCIter,     EdgeIter,     EdgeCIterCompare >     edgeOldToNew;
      map<     FaceCIter,     FaceIter,     FaceCIterCompare >     faceOldToNew;

      // copy geometry from the original mesh and create a
      // map from pointers in the original mesh to
      // those in the new mesh
      halfedges.clear(); for( HalfEdgeCIter he = mesh.halfedges.begin(); he != mesh.halfedges.end(); he++ ) halfedgeOldToNew[ he ] = halfedges.insert( halfedges.end(), *he );
       vertices.clear(); for(   VertexCIter  v =  mesh.vertices.begin();  v !=  mesh.vertices.end();  v++ )   vertexOldToNew[ v  ] =  vertices.insert(  vertices.end(), *v  );
          edges.clear(); for(     EdgeCIter  e =     mesh.edges.begin();  e !=     mesh.edges.end();  e++ )     edgeOldToNew[ e  ] =     edges.insert(     edges.end(), *e  );
          faces.clear(); for(     FaceCIter  f =     mesh.faces.begin();  f !=     mesh.faces.end();  f++ )     faceOldToNew[ f  ] =     faces.insert(     faces.end(), *f  );

      // "search and replace" old pointers with new ones
      for( HalfEdgeIter he = halfedges.begin(); he != halfedges.end(); he++ )
      {
         he->next   = halfedgeOldToNew[ he->next   ];
         he->flip   = halfedgeOldToNew[ he->flip   ];
         he->vertex =   vertexOldToNew[ he->vertex ];
         he->edge   =     edgeOldToNew[ he->edge   ];
         he->face   =     faceOldToNew[ he->face   ];
      }

      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ ) v->he = halfedgeOldToNew[ *(v->he) ];
      for(   EdgeIter e =    edges.begin(); e !=    edges.end(); e++ ) e->he = halfedgeOldToNew[ e->he ];
      for(   FaceIter f =    faces.begin(); f !=    faces.end(); f++ ) f->he = halfedgeOldToNew[ f->he ];

      return *this;
   }

   int Mesh::read( const string& filename )
   {
      inputFilename = "../data/" + filename;
      ifstream in( filename.c_str() );

      if( !in.is_open() )
      {
         cerr << "Error reading from mesh file " << filename << endl;
         return 1;
      }

      int rval;
      if( !( rval = MeshIO::read( in, *this )))
      {
         indexElements();
         initializeSmoothStructure();
         buildMassMatrix();
      }
      return rval;
   }

   int Mesh::write( const string& filename )
   // reads a mesh from a Wavefront OBJ file; return value is nonzero
   // only if there was an error
   {
      ofstream out( filename.c_str() );

      // if we computed two orthogonal coordinate functions
      // (instead of just the 1D parameterization used to
      // draw stripes), identify edges in parameter space to
      // get a coherent (i.e., continuous almost everywhere)
      // parameterization
      if( nComputedCoordinateFunctions == 2 )
      {
         glueParameterization();
      }

      if( !out.is_open() )
      {
         cerr << "Error writing to mesh file " << filename << endl;
         return 1;
      }

      MeshIO::write( out, *this );

      return 0;
   }

   void Mesh :: indexElements( void )
   {
      int nV = 0;
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->index = nV;
         nV++;
      }

      int nE = 0;
      for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
      {
         e->index = nE;
         nE++;
      }

      int nF = 0;
      for( FaceIter f = faces.begin(); f != faces.end(); f++ )
      {
         f->index = nF;
         nF++;
      }
   }

   bool Mesh::reload( void )
   {
      return read( inputFilename );
   }

   void Mesh::normalize( void )
   {
      // compute center of mass
      Vector c( 0., 0., 0. );
      for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         c += v->position;
      }
      c /= (double) vertices.size();

      // translate to origin
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->position -= c;
      }

      // rescale such that the mesh sits inside the unit ball
      double rMax = 0.;
      for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         rMax = max( rMax, v->position.norm() );
      }
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->position /= rMax;
      }
   }

   int Mesh :: eulerCharacteristic( void ) const
   {
      int nV = (int) vertices.size();
      int nE = (int) edges.size();
      int nF = (int) faces.size();

      return nV - nE + nF;
   }

   void Mesh :: initializeSmoothStructure( void )
   {
      // compute angular coordinates of each outgoing halfedge
      for( VertexIter v  = vertices.begin();
                      v != vertices.end();
                      v ++ )
      {
         // compute the cumulative angle at each outgoing
         // halfedge, relative to the initial halfedge
         double cumulativeAngle = 0.;
         HalfEdgeIter he = *(v->he);
         do
         {
            he->angularCoordinate = cumulativeAngle;
            cumulativeAngle += he->next->angle();
            he = he->next->next->flip;
         }
         while( he != *(v->he) );

         // normalize angular coordinates so that they sum to two pi
         if( !v->onBoundary() )
         {
            do
            {
               he->angularCoordinate *= 2.*M_PI/cumulativeAngle;
               he = he->flip->next;
            }
            while( he != v->he );
         }
      }
   }

   void Mesh :: buildMassMatrix( void )
   {
      size_t nV = vertices.size();
      massMatrix.resize( nV, nV );
      realMassMatrix.resize( 2*nV, 2*nV );

      for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         int i = v->index;
         double A = v->dualArea();
         massMatrix( i, i ) = A;
         realMassMatrix( i*2+0, i*2+0 ) = A;
         realMassMatrix( i*2+1, i*2+1 ) = A;
         // realMassMatrix( i*2+0, i*2+0 ) = 1.0;
         // realMassMatrix( i*2+1, i*2+1 ) = 1.0;
      }
   }

   void Mesh :: buildFieldEnergy( void )
   // build per-face
   {
      double k = fieldDegree;

      SparseMatrix<Complex>& A( energyMatrix );
      size_t nV = vertices.size();
      A.resize( nV, nV );
      for( FaceCIter f = faces.begin(); f != faces.end(); f++ )
      {
         if( f->isBoundary() ) continue;

         HalfEdgeCIter he = f->he;
         do
         {
            int i = he->vertex->index;
            int j = he->flip->vertex->index;
            double w = he->cotan() / 2.;
            double thetaI = he->angularCoordinate;
            double thetaJ = he->flip->angularCoordinate;
            double phi = k*( thetaI - thetaJ + M_PI );
            Complex r( cos(phi), sin(phi) );

            A( i, i ) += w;
            A( i, j ) -= w*r;

            A( j, j ) += w;
            A( j, i ) -= w*r.inv();

            he = he->next;
         }
         while( he != f->he );
      }

      // Some domains will admit a trivial section (e.g., a flat disk),
      // hence we shift to make the matrix strictly positive-definite for
      // the Cholesky solver.  Note, however, that a constant shift will
      // not change the eigenvectors of the matrix, hence we get exactly
      // the same solution.  I.e., it is only the eigenvectors (and not the
      // eigenvalues) that are used in the end.
      A.shift( 1e-4 );
   }

   void Mesh :: buildDualLaplacian( void )
   {
      size_t nF = faces.size();
      SparseMatrix<Real>& L( dualLaplacian );
      L = SparseMatrix<Real>( nF, nF );

      for( FaceCIter f = faces.begin(); f != faces.end(); f++ )
      {
         int i = f->index;

         HalfEdgeCIter h = f->he;
         do
         {
            int j = h->flip->face->index;
            double wij = 1.; // 2. / ( cotAlpha + cotBeta );

            L( i, i ) += wij;
            L( i, j ) -= wij;

            h = h->next;
         }
         while( h != f->he );
      }

      L.shift( 1e-10 );
   }

   void Mesh :: computeTrivialSection( void )
   {
      // store previous section
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->oldDirectionField = v->directionField;
      }

      // solve for scalar potential
      buildDualLaplacian();
      int nV = (int) vertices.size();
      int nE = (int) edges.size();
      int nF = (int) faces.size();
      int chi = nV - nE + nF;
      DenseMatrix<Real> Omega( nF );
      double indexSum = 0.;
      for( FaceCIter f = faces.begin(); f != faces.end(); f++ )
      {
         Omega( f->index ) = -f->curvature() + 2.*M_PI*f->singularIndex;
         indexSum += f->singularIndex;
      }
      if( indexSum != chi )
      {
         // Fact of life: you can't comb the hair on a billiard ball...
         cerr << "WARNING: singularity indices do not satisfy Poincaré-Hopf!" << endl;
      }

      // extract connection 1-form
      DenseMatrix<Real> u( nF );
      solvePositiveDefinite( dualLaplacian, u, Omega );

      for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
      {
         int i = e->he->face->index;
         int j = e->he->flip->face->index;
         e->omega = u(j) - u(i);
      }

      // construct parallel section
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->visited = false;
      }
      vertices.begin()->directionField = Complex( cos(fieldOffsetAngle), sin(fieldOffsetAngle) );
      vertices.begin()->visited = true;
      queue<VertexIter> Q;
      Q.push( vertices.begin() );
      while( !Q.empty() )
      {
         VertexIter vi = Q.front(); Q.pop();
         HalfEdgeIter h = *(vi->he);
         do
         {
            VertexIter vj = h->flip->vertex;
            if( !vj->visited )
            {

               double thetaI = h->angularCoordinate;
               double thetaJ = h->flip->angularCoordinate + M_PI;
               double dTheta = thetaJ - thetaI;
               Complex rij( cos(dTheta), sin(dTheta) );

               double omegaIJ = -h->omega();
               Complex sij( cos(omegaIJ), sin(omegaIJ) );

               Complex Xi = vi->directionField;
               Complex Xj = rij*sij*Xi;

               vj->directionField = Xj;
               vj->visited = true;
               Q.push( vj );
            }

            h = h->flip->next;
         }
         while( h != vi->he );
      }
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->directionField *= v->directionField;
      }
   }

   void Mesh :: alignTrivialSection( void )
   // Suppose we take the singularities of the globally smoothest section and construct the
   // corresponding trivial connection.  A section parallel with respect to this connection
   // will resemble the smoothest section only up to a constant rotation in each tangent
   // space.  This method looks for the "best fit" rotation.
   {
      // compute the mean angle difference between the old and new section
      Complex mean( 0., 0. );
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         Complex z0 = v->oldDirectionField;
         Complex z1 = v->directionField;
         Complex w = z1.inv()*z0;
         mean += w;
      }
      mean.normalize();

      // align
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->directionField *= mean;
      }

      fieldOffsetAngle = vertices.begin()->directionField.arg()/2.;
   }

   enum class CSVState {
      UnquotedField,
      QuotedField,
      QuotedQuote
   };

   std::vector<std::string> readCSVRow(const std::string &row) {
      CSVState state = CSVState::UnquotedField;
      std::vector<std::string> fields {""};
      size_t i = 0; // index of the current field
      for (char c : row) {
         switch (state) {
               case CSVState::UnquotedField:
                  switch (c) {
                     case ',': // end of field
                                 fields.push_back(""); i++;
                                 break;
                     case '"': state = CSVState::QuotedField;
                                 break;
                     default:  fields[i].push_back(c);
                                 break; }
                  break;
               case CSVState::QuotedField:
                  switch (c) {
                     case '"': state = CSVState::QuotedQuote;
                                 break;
                     default:  fields[i].push_back(c);
                                 break; }
                  break;
               case CSVState::QuotedQuote:
                  switch (c) {
                     case ',': // , after closing quote
                                 fields.push_back(""); i++;
                                 state = CSVState::UnquotedField;
                                 break;
                     case '"': // "" -> "
                                 fields[i].push_back('"');
                                 state = CSVState::QuotedField;
                                 break;
                     default:  // end of quote
                                 state = CSVState::UnquotedField;
                                 break; }
                  break;
         }
      }
      return fields;
   }

   /// Read CSV file, Excel dialect. Accept "quoted fields ""with quotes"""
   std::vector<std::vector<std::string>> readCSV(std::istream &in) {
      std::vector<std::vector<std::string>> table;
      std::string row;
      while (!in.eof()) {
         std::getline(in, row);
         if (in.bad() || in.fail()) {
               break;
         }
         auto fields = readCSVRow(row);
         table.push_back(fields);
      }
      return table;
   }

   Eigen::Vector2d projectOntoPlane(const Eigen::Vector3d &vec,
                                 const Eigen::Vector3d &normal,
                                 const Eigen::Vector3d &axis)
   {
      Eigen::Vector3d eY = normal.cross(axis).normalized();
      Eigen::Vector3d eX = eY.cross(normal).normalized();
      return {eX.dot(vec), eY.dot(vec)};
   }
   
   void Mesh :: computeSmoothestSection( void )
   {
      cout << "Computing globally smoothest direction field..." << endl;
      // srand( 1234325 );

      buildFieldEnergy();

      size_t nV = vertices.size();
      DenseMatrix<Complex> groundState( nV );
      smallestEigPositiveDefinite( energyMatrix, massMatrix, groundState );

      
      std::string csv_file = inputFilename.substr(0,inputFilename.find(".obj"));
      csv_file+=".csv";
      std::cout<<csv_file<<std::endl;

      std::ifstream is(csv_file, std::ifstream::binary);
      std::vector<std::vector<std::string>> table = readCSV(is);

      // for(int i=0; i<table.size(); i++){
      //    for(int j=0; j<table[i].size(); j++){
      //       std::cout<<table[i][j]<<" ";
      //    }
      //    std::cout<<std::endl;
      // }

      for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
         double alpha = (*v->he)->angularCoordinate;
         Complex r(cos(2. * alpha), sin(2. * alpha));
         auto n = v->normal();
         auto e = (*v->he)->vector();
         Eigen::Vector3d row_i(std::stod(table[v->index][0]),std::stod(table[v->index][1]),std::stod(table[v->index][2]));
         auto u = projectOntoPlane(row_i.transpose(), {n.x, n.y, n.z}, {e.x, e.y, e.z});
         double a = std::atan2(u.y(), u.x()) + M_PI / 2.0;
         v->directionField = r * Complex(cos(2.0 * a), sin(2.0 * a));
      }


      // for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      // {
      //    Complex cur(std::stod(table[v->index][0]), std::stod(table[v->index][1]));
      //    v->directionField = cur.unit();
      // }
   }

   void Mesh :: computeCurvatureAlignedSection( void )
   {
      cerr << "Computing curvature-aligned direction field..." << endl;

      buildFieldEnergy();

      size_t nV = vertices.size();
      DenseMatrix<Complex> principalField( nV );
      DenseMatrix<Complex> smoothedField( nV );

      for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         principalField( v->index ) = v->principalDirection().unit();
      }

      SparseMatrix<Complex> A;
      const double t = 1.;
      A = energyMatrix + Complex(t)*massMatrix;

      solvePositiveDefinite( A, smoothedField, principalField );

      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->directionField = smoothedField( v->index ).unit();
      }
   }

   void Mesh :: parameterize( void )
   {
      computeParameterization( 0 );

      // at the user's request, we can also compute a second coordinate
      // aligned with the orthogonal direction field---these two coordinates
      // together describe a 2D parameterization rather than a 1D parameterization
      if( nCoordinateFunctions == 2 )
      {
         // rotate 1-vector field by 90 degrees, which in our complex representation
         // is equivalent to rotating a 2-vector field by 180 degrees
         for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
         {
            v->directionField = -v->directionField;
         }

         // compute the second coordinate
         computeParameterization( 1 );

         // rotate back
         for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
         {
            v->directionField = -v->directionField;
         }
      }

      // keep track of how many coordinates are actually valid
      nComputedCoordinateFunctions = nCoordinateFunctions;
   }

   void Mesh :: buildDirichletEnergy( SparseMatrix<Real>& A )
   {
      size_t nV = vertices.size();
      A.resize( 2*nV, 2*nV );

      for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
      {
         // get the endpoints
         VertexCIter vi = e->he->vertex;
         VertexCIter vj = e->he->flip->vertex;

         // get the angle of the edge w.r.t. the endpoints' bases
         double thetaI = e->he->angularCoordinate;
         double thetaJ = e->he->flip->angularCoordinate + M_PI;

         // compute the parallel transport coefficient from i to j
         double dTheta = thetaJ - thetaI;
         Complex rij( cos(dTheta), sin(dTheta) );

         // compute the cotan weight
         double cotAlpha = e->he->cotan();
         double cotBeta  = e->he->flip->cotan();
         if(       e->he->face->fieldIndex(2.) != 0 ) cotAlpha = 0.;
         if( e->he->flip->face->fieldIndex(2.) != 0 ) cotBeta  = 0.;
         double w = (cotAlpha+cotBeta)/2.;

         // pick an arbitrary root at each endpoint
         Complex Xi = vi->canonicalVector();
         Complex Xj = vj->canonicalVector();

         // check if the roots point the same direction
         double s = dot( rij*Xi, Xj ) > 0. ? 1. : -1.;
         if( fieldDegree == 1 ) s = 1.;
         if( rand()%200 == 0 ) s = -1; else s = 1.;
         e->crossesSheets = ( s < 0. );

         int i = 2 * vi->index;
         int j = 2 * vj->index;

         // add the diagonal terms
         A(i+0,i+0) += w;
         A(i+1,i+1) += w;

         A(j+0,j+0) += w;
         A(j+1,j+1) += w;

         // if both vectors pointed the same direction, use a
         // 2x2 block that represents complex multiplication
         if( s > 0. )
         {
            A(i+0,j+0) = -w; A(i+0,j+1) = 0.;
            A(i+1,j+0) = 0.; A(i+1,j+1) = -w;

            A(j+0,i+0) = -w; A(j+0,i+1) = 0.;
            A(j+1,i+0) = 0.; A(j+1,i+1) = -w;
         }
         // otherwise, use a block that represents both
         // complex conjugation and multiplication
         else
         {
            A(i+0,j+0) = +w; A(i+0,j+1) = 0.;
            A(i+1,j+0) = 0.; A(i+1,j+1) = +w;

            A(j+0,i+0) = +w; A(j+0,i+1) = 0.;
            A(j+1,i+0) = 0.; A(j+1,i+1) = +w;
         }
      }

      A.shift( 1e-4 );
   }

   void Mesh :: buildEnergy( SparseMatrix<Real>& A, int coordinate )
   {
      size_t nV = vertices.size();
      A.resize( 2*nV, 2*nV );

      for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
      {
         // get the endpoints
         VertexCIter vi = e->he->vertex;
         VertexCIter vj = e->he->flip->vertex;

         // get the angle of the edge w.r.t. the endpoints' bases
         double thetaI = e->he->angularCoordinate;
         double thetaJ = e->he->flip->angularCoordinate + M_PI;

         // compute the parallel transport coefficient from i to j
         double dTheta = thetaJ - thetaI;
         Complex rij( cos(dTheta), sin(dTheta) );

         // compute the cotan weight
         double cotAlpha = e->he->cotan();
         double cotBeta  = e->he->flip->cotan();
         if(       e->he->face->fieldIndex(2.) != 0 ) cotAlpha = 0.;
         if( e->he->flip->face->fieldIndex(2.) != 0 ) cotBeta  = 0.;
         double w = (cotAlpha+cotBeta)/2.;

         // pick an arbitrary root at each endpoint
         Complex Xi = vi->canonicalVector();
         Complex Xj = vj->canonicalVector();

         // check if the roots point the same direction
         double s = dot( rij*Xi, Xj ) > 0. ? 1. : -1.;
         if( fieldDegree == 1 ) s = 1.;
         e->crossesSheets = ( s < 0. );

         // compute the 1-form value along edge ij
         double lij = e->length();
         double phiI = (Xi).arg();
         double phiJ = (s*Xj).arg();
         double omegaIJ = lambda * (lij/2.) * ( cos(phiI-thetaI) + cos(phiJ-thetaJ) );

         e->omega = omegaIJ;

         // compute the components of the new transport coefficient
         double a = w * cos(omegaIJ);
         double b = w * sin(omegaIJ);

         int i = 2 * vi->index;
         int j = 2 * vj->index;

         // add the diagonal terms
         A(i+0,i+0) += w;
         A(i+1,i+1) += w;

         A(j+0,j+0) += w;
         A(j+1,j+1) += w;

         // if both vectors pointed the same direction, use a
         // 2x2 block that represents complex multiplication
         if( s > 0. )
         {
            A(i+0,j+0) = -a; A(i+0,j+1) = -b;
            A(i+1,j+0) =  b; A(i+1,j+1) = -a;

            A(j+0,i+0) = -a; A(j+0,i+1) =  b;
            A(j+1,i+0) = -b; A(j+1,i+1) = -a;
         }
         // otherwise, use a block that represents both
         // complex conjugation and multiplication
         else
         {
            A(i+0,j+0) = -a; A(i+0,j+1) =  b;
            A(i+1,j+0) =  b; A(i+1,j+1) =  a;

            A(j+0,i+0) = -a; A(j+0,i+1) =  b;
            A(j+1,i+0) =  b; A(j+1,i+1) =  a;
         }
      }

      A.shift( 1e-4 );
   }

   void Mesh::energyGradient(SparseMatrix<Real>& A)
   {
      size_t nV = vertices.size();
      for(int i=0; i<nV; ++i){
         vertices[i].dAdX.resize(2);
         vertices[i].dAdX[0].resize(2*nV,2*nV);
         vertices[i].dAdX[0].setZero();
         vertices[i].dAdX[1].resize(2*nV,2*nV);
         vertices[i].dAdX[0].setZero();
      }
         

      for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
      {
         // get the endpoints
         VertexIter vi = e->he->vertex;
         VertexIter vj = e->he->flip->vertex;

         // get the angle of the edge w.r.t. the endpoints' bases
         double thetaI = e->he->angularCoordinate;
         double thetaJ = e->he->flip->angularCoordinate + M_PI;

         // compute the parallel transport coefficient from i to j
         double dTheta = thetaJ - thetaI;
         Complex rij( cos(dTheta), sin(dTheta) );

         // compute the cotan weight
         double cotAlpha = e->he->cotan();
         double cotBeta  = e->he->flip->cotan();
         if(       e->he->face->fieldIndex(2.) != 0 ) cotAlpha = 0.;
         if( e->he->flip->face->fieldIndex(2.) != 0 ) cotBeta  = 0.;
         double w = (cotAlpha+cotBeta)/2.;

         // pick an arbitrary root at each endpoint
         Complex Xi = vi->canonicalVector();
         Complex Xj = vj->canonicalVector();

         // check if the roots point the same direction
         double s = dot( rij*Xi, Xj ) > 0. ? 1. : -1.;
         if( fieldDegree == 1 ) s = 1.;
         e->crossesSheets = ( s < 0. );

         // compute the 1-form value along edge ij
         double lij = e->length();
         double phiI = (Xi).arg();
         double phiJ = (s*Xj).arg();

         double dphiIdia = - 0.5/(1+(vi->directionField.im/vi->directionField.re)*(vi->directionField.im/vi->directionField.re))*vi->directionField.im/(vi->directionField.re*vi->directionField.re);
         double dphiIdib =  0.5/(1+(vi->directionField.im/vi->directionField.re)*(vi->directionField.im/vi->directionField.re))*1/(vi->directionField.re);
         double dphiJdja = - 0.5/(1+(vj->directionField.im/vj->directionField.re)*(vj->directionField.im/vj->directionField.re))*vj->directionField.im/(vj->directionField.re*vj->directionField.re);
         double dphiJdjb =  0.5/(1+(vj->directionField.im/vj->directionField.re)*(vj->directionField.im/vj->directionField.re))*1/(vj->directionField.re);




         double omegaIJ = lambda * (lij/2.) * ( cos(phiI-thetaI) + cos(phiJ-thetaJ) );
         double domegaIJdia = lambda * (lij/2.) *  -sin(phiI-thetaI) *  dphiIdia;
         double domegaIJdib = lambda * (lij/2.) *  -sin(phiI-thetaI) *  dphiIdib;
         double domegaIJdja = lambda * (lij/2.) *  -sin(phiJ-thetaJ) *  dphiJdja;
         double domegaIJdjb = lambda * (lij/2.) *  -sin(phiJ-thetaJ) *  dphiJdjb;

         e->omega = omegaIJ;

         // compute the components of the new transport coefficient
         double a = w * cos(omegaIJ);
         double b = w * sin(omegaIJ);

         int i = 2 * vi->index;
         int j = 2 * vj->index;

         // if(vi->index == 0 || vj->index == 0){
         //    std::cout<<w<<" "<<omegaIJ<<" "<<dphiIdia<<" "<<dphiIdib<<std::endl;
         // }


         if( s > 0. )
         {
            // A(i+0,j+0) = -a; A(i+0,j+1) = -b;
            // A(i+1,j+0) =  b; A(i+1,j+1) = -a;

            // A(j+0,i+0) = -a; A(j+0,i+1) =  b;
            // A(j+1,i+0) = -b; A(j+1,i+1) = -a;

            vi->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+0,j+1) = -w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](i+0,j+1) = -w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+1,j+1) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](i+1,j+1) = w * sin(omegaIJ) * domegaIJdib;

            vj->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+0,j+1) = -w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](i+0,j+1) = -w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+1,j+1) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](i+1,j+1) = w * sin(omegaIJ) * domegaIJdjb;

            vi->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+1,i+0) = -w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](j+1,i+0) = -w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+1,i+1) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](j+1,i+1) = w * sin(omegaIJ) * domegaIJdib;

            vj->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+1,i+0) = -w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](j+1,i+0) = -w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+1,i+1) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](j+1,i+1) = w * sin(omegaIJ) * domegaIJdjb;
         }
         else
         {
            // A(i+0,j+0) = -a; A(i+0,j+1) =  b;
            // A(i+1,j+0) =  b; A(i+1,j+1) =  a;

            // A(j+0,i+0) = -a; A(j+0,i+1) =  b;
            // A(j+1,i+0) =  b; A(j+1,i+1) =  a;


            vi->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+0,j+1) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](i+0,j+1) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdib;

            vj->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+0,j+1) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](i+0,j+1) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdjb;

            vi->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+1,i+0) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](j+1,i+0) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdib;

            vj->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+1,i+0) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](j+1,i+0) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdjb;
         }
      }
   }

   void Mesh::checkenergyGradient(int coordinate)
   {
      SparseMatrix<Real> A;
      buildEnergy( A, coordinate );
      std::vector<SparseMatrix<Real>> dAdX;
      energyGradient(A);
      
      
      double eps = 1e-7;
      for(int i=0; i<vertices.size(); ++i)
      {
         Complex Xi = vertices[i].canonicalVector();
         double phiI = (Xi).arg();

         vertices[i].directionField.re += eps;
         SparseMatrix<Real> A_prime;
         buildEnergy( A_prime, coordinate );
         for(int j=0; j<2*vertices.size(); ++j){
            for(int k=0; k<2*vertices.size(); ++k){
               if(true){
                  // std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                  // std::cout <<vertices[i].dAdX[0](j,k)<<" "<<(A_prime(j,k)-A(j,k))/eps<<" "<<vertices[i].dAdX[0](j,k)-(A_prime(j,k)-A(j,k))/eps<<std::endl;
                  if(fabs(vertices[i].dAdX[0](j,k)-(A_prime(j,k)-A(j,k))/eps) >1e-4) std::cout<<i<<" "<<j<<" "<<k<<" "<<vertices[i].dAdX[0](j,k)<<" "<<(A_prime(j,k)-A(j,k))/eps<<std::endl;
               }
            }
         }

         vertices[i].directionField.re -= eps;


         vertices[i].directionField.im += eps;
         buildEnergy( A_prime, coordinate );
         for(int j=0; j<2*vertices.size(); ++j){
            for(int k=0; k<2*vertices.size(); ++k){
               if((A_prime(j,k)-A(j,k))/eps != 0){
                  // std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                  // std::cout <<vertices[i].dAdX[1](j,k)<<" "<<(A_prime(j,k)-A(j,k))/eps<<" "<<vertices[i].dAdX[1](j,k)-(A_prime(j,k)-A(j,k))/eps<<std::endl;
                  if(vertices[i].dAdX[1](j,k)-(A_prime(j,k)-A(j,k))/eps >1e-4) std::cout<<i<<" "<<j<<" "<<k<<" "<<vertices[i].dAdX[1](j,k)<<" "<<(A_prime(j,k)-A(j,k))/eps<<std::endl;
               }
            }
         }
         vertices[i].directionField.im -= eps;
      }
   }

   void Mesh::parameterizatioGradient(SparseMatrix<Real>& A, Eigen::MatrixXd& groundState, double l, Eigen::SparseMatrix<double>& H)
   {
      size_t nV = vertices.size();

      Eigen::MatrixXd RA(2*nV,2*nV);
      Eigen::MatrixXd C(2*nV+1,2*nV+1);
      Eigen::MatrixXd B(2*nV,2*nV);
      Eigen::VectorXd x(2*nV);
      double alpha  = 0;
      for(int i=0; i<2*nV; ++i){
         for(int j=0; j<2*nV; ++j){
            RA(i,j) = A(i,j);
            C(i,j) = A(i,j) - l*realMassMatrix(i,j);
            B(i,j) = realMassMatrix(i,j);
         }
      }


      Eigen::MatrixXd Bx = B*groundState;
      for(int i=0; i<2*nV; ++i){
         C(2*nV,i) = -Bx(i);
         C(i,2*nV) = -Bx(i);
      }


      Eigen::MatrixXd NC(4*nV+1, 4*nV+1);
      NC.setZero();
      for(int i=0; i<2*nV+1; ++i)
         NC(i,i) = 1;
      for(int i=2*nV+1; i<4*nV+1; ++i)
         NC(i,i) = alpha;
      for(int i=0; i<2*nV; ++i){
         for(int j=0; j<2*nV; ++j){
            NC(i,j+2*nV+1) = (-l*realMassMatrix(i,j)+A(i,j));
            NC(j+2*nV+1, i) = (-l*realMassMatrix(i,j)+A(i,j));
         }
      }
      for(int i=0; i<2*nV; ++i){
         NC(2*nV,i+2*nV+1) = Bx(i);
         NC(i+2*nV+1,2*nV) = Bx(i);
      }

      // Eigen::SparseMatrix<double> CS(2*nV+1, 2*nV+1);
      // CS.reserve(Eigen::VectorXi::Constant(2*nV+1,6));
      // for(int i=0; i<=2*nV; ++i){
      //    for(int j=0; j<=2*nV; ++j){
      //       // if(fabs(C(i,j)) > 1e-6)
      //          CS.insert(i,j) = C(i,j);
      //    }
      // }
      // CS.makeCompressed();


      Eigen::SparseMatrix<double> NCS(4*nV+1, 4*nV+1);
      NCS.reserve(Eigen::VectorXi::Constant(4*nV+1,6));
      for(int i=0; i<4*nV+1; ++i){
         for(int j=0; j<4*nV+1; ++j){
            // if(fabs(C(i,j)) > 1e-6)
               NCS.insert(i,j) = NC(i,j);
         }
      }
      NCS.makeCompressed();
      H = NCS;
      

      // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverLDLT;
      // solverLDLT.compute(CS);
      // solverLDLT.analyzePattern(CS);

      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverLDLT2;
      solverLDLT2.compute(NCS);
      solverLDLT2.analyzePattern(NCS);

      // int np = 0;
      // int nn = 0;
      // int nz = 0;

      // double EPS_ZERO = 1e-6;
      // double EPS_ZERO_PSD = 1e-6;
      // for(int i=0; i<solverLDLT.vectorD().size(); i++){
      //    double di = solverLDLT.vectorD()(i);
      //    if(di<EPS_ZERO_PSD) nn++;
      //    else if(fabs(di)<EPS_ZERO_PSD) nz++;
      //    else np++;
      // }

      // std::cout<<"nn: "<<nn<<std::endl;
      // std::cout<<"nz: "<<nz<<std::endl;
      // std::cout<<"np: "<<np<<std::endl;



      for(int i=0; i<nV; ++i){
         vertices[i].dxdX.resize(2);
         vertices[i].dxdX[0].resize(4*nV+1,1);
         vertices[i].dxdX[1].resize(4*nV+1,1);

         
         Eigen::VectorXd bl(4*nV+1),blp(4*nV+1);
         Eigen::VectorXd bl1(2*nV);
         for(int j=0; j<2*nV; ++j) bl1[j] = 0.0;
         bl << bl1, 0,-(vertices[i].dAdX[0]*groundState);
         vertices[i].dxdX[0] =  solverLDLT2.solve(bl);
         blp << bl1, 0,-(vertices[i].dAdX[1]*groundState);
         vertices[i].dxdX[1] =  solverLDLT2.solve(blp);
      }
   }

   double Mesh::computeSmallestEigenValueMe(SparseMatrix<Real>& A, Eigen::MatrixXd& groundState, int index)
   {
      int nV = vertices.size();

      Eigen::SparseMatrix<double> RA(2*nV, 2*nV);
      Eigen::SparseMatrix<double> B(2*nV, 2*nV);
      RA.reserve(Eigen::VectorXi::Constant(2*nV,6));
      B.reserve(Eigen::VectorXi::Constant(2*nV,6));

      for(int i=0; i<2*nV; ++i){
         for(int j=0; j<2*nV; ++j){
            if(A(i,j) != 0)
               RA.insert(i,j) = A(i,j);
            if(realMassMatrix(i,j)!= 0)
               B.insert(i,j) = realMassMatrix(i,j);
         }
      }
      RA.makeCompressed();
      B.makeCompressed();

      //Eigen::MatrixXd D = B.inverse()*RA;


      // Construct matrix operation object using the wrapper class
      SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse> op1(RA,B);
      SparseGenMatProd<double> op2(B);
   
      // Construct eigen solver object, requesting the largest
      // (in magnitude, or norm) three eigenvalues
      double sigma = 1e-5;
      SymGEigsShiftSolver<SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>,SparseGenMatProd<double>,GEigsMode::ShiftInvert> eigs(op1, op2, index, min(index+20,2*nV), sigma);
   
      // Initialize and compute
      eigs.init();
      int nconv = eigs.compute(SortRule::LargestMagn);

      auto evecs = eigs.eigenvectors().col(0);
      double evalues = eigs.eigenvalues()(0);
      std::cout<<eigs.eigenvalues()<<std::endl;

      //std::cout<<"fuck"<<std::endl;
      //std::cout<<RA*eigs.eigenvectors().col(0)-eigs.eigenvalues()(0)*B*evecs<<std::endl;

      groundState = evecs;
      return evalues;

   }

   void Mesh :: computeParameterization( int coordinate )
   {
      cerr << "Computing stripe pattern..." << endl;

      // srand( time( NULL ) );
      // srand( 1234567 );

      SparseMatrix<Real> A;
      buildEnergy( A, coordinate );
      energyGradient(A);
      checkenergyGradient(coordinate);
      

      size_t nV = vertices.size();
      int nEigs = 9;
      double l;
      std::vector<double> D;
      std::vector<DenseMatrix<Real>> V;

   
      DenseMatrix<Real> groundState(2*nV);
      Eigen::MatrixXd gs;
      //smallestEigPositiveDefinite( A, realMassMatrix, groundState );
      nSmallestEigsPositiveDefinite(A,realMassMatrix,V,D,nEigs);
      l = D[0];
      groundState = V[0];

      Eigen::MatrixXd gs1(2*nV,1);
      for(int i=0; i<2*nV; ++i)
         gs1(i) = groundState(i);

      Eigen::SparseMatrix<double> H(4*nV+1, 4*nV+1);
      parameterizatioGradient(A,gs1,D[0], H);

      

      double eps = 1e-5;
      for(int iter=0; iter<nV; ++iter){
         vertices[iter].directionField.re += eps;

         SparseMatrix<Real> A_prime;
         buildEnergy( A_prime, coordinate );

         Eigen::MatrixXd B(2*nV,2*nV);
         Eigen::MatrixXd dA(2*nV,2*nV);
         Eigen::MatrixXd A_p(2*nV,2*nV);
            

         for(int i=0; i<2*nV; ++i)
            for(int j=0; j<2*nV; ++j){
               A_p(i,j) = A_prime(i,j);
               dA(i,j) = A_prime(i,j)-A(i,j);
               B(i,j) = realMassMatrix(i,j);
            }

         Eigen::MatrixXd Bx = B*gs1;

         Eigen::MatrixXd NC(4*nV+1, 4*nV+1);
         NC.setZero();
         for(int i=0; i<2*nV+1; ++i)
            NC(i,i) = 1;
         for(int i=2*nV+1; i<4*nV+1; ++i)
            NC(i,i) = 0;
         for(int i=0; i<2*nV; ++i){
            for(int j=0; j<2*nV; ++j){
               NC(i,j+2*nV+1) = (-l*realMassMatrix(i,j)+A(i,j));
               NC(j+2*nV+1, i) = (-l*realMassMatrix(i,j)+A(i,j));
            }
         }
         for(int i=0; i<2*nV; ++i){
            NC(2*nV,i+2*nV+1) = Bx(i);
            NC(i+2*nV+1,2*nV) = Bx(i);
         }

         std::cout<<std::endl;

         Eigen::SparseMatrix<double> NCS(4*nV+1, 4*nV+1);
         NCS.reserve(Eigen::VectorXi::Constant(4*nV+1,6));
         for(int i=0; i<4*nV+1; ++i){
            for(int j=0; j<4*nV+1; ++j){
               // if(fabs(C(i,j)) > 1e-6)
                  NCS.insert(i,j) = NC(i,j);
            }
         }
         NCS.makeCompressed();

               
         Eigen::VectorXd bl1(2*nV), bl3(2*nV), bl4(2*nV);
         for(int j=0; j<2*nV; ++j) bl1[j] = 0.0;
         bl3 = dA*gs1;
         bl4 = vertices[iter].dAdX[0]*gs1*eps;


         Eigen::VectorXd bl(4*nV+1);
         bl << -bl1, 0,-bl3;
         Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverLDLT2;
         solverLDLT2.compute(NCS);
         solverLDLT2.analyzePattern(NCS);
         Eigen::VectorXd dxdXs_n = solverLDLT2.solve(bl);

         vertices[iter].directionField.re -= 2*eps;

         buildEnergy( A_prime, coordinate );

         for(int i=0; i<2*nV; ++i)
            for(int j=0; j<2*nV; ++j){
               A_p(i,j) = A_prime(i,j);
               dA(i,j) = A_prime(i,j)-A(i,j);
               B(i,j) = realMassMatrix(i,j);
            }

         NC.setZero();
         for(int i=0; i<2*nV+1; ++i)
            NC(i,i) = 1;
         for(int i=2*nV+1; i<4*nV+1; ++i)
            NC(i,i) = 0;
         for(int i=0; i<2*nV; ++i){
            for(int j=0; j<2*nV; ++j){
               NC(i,j+2*nV+1) = (-l*realMassMatrix(i,j)+A(i,j));
               NC(j+2*nV+1, i) = (-l*realMassMatrix(i,j)+A(i,j));
            }
         }
         for(int i=0; i<2*nV; ++i){
            NC(2*nV,i+2*nV+1) = Bx(i);
            NC(i+2*nV+1,2*nV) = Bx(i);
         }

         Eigen::SparseMatrix<double> NCSp(4*nV+1, 4*nV+1);
         NCSp.reserve(Eigen::VectorXi::Constant(4*nV+1,4*nV+1));
         for(int i=0; i<4*nV+1; ++i){
            for(int j=0; j<4*nV+1; ++j){
               // if(fabs(C(i,j)) > 1e-6)
                  NCSp.insert(i,j) = NC(i,j);
            }
         }
         NCSp.makeCompressed();

         bl3 = (A_p-l*B)*gs1;
         Eigen::VectorXd blp(4*nV+1);
         blp << -bl1, 0,-bl3;
         Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverLDLT2p;
         solverLDLT2p.compute(NCSp);
         solverLDLT2p.analyzePattern(NCSp);

         Eigen::VectorXd dxdXs_np = solverLDLT2p.solve(blp);


         for(int j=0; j<2*nV; ++j){
            double dxdX_a = vertices[iter].dxdX[0](j);
            double dxdX_n = (dxdXs_n[j]-dxdXs_np[j])/(2*eps);
            if(dxdX_a != 0)
               std::cout<<iter<<" "<<j<<" analytical: "<<dxdX_a<<" numerical: "<<dxdX_n<<" diff: "<<(dxdX_n-dxdX_a)/dxdX_a<<std::endl;
         }

         vertices[iter].directionField.re += eps;
      }

      
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         int i = v->index;
         v->parameterization = Complex( groundState( i*2+0 ),
                                        groundState( i*2+1 ) ).unit();
      }

      
      assignTextureCoordinates( coordinate );

            

   }

   double Mesh :: energy( const SparseMatrix<Real>& A, const DenseMatrix<Real>& x, double eps )
   {
      // evaluate quadratic part of energy
      double mu = inner( A*x, x );

      // evaluate quartic part of energy
      int nV = (int) vertices.size();
      for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         if( v->index != nV-1 )
         {
            Complex psi = v->parameterization;
            mu += eps * sqr( psi.norm2() - 1. ) / 4.;
         }
      }

      return mu;
   }

   void Mesh :: assignTextureCoordinates( int p )
   // This method computes the final texture coordinates that are actually
   // used by OpenGL to draw the stripes, starting with the solution to
   // the eigenvalue problem.
   {
      for( FaceIter f = faces.begin(); f != faces.end(); f++ )
      {
         if( f->isBoundary() ) continue;

         // grab the halfedges
         HalfEdgeIter hij = f->he;
         HalfEdgeIter hjk = hij->next;
         HalfEdgeIter hki = hjk->next;

         // grab the parameter values at vertices
         Complex psiI = hij->vertex->parameterization;
         Complex psiJ = hjk->vertex->parameterization;
         Complex psiK = hki->vertex->parameterization;

         double cIJ = ( hij->edge->he != hij ? -1. : 1. );
         double cJK = ( hjk->edge->he != hjk ? -1. : 1. );
         double cKI = ( hki->edge->he != hki ? -1. : 1. );

         // grab the connection coeffients
         double omegaIJ = cIJ * hij->edge->omega;
         double omegaJK = cJK * hjk->edge->omega;
         double omegaKI = cKI * hki->edge->omega;

         if( hij->edge->crossesSheets )
         {
            psiJ = psiJ.bar();
            omegaIJ =  cIJ * omegaIJ;
            omegaJK = -cJK * omegaJK;
         }

         if( hki->edge->crossesSheets )
         {
            psiK = psiK.bar();
            omegaKI = -cKI * omegaKI;
            omegaJK =  cJK * omegaJK;
         }

         // construct complex transport coefficients
         Complex rij( cos(omegaIJ), sin(omegaIJ) );
         Complex rjk( cos(omegaJK), sin(omegaJK) );
         Complex rki( cos(omegaKI), sin(omegaKI) );

         // compute the angles at the triangle corners closest to the target omegas
         double alphaI = psiI.arg();
         double alphaJ = alphaI + omegaIJ - (rij*psiI/psiJ).arg(); //fmodPI((varphiI + omegaIJ) - varphiJ); // could do this in terms of angles instead of complex numbers...
         double alphaK = alphaJ + omegaJK - (rjk*psiJ/psiK).arg(); //fmodPI((varphiJ + omegaJK) - varphiK); // mostly a matter of taste---possibly a matter of performance?
         double alphaL = alphaK + omegaKI - (rki*psiK/psiI).arg(); //fmodPI((varphiK + omegaKI) - varphiI);

         // adjust triangles containing zeros
         double n = lround((alphaL-alphaI)/(2.*M_PI));
         alphaJ -= 2.*M_PI*n/3.;
         alphaK -= 4.*M_PI*n/3.;

         // store the coordinates
         hij->texcoord[p] = alphaI;
         hjk->texcoord[p] = alphaJ;
         hki->texcoord[p] = alphaK;
         f->paramIndex[p] = n;
      }
   }

   void Mesh :: glueParameterization( void )
   // Running this method is completely optional - if you compute two orthogonal
   // stripe coordinates, it simply glues together the triangles in the parameter
   // domain so that they describe a continuous parameterization (otherwise, each
   // triangle can be different from its neighbor by appropriate reflections/
   // rotations/translations).  It is not usually needed to visualize the stripe
   // patterns.  But in some rendering packages, for certain effects like
   // displacement mapping, it can reduce rendering artifacts to have a continuous
   // parameterization, rather than one where each triangle is its own little island.
   {
      queue<FaceIter> Q;

      // start at an arbitrary nonsingular face
      FaceIter f0 = faces.begin();
      while( f0->paramIndex[0] != 0. ||
             f0->paramIndex[1] != 0. ||
             f0->fieldIndex(2.) != 0. )
      {
         f0++;
      }
      f0->visited = true;
      Q.push( f0 );

      while( !Q.empty() )
      {
         FaceIter f = Q.front(); Q.pop();

         HalfEdgeIter he = f->he;
         do
         {
            FaceIter fj = he->flip->face;
            // traverse neighbors only if they're nonsingular
            if( !f->isBoundary() &&
                !fj->visited &&
                fj->paramIndex[0]  == 0. &&
                fj->paramIndex[1]  == 0. &&
                fj->fieldIndex(2.) == 0. )
            {
               // grab handles to the three neighboring vertex coordinates
               Complex& bi = he->flip->next->texcoord;
               Complex& bj = he->flip->texcoord;
               Complex& bk = he->flip->next->next->texcoord;

               // compute the two parameter-space vectors along the shared edge
               Complex ai = he->texcoord;
               Complex aj = he->next->texcoord;
               Complex u = aj-ai;
               Complex v = bj-bi;

               // compute the rotation between these two edges
               double theta = (u*v.inv()).arg();
               Complex z( cos(theta), sin(theta) );

               // if the angle is anything other than 0 or 180...
               if( fabs(z.im) > 1e-9 )
               {
                  // ...apply a reflection to the neighboring triangle
                  bi=bi.bar();
                  bj=bj.bar();
                  bk=bk.bar();
               }

               // now recompute the rotation
               v = bj-bi;
               theta = (u*v.inv()).arg();
               z = Complex( cos(theta), sin(theta) );

               // as long as we now have a valid rotation...
               if( fabs(z.im) < 1e-9 )
               {
                  // ...rotate and translate the neighbor so that
                  // the shared edge matches up with the parent
                  Complex b0 = bi;
                  bi = z*(bi-b0) + ai;
                  bj = z*(bj-b0) + ai;
                  bk = z*(bk-b0) + ai;
               }

               // enqueue the neighbor
               fj->visited = true;
               Q.push( fj );
            }

            he = he->next;
         }
         while( he != f->he );
      }
   }

   void Mesh :: extractSingularities( void )
   {
      for( FaceIter f = faces.begin(); f != faces.end(); f++ )
      {
         f->singularIndex = f->fieldIndex( fieldDegree ) / (double) fieldDegree;
      }
   }
}

