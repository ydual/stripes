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
     nCoordinateFunctions( 1 ),
     target(true)
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

   int Mesh::read( const string& filename , const string& filename2)
   {
      inputFilename = "../data/" + filename;
      inputFilename2 = "../data/" + filename2;
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
      return read( inputFilename, inputFilename2 );
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

      
      // std::string csv_file = inputFilename.substr(0,inputFilename.find(".obj"));
      // csv_file+=".csv";
      std::string csv_file = inputFilename2;
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

      std::cout<<"Succeed Reading!"<<std::endl;


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

   double ifBranchDiff(double x, double& derivative)
   {
      double rou = 10000.;
      double t1 = 1.+exp(rou*x);
      double dt1dx = rou*exp(rou*x);
      double t2 = 2./t1;
      double dt2dt1 = -2./t1*t1;
      double result = -t2+1.;  
      derivative = -dt2dt1*dt1dx;
      return result;
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
         // double s = dot( rij*Xi, Xj ) > 0. ? 1. : -1.;
         // if( fieldDegree == 1 ) s = 1.;
         // e->crossesSheets = ( s < 0. );

         double dotc = dot( rij*Xi, Xj );
         double amp = -100.;
         double s = 1.-2.*(exp(amp*dotc)/(1.+exp(amp*dotc)));
         
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

         //std::cout<<i<<" "<<j<<" "<<w<<std::endl; 
         // add the diagonal terms
         
         A(i+0,i+0) += w;
         A(i+1,i+1) += w;

         A(j+0,j+0) += w;
         A(j+1,j+1) += w;

         // if both vectors pointed the same direction, use a
         // 2x2 block that represents complex multiplication
         if( true )
         {
            A(i+0,j+0) = -a; A(i+0,j+1) = -b*s;
            A(i+1,j+0) =  b; A(i+1,j+1) = -a*s;

            A(j+0,i+0) = -a; A(j+0,i+1) =  b;
            A(j+1,i+0) = -b*s; A(j+1,i+1) = -a*s;
         }
         // otherwise, use a block that represents both
         // complex conjugation and multiplication
         // else
         // {
         //    A(i+0,j+0) = -a; A(i+0,j+1) =  b;
         //    A(i+1,j+0) =  b; A(i+1,j+1) =  a;

         //    A(j+0,i+0) = -a; A(j+0,i+1) =  b;
         //    A(j+1,i+0) =  b; A(j+1,i+1) =  a;
         // }
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
         vertices[i].dAdX[1].setZero();
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

         // !!!! Have to do that
         if(       e->he->face->fieldIndex(2.) != 0 ) cotAlpha = 0.;
         if( e->he->flip->face->fieldIndex(2.) != 0 ) cotBeta  = 0.;
         double w = (cotAlpha+cotBeta)/2.;

         // pick an arbitrary root at each endpoint
         Complex Xi = vi->canonicalVector();
         Complex Xj = vj->canonicalVector();

         double ri = Xi.norm();
         double thetai = Xi.arg();
         double rj = Xj.norm();
         double thetaj = Xj.arg();


         // check if the roots point the same direction
         // double s = dot( rij*Xi, Xj ) > 0. ? 1. : -1.;
         // if( fieldDegree == 1 ) s = 1.;
         // e->crossesSheets = ( s < 0. );

         double dotc = dot( rij*Xi, Xj );
         double amp = -100.;
         double s = 1.-2.*(exp(amp*dotc)/(1.+exp(amp*dotc)));

         double dsdt1 = -2.*exp(-amp*dotc)*amp/((1.+exp(-amp*dotc))*(1.+exp(-amp*dotc)));
         
         double dt1dri = rj*cos(thetai+dTheta)*cos(thetaj) + rj*sin(thetai+dTheta)*sin(thetaj);
         double dt1drj = ri*cos(thetai+dTheta)*cos(thetaj) + ri*sin(thetai+dTheta)*sin(thetaj);
         double dt1dthetai = -ri*rj*sin(thetai+dTheta)*cos(thetaj)+ ri*rj*cos(thetai+dTheta)*sin(thetaj);
         double dt1dthetaj = -ri*rj*cos(thetai+dTheta)*sin(thetaj)+ ri*rj*sin(thetai+dTheta)*cos(thetaj);

         double ia = vi->directionField.re;
         double ib = vi->directionField.im;
         double ja = vj->directionField.re;
         double jb = vj->directionField.im;

         double dridia = ia/sqrt(ia*ia+ib*ib);
         double dridib = ib/sqrt(ia*ia+ib*ib);
         double dthetaidia = 0.5*1./(1+(ib/ia)*(ib/ia))*-ib/(ia*ia);
         double dthetaidib = 0.5*1./(1+(ib/ia)*(ib/ia))*1./ia;

         double drjdja = ja/sqrt(ja*ja+jb*jb);
         double drjdjb = jb/sqrt(ja*ja+jb*jb);
         double dthetajdja = 0.5*1./(1+(jb/ja)*(jb/ja))*-jb/(ja*ja);
         double dthetajdjb = 0.5*1./(1+(jb/ja)*(jb/ja))*1./ja;

         double dsdia = dsdt1*(dt1dri*dridia+dt1dthetai*dthetaidia);
         double dsdib = dsdt1*(dt1dri*dridib+dt1dthetai*dthetaidib);
         double dsdja = dsdt1*(dt1drj*drjdja+dt1dthetaj*dthetajdja);
         double dsdjb = dsdt1*(dt1drj*drjdjb+dt1dthetaj*dthetajdjb);

         // double s = ifBranchDiff(t1,dsdt1);
         // if( fieldDegree == 1 ) s = 1.;
         // e->crossesSheets = ( s < 0. );

         // compute the 1-form value along edge ij
         double lij = e->length();
         double phiI = (Xi).arg();
         double phiJ = (s*Xj).arg();

         // if(s < 0){
         //    std::cout<<phiJ<<" "<<Xj.arg()<<" "<<sqrt(Xj.arg()*Xj.arg()+0.000000001)<<std::endl;
         // }

         
         // std::cout<<s<<" "<<(s*Xj).arg()<<" "<<Xj.arg()<<std::endl;
         // std::cout<<(s*Xj).re<<" "<<(s*Xj).im<<" "<<Xj.re<<" "<<Xj.im<<std::endl;

         double dphiIdia = - 0.5/(1+(vi->directionField.im/vi->directionField.re)*(vi->directionField.im/vi->directionField.re))*vi->directionField.im/(vi->directionField.re*vi->directionField.re);
         double dphiIdib =  0.5/(1+(vi->directionField.im/vi->directionField.re)*(vi->directionField.im/vi->directionField.re))*1/(vi->directionField.re);
         double dphiJdja = - 0.5/(1+(vj->directionField.im/vj->directionField.re)*(vj->directionField.im/vj->directionField.re))*vj->directionField.im/(vj->directionField.re*vj->directionField.re);
         double dphiJdjb =  0.5/(1+(vj->directionField.im/vj->directionField.re)*(vj->directionField.im/vj->directionField.re))*1/(vj->directionField.re);

         // double ire = vi->directionField.re;
         // double iim = vi->directionField.im;
         // double jre = vj->directionField.re;
         // double jim = vj->directionField.im;

         // double dphiIdia = 0.5*(-iim/(ire*ire+iim*iim));
         // double dphiIdib =  0.5*(ire/(ire*ire+iim*iim));
         // double dphiJdja = 0.5*(-jim/(jre*jre+jim*jim));
         // double dphiJdjb =  0.5*(jre/(jre*jre+jim*jim));




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

         // // add the diagonal terms
         // A(i+0,i+0) += w;
         // A(i+1,i+1) += w;

         // A(j+0,j+0) += w;
         // A(j+1,j+1) += w;


         if( true )
         {
            // A(i+0,j+0) = -a; A(i+0,j+1) = -b*s;
            // A(i+1,j+0) =  b; A(i+1,j+1) = -a*s;

            // A(j+0,i+0) = -a; A(j+0,i+1) =  b;
            // A(j+1,i+0) = -b*s; A(j+1,i+1) = -a*s;

            vi->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+0,j+1) = - s * w * cos(omegaIJ) * domegaIJdia - b * dsdia; vi->dAdX[1](i+0,j+1) = - s * w * cos(omegaIJ) * domegaIJdib - b * dsdib;
            vi->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](i+1,j+1) = s * w * sin(omegaIJ) * domegaIJdia - a * dsdia; vi->dAdX[1](i+1,j+1) = s * w * sin(omegaIJ) * domegaIJdib - a * dsdib;

            vj->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+0,j+1) = - s * w * cos(omegaIJ) * domegaIJdja - b * dsdja; vj->dAdX[1](i+0,j+1) = - s * w * cos(omegaIJ) * domegaIJdjb - b * dsdjb;
            vj->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](i+1,j+1) = s * w * sin(omegaIJ) * domegaIJdja - a * dsdja; vj->dAdX[1](i+1,j+1) = s * w * sin(omegaIJ) * domegaIJdjb - a * dsdjb;

            vi->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdib;
            vi->dAdX[0](j+1,i+0) = - s * w * cos(omegaIJ) * domegaIJdia - b * dsdia; vi->dAdX[1](j+1,i+0) = - s * w * cos(omegaIJ) * domegaIJdib - b * dsdib;
            vi->dAdX[0](j+1,i+1) = s * w * sin(omegaIJ) * domegaIJdia - a * dsdia; vi->dAdX[1](j+1,i+1) = s * w * sin(omegaIJ) * domegaIJdib - a * dsdib;

            vj->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdjb;
            vj->dAdX[0](j+1,i+0) = - s * w * cos(omegaIJ) * domegaIJdja - b * dsdja; vj->dAdX[1](j+1,i+0) = - s * w * cos(omegaIJ) * domegaIJdjb - b * dsdjb;
            vj->dAdX[0](j+1,i+1) = s * w * sin(omegaIJ) * domegaIJdja - a * dsdja; vj->dAdX[1](j+1,i+1) = s * w * sin(omegaIJ) * domegaIJdjb - a * dsdjb;
         }
         // else
         // {
         //    // A(i+0,j+0) = -a; A(i+0,j+1) =  b;
         //    // A(i+1,j+0) =  b; A(i+1,j+1) =  a;

         //    // A(j+0,i+0) = -a; A(j+0,i+1) =  b;
         //    // A(j+1,i+0) =  b; A(j+1,i+1) =  a;


         //    vi->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdib;
         //    vi->dAdX[0](i+0,j+1) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](i+0,j+1) = w * cos(omegaIJ) * domegaIJdib;
         //    vi->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdib;
         //    vi->dAdX[0](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdib;

         //    vj->dAdX[0](i+0,j+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](i+0,j+0) = w * sin(omegaIJ) * domegaIJdjb;
         //    vj->dAdX[0](i+0,j+1) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](i+0,j+1) = w * cos(omegaIJ) * domegaIJdjb;
         //    vj->dAdX[0](i+1,j+0) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](i+1,j+0) = w * cos(omegaIJ) * domegaIJdjb;
         //    vj->dAdX[0](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](i+1,j+1) = -w * sin(omegaIJ) * domegaIJdjb;

         //    vi->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdib;
         //    vi->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdib;
         //    vi->dAdX[0](j+1,i+0) = w * cos(omegaIJ) * domegaIJdia; vi->dAdX[1](j+1,i+0) = w * cos(omegaIJ) * domegaIJdib;
         //    vi->dAdX[0](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdia; vi->dAdX[1](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdib;

         //    vj->dAdX[0](j+0,i+0) = w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+0) = w * sin(omegaIJ) * domegaIJdjb;
         //    vj->dAdX[0](j+0,i+1) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](j+0,i+1) = w * cos(omegaIJ) * domegaIJdjb;
         //    vj->dAdX[0](j+1,i+0) = w * cos(omegaIJ) * domegaIJdja; vj->dAdX[1](j+1,i+0) = w * cos(omegaIJ) * domegaIJdjb;
         //    vj->dAdX[0](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdja; vj->dAdX[1](j+1,i+1) = -w * sin(omegaIJ) * domegaIJdjb;
         // }
      }
   }

   void Mesh::checkenergyGradient(int coordinate)
   {
      SparseMatrix<Real> A;
      buildEnergy( A, coordinate );
      std::vector<SparseMatrix<Real>> dAdX;
      //energyGradient(A);
      
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
                  //std::cout<<(A_prime(j,k)-A(j,k))/eps<<" ";
                  if(fabs(vertices[i].dAdX[0](j,k)-(A_prime(j,k)-A(j,k))/eps) >1e-4) std::cout<<i<<" "<<j<<" "<<k<<" "<<vertices[i].dAdX[0](j,k)<<" "<<(A_prime(j,k)-A(j,k))/eps<<std::endl;
               }
            }
            //std::cout<<endl;
         }

         //std::cout<<endl;
         vertices[i].directionField.re -= eps;


         vertices[i].directionField.im += eps;
         buildEnergy( A_prime, coordinate );
         for(int j=0; j<2*vertices.size(); ++j){
            for(int k=0; k<2*vertices.size(); ++k){
               if(true){
                  // std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                  // std::cout <<vertices[i].dAdX[1](j,k)<<" "<<(A_prime(j,k)-A(j,k))/eps<<" "<<vertices[i].dAdX[1](j,k)-(A_prime(j,k)-A(j,k))/eps<<std::endl;
                  if(fabs(vertices[i].dAdX[1](j,k)-(A_prime(j,k)-A(j,k))/eps) >1e-4) std::cout<<i<<" "<<j<<" "<<k<<" "<<vertices[i].dAdX[1](j,k)<<" "<<(A_prime(j,k)-A(j,k))/eps<<std::endl;
                  //std::cout<<(A_prime(j,k)-A(j,k))/eps<<" ";
               }
            }
            //std::cout<<endl;
         }
         //std::cout<<endl;
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
      NCS.reserve(Eigen::VectorXi::Constant(4*nV+1,4*nV+1));
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

      Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solverLDLT2;
      solverLDLT2.compute(NCS);
      //solverLDLT2.analyzePattern(NCS);

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

      auto NC_inv = NC.inverse();

      // std::cout<<NC_inv*NC<<std::endl;
      // std::cout<<std::endl;


      //std::cout<<target<<std::endl;

      for(int i=0; i<nV; ++i){
         vertices[i].dxdX.resize(2);
         vertices[i].dxdX[0].resize(2*nV,1);
         vertices[i].dxdX[1].resize(2*nV,1);
         
         Eigen::VectorXd bl(4*nV+1),blp(4*nV+1);
         Eigen::VectorXd bl1(2*nV);
         for(int j=0; j<2*nV; ++j) bl1[j] = 0.0;
         bl << bl1, 0,-(vertices[i].dAdX[0]*groundState);

         // std::cout<<(vertices[i].dAdX[0]*gst).transpose()-(vertices[i].dAdX[0]*groundState).transpose()<<std::endl;
         // std::cout<<std::endl;
         Eigen::MatrixXd temp1 = solverLDLT2.solve(bl);

         //std::cout<<i<<std::endl;
         for(int j=0; j<2*nV; ++j)
            vertices[i].dxdX[0](j)=  temp1(j);
         
         //std::cout<<vertices[i].dxdX[0].transpose()<<std::endl;

         blp << bl1, 0,-(vertices[i].dAdX[1]*groundState);
         Eigen::MatrixXd temp2 = solverLDLT2.solve(blp);
         for(int j=0; j<2*nV; ++j)
            vertices[i].dxdX[1](j)=  temp2(j);
      }
   }

   double Mesh::computeSmallestEigenValueMe(SparseMatrix<Real>& A, Eigen::MatrixXd& groundState, int index)
   {
      int nV = vertices.size();
      groundState.resize(2*nV, index);

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
      //std::cout<<evecs<<std::endl;

      //std::cout<<eigs.eigenvectors().col(1).transpose()<<std::endl;
      // std::cout<<(RA*eigs.eigenvectors().col(0)-eigs.eigenvalues()(0)*B*eigs.eigenvectors().col(0)).transpose()<<std::endl;
      // std::cout<<std::endl;
      // if(index > 1)
      //    std::cout<<(RA*eigs.eigenvectors().col(1)-eigs.eigenvalues()(0)*B*eigs.eigenvectors().col(1)).transpose()<<std::endl;
      // std::cout<<std::endl;

      for(int i=0; i<index; ++i)
         for(int j=0; j<2*nV; ++j)
            groundState(j,i) = eigs.eigenvectors().col(i)(j);
      return evalues;

   }

   void Mesh :: computeParameterization( int coordinate )
   {
      //cerr << "Computing stripe pattern..." << endl;

      // srand( time( NULL ) );
      // srand( 1234567 );

      SparseMatrix<Real> A;
      buildEnergy( A, coordinate );
      energyGradient(A);
      checkenergyGradient(coordinate);

      std::cout<<"checkenergyGradient Complete"<<std::endl;
      

      size_t nV = vertices.size();
      int nEigs = 9;
      double l;
      std::vector<double> D;
      std::vector<DenseMatrix<Real>> V;

   
      DenseMatrix<Real> groundState(2*nV);
      Eigen::MatrixXd gs1(2*nV,1);

      if(target){
         //smallestEigPositiveDefinite( A, realMassMatrix, groundState );
         //nSmallestEigsPositiveDefinite(A,realMassMatrix,V,D,nEigs);
         // l = D[0];
         // groundState = V[0];

         l = computeSmallestEigenValueMe(A,gs1,1);
         for(int i=0; i<2*nV; ++i)
            groundState(i) = gs1(i);
      }
      else
      {
         //std::cout<<v_tar<<std::endl;
         Eigen::MatrixXd gs_p;
         l = computeSmallestEigenValueMe(A,gs_p,2);
         Eigen::VectorXd v1(2*nV),v2(2*nV);
         for(int i=0; i<2*nV; ++i){
            v1(i) = gs_p(i,0);
            v2(i) = gs_p(i,1);
         }

         double a = v1.transpose()*v1;
         double b = v1.transpose()*v2;
         double c = v2.transpose()*v2;
         double m = v_tar.transpose()*v1;
         double n = v_tar.transpose()*v2;

         double x1 = (c*m-b*n)/(a*c-b*b);
         double x2 = -(b*m-a*n)/(a*c-b*b);

         // std::cout<<x1<<" "<<x2<<std::endl;

         gs1 = x1*v1+x2*v2;

      }
      
      Eigen::SparseMatrix<double> H(4*nV+1, 4*nV+1);
      parameterizatioGradient(A,gs1,l,H);
      

      if(target) v_tar = gs1;

      //std::cout<<gs1<<std::endl;

      // double eps = 1e-5;
      // for(int iter=0; iter<nV; ++iter){
      //    vertices[iter].directionField.re += eps;

      //    SparseMatrix<Real> A_prime;
      //    buildEnergy( A_prime, coordinate );

      //    Eigen::MatrixXd gs_p;
      //    double lp = computeSmallestEigenValueMe(A_prime,gs_p,2);
      //    Eigen::VectorXd v0(2*nV),v1(2*nV),v2(2*nV);
      //    for(int i=0; i<2*nV; ++i){
      //       v0(i) = gs1(i);
      //       v1(i) = gs_p(i,0);
      //       v2(i) = gs_p(i,1);
      //    }

      //    double a = v1.transpose()*v1;
      //    double b = v1.transpose()*v2;
      //    double c = v2.transpose()*v2;
      //    double m = v0.transpose()*v1;
      //    double n = v0.transpose()*v2;

      //    // std::cout<<a<<" "<<b<<" "<<c<<std::endl;

      //    // double x1 = (c*m-b*n)/(a*c-b*b);
      //    // double x2 = -(b*m-a*n)/(a*c-b*b);

      //    double x1 = m/a;
      //    double x2 = n/c;

      //    Eigen::VectorXd vp;

      //    vp = x1*v1+x2*v2;

      //    vertices[iter].directionField.re -= 2*eps;
      //    SparseMatrix<Real> A_primep;
      //    buildEnergy( A_primep, coordinate );

      //    Eigen::MatrixXd gs_pp;
      //    double lpp = computeSmallestEigenValueMe(A_primep,gs_pp,2);
      //    Eigen::VectorXd v0p(2*nV),v1p(2*nV),v2p(2*nV);
      //    for(int i=0; i<2*nV; ++i){
      //       v0p(i) = gs1(i);
      //       v1p(i) = gs_pp(i,0);
      //       v2p(i) = gs_pp(i,1);
      //    }

      //    double ap = v1p.transpose()*v1p;
      //    double bp = v1p.transpose()*v2p;
      //    double cp = v2p.transpose()*v2p;
      //    double mp = v0p.transpose()*v1p;
      //    double np = v0p.transpose()*v2p;

      //    double x1p = (cp*mp-bp*np)/(ap*cp-bp*bp);
      //    double x2p = -(bp*mp-ap*np)/(ap*cp-bp*bp);

      //    Eigen::VectorXd vpp;

      //    vpp = x1p*v1p+x2p*v2p;


      //    for(int j=0; j<2*nV; ++j){
      //       double dxdX_a = vertices[iter].dxdX[0](j);
      //       double dxdX_n = (vp[j]-vpp[j])/(2*eps);
      //       if(dxdX_a != 0)
      //          std::cout<<iter<<" "<<j<<" analytical: "<<dxdX_a<<" numerical: "<<dxdX_n<<" diff: "<<(dxdX_n-dxdX_a)/dxdX_a<<std::endl;
      //    }

      //    vertices[iter].directionField.re += eps;
      // }

      std::vector<Eigen::MatrixXd> duxdx(nV);
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         int i = v->index;
         v->parameterization = Complex( gs1( i*2+0 ),
                                        gs1( i*2+1 ) );
         // Eigen::VectorXd f(2);
         // f<<gs1( i*2+0 ),gs1( i*2+1 );
         // Eigen::VectorXd n = f.normalized();
         // double s = f.norm();
         // Eigen::MatrixXd I(2,2);
         // I.setZero();
         // I(0,0) = 1; I(1,1) = 1;
         // Eigen::MatrixXd duxdx_v = (I-n*n.transpose())/s;

         // vertices[i].dnxdX.resize(2);
         // vertices[i].dnxdX[0].resize(2*nV,1);
         // vertices[i].dnxdX[1].resize(2*nV,1);

         // duxdx[i] = duxdx_v;
      }

      // for(int i=0; i<nV; ++i){
      //    for(int j=0; j<nV; ++j){
      //       for(int k=0; k<2; ++k){
      //          vertices[i].dnxdX[k](2*j+0) = vertices[i].dxdX[k](2*j+0)*duxdx[j](0,0)+vertices[i].dxdX[k](2*j+1)*duxdx[j](0,1);
      //          vertices[i].dnxdX[k](2*j+1) = vertices[i].dxdX[k](2*j+0)*duxdx[j](1,0)+vertices[i].dxdX[k](2*j+1)*duxdx[j](1,1);
      //       }
      //    }
      // }

      
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         if(v->dXdD.size() > 0){
            Eigen::MatrixXd temp = v->dxdX[0]*v->dXdD[0].transpose() + v->dxdX[1]*v->dXdD[1].transpose();
            v->dnxdD = temp;
            //std::cout<<temp.rows()<<" "<<temp.cols()<<std::endl;
         }
      }

      gs = gs1;

      
      assignTextureCoordinatesDerivative( coordinate );

      int nF = faces.size();
      Eigen::MatrixXd dtdnx(3*nF,2*nV);
      Eigen::MatrixXd dndnx(nF,2*nV);
      for(int i=0; i<nV; ++i){
         for(int j=0; j<3*nF; ++j){
            dtdnx(j,2*i) = vertices[i].dtexturecoordnx[0](j);
            dtdnx(j,2*i+1) = vertices[i].dtexturecoordnx[1](j);
         }
         for(int j=0; j<nF; ++j){
            dndnx(j,2*i) = vertices[i].dndnx[0](j);
            dndnx(j,2*i+1) = vertices[i].dndnx[1](j);
         }
      }
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         if(v->dXdD.size() > 0){
            v->dtexturecoordD = dtdnx*v->dnxdD;
            v->dndD = dndnx*v->dnxdD;
         }
         
      }

      generateConformalMesh();
      // std::vector<Eigen::MatrixXd> dcolorsdx_a = dcolorsdx;
      // Eigen::VectorXd colors_o = colors;
      

      // for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      // {
      //    double theta_o = v->parameterization.arg();
      //    double rou_o = v->parameterization.norm();
      //    v->parameterization.re += eps;
      //    assignTextureCoordinatesDerivative( coordinate );
      //    generateConformalMesh();
      //    Eigen::VectorXd colors_p = colors;
      //    Eigen::VectorXd dcolorsdx_n = (colors_p-colors_o)/eps;
      //    int j = v->index;
      //    for(int i=0; i<colors.size(); ++i)
      //    {
      //       //std::cout<<(v->parameterization.arg()-theta_o)/eps_theta<<" "<<(v->parameterization.norm()-rou_o)/eps_theta<<std::endl;
      //       //if(abs(dcolorsdx_n(i)) > 1e-10)
      //          std::cout<<i<<" "<<j<<" "<<dcolorsdx_a[0](i,j)<<" "<<dcolorsdx_n(i)<<std::endl;
      //    }
      //    v->parameterization.re -=eps;  
      // }

   }

   double fmodPIDiff(double theta, double& derivative)
   {
      //theta - (2.*M_PI) * std::floor( (theta+M_PI) / (2.*M_PI) );

      double x = (theta+M_PI) / (2.*M_PI);
      double dxdtheta = 1./(2.*M_PI);

      double phi = 0.00000001;

      double y1 = (2.*x-1)/4.;
      double y2 = (x)/2.;
      
      double m1 = (1.-phi)*sin(2.*M_PI*y1);
      double dm1dy1 = (1.-phi)*2.*M_PI*cos(2.*M_PI*y1);
      double T1 = 1.-(2.*acos(m1))/M_PI;
      double dT1dm1 = 2./(M_PI*sqrt(1.-m1*m1));
      
      double m2 = sin(2.*M_PI*y2)/phi;
      double dm2dy2  = 2.*M_PI*cos(2.*M_PI*y2);
      double T2 = 2.*(atan(m2))/M_PI;
      double dT2dm2 = 2./(M_PI*(1+m2*m2));

      double T3 = (1.+T1*T2)/2.;
      double dT3dx = 1./4.0*(T2*dT1dm1*dm1dy1+T1*dT2dm2*dm2dy2);

      double n =  x-T3;
      double dndx = 1.-dT3dx;
      
      derivative = 1.-(2.*M_PI)*dndx*dxdtheta;
      return theta-(2.*M_PI)*n;
   }

   void checkfmodPIDiff(double theta)
   {
      double doutdtheta = 0;
      double dummy = 0;
      double out1 = fmodPIDiff(theta, doutdtheta);
      double eps = 1e-5;

      double thetap = theta+eps;
      double thetapp = theta-eps;
      double outp = fmodPIDiff(thetap, dummy);
      double outpp = fmodPIDiff(thetapp, dummy);

      std::cout<<fmodPI(theta)<<std::endl;
      std::cout<<theta<<" "<<out1<<" "<<doutdtheta<<" "<<(outp-outpp)/(2*eps)<<std::endl;
   }

   double lroundDiff(double x, double& derivative)
   {

      double phi = 0.00000001;

      double y1 = (2.*(x-0.5)-1)/4.;
      double y2 = (x-0.5)/2.;
      
      double m1 = (1.-phi)*sin(2.*M_PI*y1);
      double dm1dy1 = (1.-phi)*2.*M_PI*cos(2.*M_PI*y1);
      double T1 = 1-(2.*acos((1.-phi)*sin(2.*M_PI*y1)))/M_PI;
      double dT1dm1 = 2./(M_PI*sqrt(1.-m1*m1));
      
      double m2 = sin(2.*M_PI*y2)/phi;
      double dm2dy2  = 2.*M_PI*cos(2.*M_PI*y2);
      double T2 = 2.*(atan(sin(2.*M_PI*y2)/phi))/M_PI;
      double dT2dm2 = 2./(M_PI*(1.+m2*m2));

      double T3 = (1.+T1*T2)/2.;
      double dT3dx = 1./4.*(T2*dT1dm1*dm1dy1+T1*dT2dm2*dm2dy2);

      double n =  x-T3+0.5;
      derivative = 1-dT3dx;
      
      return n;
   }

   void checkflroundDiff(double x)
   {
      double doutdx = 0;
      double dummy = 0;
      double out1 = lroundDiff(x, doutdx);
      double eps = 1e-5;

      double xp = x+eps;
      double xpp = x-eps;
      double outp = lroundDiff(xp, dummy);
      double outpp = lroundDiff(xpp, dummy);

      std::cout<<x<<" "<<out1<<" "<<doutdx<<" "<<(outp-outpp)/(2*eps)<<std::endl;
   }

   bool isValid (double a, double b, double c)
   {
      if(a<0 || a>1 || b<0 || b>1 || c<0 || c>1) return false;
      else return true;
   }

   double Mesh::SampleTexturePointDerivative(int num, int triangle, Eigen::Vector2d& bary)
   {
      Face fa = faces[triangle];
      double nu = fa.paramIndex[0];
      double k = fa.fieldIndex(2.);

      double s = 30;
      double w = 0.6;
      double f = 1.0;
      double tx;

      Eigen::Vector3d barycentric(bary(0), bary(1), 1.0-bary(0)-bary(1));

      // Get the three half edges.
      HalfEdgeCIter hij = fa.he;
      HalfEdgeCIter hjk = hij->next;
      HalfEdgeCIter hkl = hjk->next;

      // Get the three vertices.
      VertexCIter vi = hij->vertex;
      VertexCIter vj = hjk->vertex;
      VertexCIter vk = hkl->vertex;

      // Get the index of corners
      int indexI = vi->index;
      int indexJ = vj->index;
      int indexK = vk->index;

      Eigen::Vector3d dtxdtexture;
      std::vector<Eigen::VectorXd> dtxdx(2);
      dtxdx[0].setZero(3);
      dtxdx[1].setZero(3);

      //std::cout<<num<<" "<<indexI<<" "<<indexJ<<" "<<indexK<<std::endl;
      // std::cout<<std::endl;

      k = 0;

      if(k == 0)
      {
         Eigen::Vector3d txv(textureCoor(3*triangle), textureCoor(3*triangle+1), textureCoor(3*triangle+2));
         tx = barycentric.dot(txv);

         dtxdtexture = barycentric;
      }
      else
      {
         nu = 0.0;

         // Get the three parameter values---for clarity, let "l"
         // denote the other point in the same fiber as "i".  Purely
         // for clarity, we will explicitly define the value of psi
         // at l, which of course is always just the conjugate of the
         // value of psi at i.
         Complex psiI = vi->parameterization;
         Complex psiJ = vj->parameterization;
         Complex psiK = vk->parameterization;
         Complex psiL = psiI.bar();

         double cIJ = ( hij->edge->he != hij ? -1. : 1. );
         double cJK = ( hjk->edge->he != hjk ? -1. : 1. );
         double cKL = ( hkl->edge->he != hkl ? -1. : 1. );

         // Get the three omegas, which were used to define our energy.
         double omegaIJ = hij->omega();
         double omegaJK = hjk->omega();
         double omegaKL = hkl->omega();

         // Here's the trickiest part.  If the canonical orientation of
         // this last edge is from l to k (rather than from k to l)...
         omegaKL *= cKL;

         int revert_d = 1;
         int revert_f = 1;

         // Now we just get consecutive values along the curve from i to j to k to l.
         // (The following logic was already described in our routine for finding
         // zeros of the parameterization.)
         if( hij->crossesSheets() )
         {
            psiJ = psiJ.bar();
            revert_d = -1;
            omegaIJ =  cIJ * omegaIJ;
            omegaJK = -cJK * omegaJK;
         }

         // Note that the flag hkl->crossesSheets() is the opposite of what we want here:
         // based on the way it was originally computed, it flags whether the vectors at
         // Xk and Xi have a negative dot product.  But here, we instead want to know if
         // the vectors at Xk and Xl have a negative dot product.  (And since Xi=-Xl, this
         // flag will be reversed.)
         if( !hkl->crossesSheets() )
         {
            psiK = psiK.bar();
            revert_f = -1;
            omegaKL = -cKL * omegaKL;
            omegaJK =  cJK * omegaJK;
         }

         double thetaI = psiI.arg();
         double thetaJ = psiJ.arg();
         double thetaK = psiK.arg();
         double thetaL = psiL.arg();

         // From here, everything gets computed as usual.
         Complex rij( cos(omegaIJ), sin(omegaIJ) );
         Complex rjk( cos(omegaJK), sin(omegaJK) );
         Complex rkl( cos(omegaKL), sin(omegaKL) );

         double de1,de2,de3;

         // double sigmaIJ = omegaIJ - ((rij*psiI)/psiJ).arg();
         // double sigmaJK = omegaJK - ((rjk*psiJ)/psiK).arg();
         // double sigmaKL = omegaKL - ((rkl*psiK)/psiL).arg();
         //double xi = sigmaIJ + sigmaJK + sigmaKL;

         double sigmaIJ = omegaIJ - fmodPIDiff((omegaIJ+thetaI-thetaJ),de1);
         double sigmaJK = omegaJK - fmodPIDiff((omegaJK+thetaJ-thetaK),de2);
         double sigmaKL = omegaKL - fmodPIDiff((omegaKL+thetaK-thetaL),de3);

         double betaI = psiI.arg();
         double betaJ = betaI + sigmaIJ;
         double betaK = betaJ + sigmaJK;
         double betaL = betaK + sigmaKL;
         double betaM = betaI + (betaL-betaI)/2.;

         Eigen::MatrixXd dbetasdthetas(5,3);
         dbetasdthetas.setZero();
         dbetasdthetas(0,0) = 1.;
         dbetasdthetas(1,0) = dbetasdthetas(0,0) - de1;
         dbetasdthetas(1,1) = de1;
         dbetasdthetas(2,0) = dbetasdthetas(1,0);
         dbetasdthetas(2,1) = dbetasdthetas(1,1) - de2;
         dbetasdthetas(2,2) = de2;
         dbetasdthetas(3,0) = dbetasdthetas(2,0) - de3;  // Is it correct?
         dbetasdthetas(3,1) = dbetasdthetas(2,1);
         dbetasdthetas(3,2) = dbetasdthetas(2,2) - de3;
         dbetasdthetas(4,0) = 0.5*(dbetasdthetas(0,0) + dbetasdthetas(3,0));
         dbetasdthetas(4,1) = 0.5*(dbetasdthetas(0,1) + dbetasdthetas(3,1));
         dbetasdthetas(4,2) = 0.5*(dbetasdthetas(0,2) + dbetasdthetas(3,2));

         //std::cout<<dbetasdthetas<<std::endl;

         Eigen::Vector3d dtxdthetas;

         if(barycentric(2) <= barycentric(1) && barycentric(2) <= barycentric(0))
         {
            Eigen::Matrix2d trans;
            trans<<2,1,1,2;

            Eigen::Vector2d new_bary = trans*bary+Eigen::Vector2d(-1.0,-1.0);
            Eigen::Vector3d new_barycentric(new_bary(0), new_bary(1), 1.0-new_bary(0)-new_bary(1));

            Eigen::Vector3d txv(betaI, betaJ, betaM);
            tx = new_barycentric.dot(txv);

            for(int i=0; i<3; ++i)
               dtxdthetas(i) = new_barycentric(0)*dbetasdthetas(0,i)+new_barycentric(1)*dbetasdthetas(1,i)+new_barycentric(2)*dbetasdthetas(4,i);
         }
         else if(barycentric(0) <= barycentric(1) && barycentric(0) <= barycentric(2))
         {
            Eigen::Matrix2d trans;
            trans<<-1,1,-2,-1;

            Eigen::Vector2d new_bary = trans*bary+Eigen::Vector2d(0.0,1.0);
            Eigen::Vector3d new_barycentric(new_bary(0), new_bary(1), 1.0-new_bary(0)-new_bary(1));

            Eigen::Vector3d txv(betaJ, betaK, betaM);
            tx = new_barycentric.dot(txv);

            for(int i=0; i<3; ++i)
               dtxdthetas(i) = new_barycentric(0)*dbetasdthetas(1,i)+new_barycentric(1)*dbetasdthetas(2,i)+new_barycentric(2)*dbetasdthetas(4,i);
         }
         else if(barycentric(1) <= barycentric(0) && barycentric(1) <= barycentric(2))
         {
            Eigen::Matrix2d trans;
            trans<<-1,-2,1,-1;

            Eigen::Vector2d new_bary = trans*bary+Eigen::Vector2d(1.0,0.0);
            Eigen::Vector3d new_barycentric(new_bary(0), new_bary(1), 1.0-new_bary(0)-new_bary(1));

            Eigen::Vector3d txv(betaK, betaL, betaM);
            tx = new_barycentric.dot(txv);

            for(int i=0; i<3; ++i)
               dtxdthetas(i) = new_barycentric(0)*dbetasdthetas(2,i)+new_barycentric(1)*dbetasdthetas(3,i)+new_barycentric(2)*dbetasdthetas(4,i);
         }

         dtxdx[0](0) = -psiI.im/psiI.norm2()*dtxdthetas(0);
         dtxdx[1](0) = psiI.re/psiI.norm2()*dtxdthetas(0);
         dtxdx[0](1) = -psiJ.im/psiJ.norm2()*dtxdthetas(1);
         dtxdx[1](1) = psiJ.re/psiJ.norm2()*dtxdthetas(1)*revert_d;
         dtxdx[0](2) = -psiK.im/psiK.norm2()*dtxdthetas(2);
         dtxdx[1](2) = psiK.re/psiK.norm2()*dtxdthetas(2)*revert_f;
      }

      // Sample
      double theta = tx;

      double lArg = 0.0;
      if(barycentric(2) <= barycentric(1) && barycentric(2) <= barycentric(0))
         lArg = (EIGEN_PI/3.0)*(1.0+(barycentric(1)-barycentric(0))/(1-3.0*barycentric(2)));
      if(barycentric(0) <= barycentric(1) && barycentric(0) <= barycentric(2))
         lArg = (EIGEN_PI/3.0)*(3.0+(barycentric(2)-barycentric(1))/(1-3.0*barycentric(0)));
      else
         lArg = (EIGEN_PI/3.0)*(5.0+(barycentric(0)-barycentric(2))/(1-3.0*barycentric(1)));

      // Adjust texture coordinates
      theta += lArg*nu;

      double t1 = 1.0+exp(s*(cos(f*(theta))-w));
      double dt1dtheta = exp(s*(cos(f*(theta))-w))*s*-sin(f*theta)*f;
      double result = 1.0/t1-0.5;
      double dresultdt1 = -1.0/(t1*t1);

      double dresultdtheta = dresultdt1*dt1dtheta;

      if(k == 0)
      {
         double dthetadtx = 1;
         double dthetadnu = lArg;

         for(int i=0; i<2; ++i)
         {
            for(int j=0; j<3; ++j)
               dcolorsdx[i](num, indexI) += dresultdtheta*dthetadtx*dtxdtexture(j)*vertices[indexI].dtexturecoordnx[i](3*triangle+j); 
            dcolorsdx[i](num, indexI) += dresultdtheta*dthetadnu*vertices[indexI].dndnx[i](triangle);

            for(int j=0; j<3; ++j)
               dcolorsdx[i](num, indexJ) += dresultdtheta*dthetadtx*dtxdtexture(j)*vertices[indexJ].dtexturecoordnx[i](3*triangle+j); 
            dcolorsdx[i](num, indexJ) += dresultdtheta*dthetadnu*vertices[indexJ].dndnx[i](triangle);

            for(int j=0; j<3; ++j)
               dcolorsdx[i](num, indexK) += dresultdtheta*dthetadtx*dtxdtexture(j)*vertices[indexK].dtexturecoordnx[i](3*triangle+j); 
            dcolorsdx[i](num, indexK) += dresultdtheta*dthetadnu*vertices[indexK].dndnx[i](triangle);
         }
      }
      else
      {
         for(int i=0; i<2; ++i)
         {
            dcolorsdx[i](num, indexI) = dresultdtheta*dtxdx[i](0);
            dcolorsdx[i](num, indexJ) = dresultdtheta*dtxdx[i](1);
            dcolorsdx[i](num, indexK) = dresultdtheta*dtxdx[i](2);
         }
      }

      // std::cout<<num<<" "<<indexI<<" "<<dcolorsdx[0](num, indexI)<<std::endl;
      //    std::cout<<num<<" "<<indexJ<<" "<<dcolorsdx[0](num, indexJ)<<std::endl;
      //    std::cout<<num<<" "<<indexK<<" "<<dcolorsdx[0](num, indexK)<<std::endl;

      return result;
   }

   bool marchingCube(Eigen::Matrix3d& corners, Eigen::Vector3d& values, Eigen::Vector3d& endpts1, Eigen::Vector3d& endpts2)
   {
      Eigen::Vector3d v1,v2,v3;
      v1 = corners.row(0);
      v2 = corners.row(1);
      v3 = corners.row(2);

      double t1,t2,t3;
      t1 = values(0);
      t2 = values(1);
      t3 = values(2);

      if(t1<=0 && t2<=0 && t3<=0) return false;
      if(t1>=0 && t2>=0 && t3>=0) return false;

      if((t1 >=0 && t2>=0 && t3<=0) || (t1 <=0 && t2<=0 && t3>=0))
      {
         double l1 = (abs(t3))/(abs(t1)+abs(t3));
         endpts1 = v1*l1+v3*(1-l1);
         double l2 = (abs(t3))/(abs(t2)+abs(t3));
         endpts2 = v2*l2+v3*(1-l2);
      }
      else if((t2 >=0 && t3>=0 && t1<=0) || (t2 <=0 && t3<=0 && t1>=0))
      {
         double l1 = (abs(t1))/(abs(t1)+abs(t2));
         endpts1 = v2*l1+v1*(1-l1);
         double l2 = (abs(t1))/(abs(t1)+abs(t3));
         endpts2 = v3*l2+v1*(1-l2);
      }
      else if((t3 >=0 && t1>=0 && t2<=0) || (t3 <=0 && t1<=0 && t2>=0))
      {
         double l1 = (abs(t2))/(abs(t1)+abs(t2));
         endpts1 = v1*l1+v2*(1-l1);
         double l2 = (abs(t2))/(abs(t2)+abs(t3));
         endpts2 = v3*l2+v2*(1-l2);
      }
      return true;
   }

   void Mesh::generateConformalMesh()
   {
      Eigen::MatrixXd V,NV;
      Eigen::MatrixXi F,NF;
      std::vector<std::vector<double>> poly_lines;

      //igl::readOBJ(inputFilename,V,F);
      igl::readOBJ("../data/square50.obj",V,F);
      igl::upsample(V,F,NV,NF,1);
      //igl::false_barycentric_subdivision(V,F,NV,NF);
      //igl::loop(V,F,NV,NF,5);
      Eigen::MatrixXd bs(NV.rows(),3);
      Eigen::VectorXd face_ids(NV.rows());
      for(int i=0; i<NV.rows(); ++i)
      {
         Eigen::MatrixXd P = NV.row(i);
         for(int j=0; j<F.rows(); j++)
         {
            Face fa = faces[j];
            // Get the three half edges.
            HalfEdgeCIter hij = fa.he;
            HalfEdgeCIter hjk = hij->next;
            HalfEdgeCIter hkl = hjk->next;

            // Get the three vertices.
            VertexCIter vi = hij->vertex;
            VertexCIter vj = hjk->vertex;
            VertexCIter vk = hkl->vertex;
            
            Eigen::MatrixXd Va = V.row(vi->index);
            Eigen::MatrixXd Vb = V.row(vj->index);
            Eigen::MatrixXd Vc = V.row(vk->index);

            Eigen::MatrixXd l(1,3);
            igl::barycentric_coordinates(P,Va,Vb,Vc,l);
            if(isValid(l(0), l(1), l(2))){
               bs.row(i) = l;
               face_ids(i) = j;
               break;
            }
         }
      }

      colors.resize(NV.rows());
      dcolorsdx.resize(2);
      dcolorsdx[0].setZero(NV.rows(), V.rows());
      dcolorsdx[1].setZero(NV.rows(), V.rows());
      for(int i=0; i<NV.rows(); ++i)
      {
         int tri = face_ids(i);
         Eigen::Vector2d bary(bs(i,0), bs(i,1));
         colors(i) = SampleTexturePointDerivative(i, tri,bary);
      }


      Eigen::MatrixXd dcolorsdnx(NV.rows(),2*V.rows());
      for(int i=0; i<NV.rows(); ++i){
         for( int j=0; j<V.rows(); ++j){
            dcolorsdnx(i,2*j) = dcolorsdx[0](i,j);
            dcolorsdnx(i,2*j+1) = dcolorsdx[1](i,j);
         }
      }
      

      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         if(v->dXdD.size() > 0){
            v->dcolorsdD = dcolorsdnx*v->dnxdD;
            //std::cout<<v->dnxdD<<std::endl;
         }
      }

      // for(int i=0; i<NF.rows(); ++i)
      // {
      //    //std::cout<<NF.size()<<" "<<i<<std::endl;
      //    Eigen::Matrix3d corners;
      //    corners.row(0) = NV.row(NF(i,0));
      //    corners.row(1) = NV.row(NF(i,1));
      //    corners.row(2) = NV.row(NF(i,2));

      //    Eigen::Vector3d values(colors(NF(i,0)), colors(NF(i,1)), colors(NF(i,2)));
      //    Eigen::Vector3d pt1, pt2;
      //    bool exist_polyline = marchingCube(corners,values, pt1, pt2);

      //    if(exist_polyline && (pt1-pt2).norm() > 0.001){
      //       poly_lines.push_back({pt1(0), pt1(1), pt2(0), pt2(1)});
      //    }
      // }

      // std::ofstream myfile;
      // myfile.open("polyline_conformal_mag.txt");
      // myfile<<poly_lines.size()<<std::endl;
      // for(int i=0; i<poly_lines.size(); ++i)
      // {
      //    myfile<<poly_lines[i][0]<<" "<<poly_lines[i][1]<<" "<<poly_lines[i][2]<<" "<<poly_lines[i][3]<<std::endl;
      // }
      // myfile.close();

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

   void Mesh :: assignTextureCoordinatesDerivative( int p )
   // This method computes the final texture coordinates that are actually
   // used by OpenGL to draw the stripes, starting with the solution to
   // the eigenvalue problem.
   {
      int nF = faces.size();
      int nV = vertices.size();
      textureCoor.setZero(3*nF);
      ns.setZero(nF);
      dtexturedthetas.setZero(3*nF,nV);
      for(int i=0; i<vertices.size(); ++i)
      {
         vertices[i].dtexturecoordnx.resize(2);
         vertices[i].dtexturecoordnx[0].setZero(3*nF,1);
         vertices[i].dtexturecoordnx[1].setZero(3*nF,1);
         
         vertices[i].dndnx.resize(2);
         vertices[i].dndnx[0].setZero(nF,1);
         vertices[i].dndnx[1].setZero(nF,1);
      }

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

         int indexI = hij->vertex->index; 
         int indexJ = hjk->vertex->index; 
         int indexK = hki->vertex->index; 

         int revert_d = 1;
         int revert_f = 1;

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
            revert_d = -1;
            omegaIJ =  cIJ * omegaIJ;
            omegaJK = -cJK * omegaJK;
         }

         if( hki->edge->crossesSheets )
         {
            psiK = psiK.bar();
            revert_f = -1;
            omegaKI = -cKI * omegaKI;
            omegaJK =  cJK * omegaJK;
         }

         double thetaI = psiI.arg();
         double thetaJ = psiJ.arg();
         double thetaK = psiK.arg();


         // construct complex transport coefficients
         Complex rij( cos(omegaIJ), sin(omegaIJ) );
         Complex rjk( cos(omegaJK), sin(omegaJK) );
         Complex rki( cos(omegaKI), sin(omegaKI) );

         // compute the angles at the triangle corners closest to the target omegas
         double alphaI = psiI.arg();
         // double alphaJ = alphaI + omegaIJ - (rij*psiI/psiJ).arg(); //fmodPI((varphiI + omegaIJ) - varphiJ); // could do this in terms of angles instead of complex numbers...
         // double alphaK = alphaJ + omegaJK - (rjk*psiJ/psiK).arg(); //fmodPI((varphiJ + omegaJK) - varphiK); // mostly a matter of taste---possibly a matter of performance?
         // double alphaL = alphaK + omegaKI - (rki*psiK/psiI).arg(); //fmodPI((varphiK + omegaKI) - varphiI);

         double de1,de2,de3;

         double alphaJ = alphaI + omegaIJ - fmodPIDiff(thetaI+omegaIJ-thetaJ,de1);
         double alphaK = alphaJ + omegaJK - fmodPIDiff(thetaJ+omegaJK-thetaK,de2);
         double alphaL = alphaK + omegaKI - fmodPIDiff(thetaK+omegaKI-thetaI,de3);

         // checkfmodPIDiff(thetaI+omegaIJ-thetaJ);
         // checkfmodPIDiff(thetaJ+omegaJK-thetaK);
         // checkfmodPIDiff(thetaK+omegaKI-thetaI);

         //std::cout<<omegaIJ<<" "<<omegaJK<<" "<<omegaKI<<std::endl;

         Eigen::MatrixXd DalphasDthetas(4,3);
         DalphasDthetas.setZero();
         DalphasDthetas(0,0) = 1.0;
         DalphasDthetas(1,0) = DalphasDthetas(0,0) - de1;
         DalphasDthetas(1,1) = de1;
         DalphasDthetas(2,0) = DalphasDthetas(1,0);
         DalphasDthetas(2,1) = DalphasDthetas(1,1) - de2;
         DalphasDthetas(2,2) = de2;
         DalphasDthetas(3,0) = DalphasDthetas(2,0) + de3;
         DalphasDthetas(3,1) = DalphasDthetas(2,1);
         DalphasDthetas(3,2) = DalphasDthetas(2,2) - de3;



         // adjust triangles containing zeros
         //double n = lround((alphaL-alphaI)/(2.*M_PI));
         
         // Differentiable left round
         double dndL = 0;
         double L = (alphaL-alphaI)/(2.*M_PI);
         double n = lroundDiff(L,dndL);

         Eigen::MatrixXd DnDthetas(1,3);
         DnDthetas.setZero();
         for(int i=0; i<3; ++i)
            DnDthetas(i) = -dndL/(2.*M_PI)*DalphasDthetas(0,i)+dndL/(2.*M_PI)*DalphasDthetas(3,i);

         alphaJ -= 2.*M_PI*n/3.;
         alphaK -= 4.*M_PI*n/3.;

         for(int i=0; i<3; ++i){
            DalphasDthetas(1,i) -= 2.*M_PI*DnDthetas(i)/3.;
            DalphasDthetas(2,i) -= 4.*M_PI*DnDthetas(i)/3.;
         }

         // store the coordinates
         hij->texcoord[p] = alphaI;
         hjk->texcoord[p] = alphaJ;
         hki->texcoord[p] = alphaK;
         f->paramIndex[p] = n;

         int i = f->index;
         textureCoor(3*i) = alphaI;
         textureCoor(3*i+1) = alphaJ;
         textureCoor(3*i+2) = alphaK;
         ns(i) = n;

         //if(i==0) std::cout<<omegaIJ<<" "<<omegaJK<<" "<<omegaKI<<std::endl;

         

         for(int a=0; a<3; ++a){
            vertices[indexI].dtexturecoordnx[0](3*i+a) = -psiI.im/psiI.norm2()*DalphasDthetas(a,0);
            vertices[indexI].dtexturecoordnx[1](3*i+a) = psiI.re/psiI.norm2()*DalphasDthetas(a,0);
            vertices[indexJ].dtexturecoordnx[0](3*i+a) = -psiJ.im/psiJ.norm2()*DalphasDthetas(a,1);
            vertices[indexJ].dtexturecoordnx[1](3*i+a) = psiJ.re/psiJ.norm2()*DalphasDthetas(a,1)*revert_d;
            vertices[indexK].dtexturecoordnx[0](3*i+a) = -psiK.im/psiK.norm2()*DalphasDthetas(a,2);
            vertices[indexK].dtexturecoordnx[1](3*i+a) = psiK.re/psiK.norm2()*DalphasDthetas(a,2)*revert_f;

            dtexturedthetas(3*i+a,indexI) = DalphasDthetas(a,0);
            dtexturedthetas(3*i+a,indexJ) = DalphasDthetas(a,1)*revert_d;
            dtexturedthetas(3*i+a,indexK) = DalphasDthetas(a,2)*revert_f;
         }
         vertices[indexI].dndnx[0](i) = -psiI.im*DnDthetas(0);
         vertices[indexI].dndnx[1](i) = psiI.re*DnDthetas(0);
         vertices[indexJ].dndnx[0](i) = -psiJ.im*DnDthetas(1);
         vertices[indexJ].dndnx[1](i) = psiJ.re*DnDthetas(1)*revert_d;
         vertices[indexK].dndnx[0](i) = -psiK.im*DnDthetas(2);
         vertices[indexK].dndnx[1](i) = psiK.re*DnDthetas(2)*revert_f;

         //std::cout<<nF<<" "<<i<<std::endl;
         
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

