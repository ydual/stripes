#include "Face.h"
#include "Mesh.h"
#include "Vector.h"
#include "Utility.h"

namespace DDG
{
   Vector Face::normal( void ) const
   {
      Vector p0 = he->vertex->position;
      Vector p1 = he->next->vertex->position;
      Vector p2 = he->next->next->vertex->position;

      return cross( p1-p0, p2-p0 ).unit();
   }

   double Face::area( void ) const
   {
      Vector p0 = he->vertex->position;
      Vector p1 = he->next->vertex->position;
      Vector p2 = he->next->next->vertex->position;

      return cross( p1-p0, p2-p0 ).norm() / 2.;
   }

   bool Face::isBoundary( void ) const
   {
      return he->onBoundary;
   }

   double Face :: curvature( void ) const
   {
      double Omega = 0.;

      HalfEdgeCIter h = he;
      do
      {
         double thetaI = h->angularCoordinate;
         double thetaJ = h->flip->angularCoordinate + M_PI;

         Omega += thetaJ - thetaI;

         h = h->next;
      }
      while( h != he );

      return fmodPI( Omega );
   }

   double fmodPIDiff2(double theta, double& derivative)
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

   double lroundDiff2(double x, double& derivative)
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

   double Face :: fieldIndex( double k ) const
   {
      if( isBoundary() ) return 0.;

      double Omega = 0.;
      double idx = 0.;
      HalfEdgeCIter hi = he;
      do
      {
         double phiI = hi->vertex->directionField.arg();
         double phiJ = hi->flip->vertex->directionField.arg();
         double thetaI = hi->angularCoordinate;
         double thetaJ = hi->flip->angularCoordinate + M_PI;
         double dTheta = thetaI - thetaJ;

         Omega += dTheta;
         idx += fmodPI( phiJ - phiI + k*dTheta );

         hi = hi->next;
      }
      while( hi != he );

      idx -= k * fmodPI(Omega);

      return lround( idx / (2.*M_PI) );
      //return 0;
   }

   // double Face :: paramIndex( void ) const
   // {
   //    // TODO
   //    HalfEdgeCIter hij = he;
   //    HalfEdgeCIter hjk = hij->next;
   //    HalfEdgeCIter hkl = hjk->next;

   //    double omegaIJ = hij->edge->omega;
   //    double omegaJK = hjk->edge->omega;
   //    double omegaKL = hkl->edge->omega;

   //    Complex phiI = hij->vertex->parameterization;
   //    Complex phiJ = hjk->vertex->parameterization;
   //    Complex phiK = hkl->vertex->parameterization;
   //    Complex phiL = phiI;

   //    // TODO clean up the logic here...
   //    if( hij->edge->crossesSheets ) { phiJ = phiJ.bar(); phiK = phiK.bar(); omegaJK = -omegaJK; phiL = phiL.bar(); omegaKL = -omegaKL; } else if( hij->edge->he == hij ) { omegaIJ = -omegaIJ; }
   //    if( hjk->edge->crossesSheets ) { phiK = phiK.bar(); phiL = phiL.bar(); omegaKL = -omegaKL; } else if( hjk->edge->he == hjk ) { omegaJK = -omegaJK; }
   //    if( hkl->edge->crossesSheets ) { phiL = phiL.bar(); } else if( hkl->edge->he == hkl ) { omegaKL = -omegaKL; }

   //    Complex rij( cos(omegaIJ), sin(omegaIJ) );
   //    Complex rjk( cos(omegaJK), sin(omegaJK) );
   //    Complex rkl( cos(omegaKL), sin(omegaKL) );

   //    double sigmaIJ = omegaIJ + ((rij*phiI)*(phiJ.inv())).arg();
   //    double sigmaJK = omegaJK + ((rjk*phiJ)*(phiK.inv())).arg();
   //    double sigmaKL = omegaKL + ((rkl*phiK)*(phiL.inv())).arg();
   //    double sigmaIJK = sigmaIJ + sigmaJK + sigmaKL;

   //    return lround( sigmaIJK / (2.*M_PI) );
   // }

   Vector Face :: barycenter( void ) const
   {
      return ( he->vertex->position +
               he->next->vertex->position +
               he->next->next->vertex->position ) / 3.;
   }

   void Face :: orientTexCoords( void )
   {
      Complex& p0 = he->texcoord;
      Complex& p1 = he->next->texcoord;
      Complex& p2 = he->next->next->texcoord;

      Complex u = p1-p0;
      Complex v = p2-p0;

      if( cross( u, v ) < 0. )
      {
         p0 = p0.bar();
         p1 = p1.bar();
         p2 = p2.bar();
      }
   }

   void Face::getLocalSheet( std::vector<Complex>& psi, std::vector<double>& omega )
   {
      psi.resize( 3 );
      omega.resize( 3 );

      HalfEdgeCIter hij = he;
      HalfEdgeCIter hjk = hij->next;
      HalfEdgeCIter hki = hjk->next;

      Complex psiI = hij->vertex->parameterization;
      Complex psiJ = hjk->vertex->parameterization;
      Complex psiK = hki->vertex->parameterization;

      double omegaIJ = hij->omega();
      double omegaJK = hjk->omega();
      double omegaKI = hki->omega();

      if( hij->edge->crossesSheets )
      {
         psiJ = -psiJ;

         if( hij->edge->he != hij ) omegaIJ = -omegaIJ;
         if( hjk->edge->he == hjk ) omegaJK = -omegaJK;
      }

      if( hjk->edge->crossesSheets )
      {
         psiK = -psiK;

         if( hki->edge->he == hki ) omegaKI = -omegaKI;
         if( hjk->edge->he != hjk ) omegaJK = -omegaJK;
      }

      psi[0] = psiI;
      psi[1] = psiJ;
      psi[2] = psiK;

      omega[0] = omegaIJ;
      omega[1] = omegaJK;
      omega[2] = omegaKI;
   }
}

