#include "Stripes.h"
#include <stripes/Mesh.h>
#include <stripes/MeshIO.h>

namespace DDG {

////////////////////////////////////////////////////////////////////////////////

namespace {

///
/// Project a vector onto a given plane.
///
/// @param[in]  vec     The vector to project onto the plane.
/// @param[in]  normal  Normal to the plane to project onto.
/// @param[in]  axis    Axis vector whose projection onto the plane defines the local X axis.
///
/// @return     Projected vector.
///
Eigen::Vector2d 
projectOntoPlane(const Eigen::Vector3d &vec,
                                 const Eigen::Vector3d &normal,
                                 const Eigen::Vector3d &axis)
{
    Eigen::Vector3d eY = normal.cross(axis).normalized();
    Eigen::Vector3d eX = eY.cross(normal).normalized();
    return {eX.dot(vec), eY.dot(vec)};
}

// -----------------------------------------------------------------------------

class MeshEigen : public Mesh {
public:
    int init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
    void setDirectionField(const Eigen::MatrixXd &D, bool useSmoothestField, double theta);
    void setDirectionFieldDerivative(const Eigen::MatrixXd &D, bool useSmoothestField, double theta);

};

// Return 0 for success
int MeshEigen::init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    MeshData meshData;
    meshData.positions.reserve(V.rows());
    for (int i = 0; i < V.rows(); ++i) {
        meshData.positions.emplace_back(V(i, 0), V(i, 1), V(i, 2));
    }
    meshData.indices.reserve(F.rows());
    for (int i = 0; i < F.rows(); ++i) {
        meshData.indices.emplace_back();
        meshData.indices.back().reserve(F.cols());
        for (int lv = 0; lv < F.cols(); ++lv) {
            meshData.indices.back().emplace_back(F(i, lv), -1, -1);
        }
    }
    if (auto ret = MeshIO::buildMesh(meshData, *this)) {
        return ret;
    }
    indexElements();
    initializeSmoothStructure();
    buildMassMatrix();

    return 0;
}
void MeshEigen::setDirectionField(const Eigen::MatrixXd &D, bool useSmoothestField, double theta)
{
    if ((size_t)D.rows() == vertices.size()) {
        for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
            double alpha = (*v->he)->angularCoordinate;
            Complex r(cos(2. * alpha), sin(2. * alpha));

            auto i = v->index;
            auto n = v->normal();
            auto e = (*v->he)->vector();
            auto u = projectOntoPlane(D.row(i).transpose(), {n.x, n.y, n.z}, {e.x, e.y, e.z});
            //std::cout<<e.x<<" "<<e.y<<" "<<e.z<<std::endl;
            double a = std::atan2(u.y(), u.x()) + M_PI / 2.0;
            v->directionField = r * Complex(cos(2.0 * a), sin(2.0 * a));
        }
    }
    else {
        if (useSmoothestField) {
            computeSmoothestSection();
        }
        else {
            computeCurvatureAlignedSection();
        }
    }

    // rotate 1-vector field by θ radians, which in our complex representation
    // is equivalent to rotating a 2-vector field by 2*θ radians
    Complex z(cos(2.0 * theta), sin(2.0 * theta));
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->directionField *= z;
    }
}

void MeshEigen::setDirectionFieldDerivative(const Eigen::MatrixXd &D, bool useSmoothestField, double theta)
{
    if ((size_t)D.rows() == vertices.size()) {
        for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
            double alpha = (*v->he)->angularCoordinate;
            Complex r(cos(2. * alpha), sin(2. * alpha));

            auto i = v->index;
            auto n = v->normal();
            auto e = (*v->he)->vector();
            auto u = projectOntoPlane(D.row(i).transpose(), {n.x, n.y, n.z}, {e.x, e.y, e.z});
            double l = u.y() / u.x();
            double a = std::atan2(u.y(), u.x()) + M_PI / 2.0;
            v->directionField = r * Complex(cos(2.0 * a), sin(2.0 * a));


            Eigen::Vector3d normal = {n.x, n.y, n.z};
            Eigen::Vector3d axis = {e.x, e.y, e.z};
            Eigen::Vector3d eY = normal.cross(axis).normalized();
            Eigen::Vector3d eX = eY.cross(normal).normalized();

            double dreda = -2.0*cos(2.0*alpha)*sin(2.0*a)-2*sin(2.0*alpha)*cos(2.0*a);
            double dimda = -2.0*sin(2.0*alpha)*sin(2.0*a)+2*cos(2.0*alpha)*cos(2.0*a);
            double dadl = 1./(1+l*l);
            double dldux = -u.y()/(u.x()*u.x());
            double dlduy = 1./u.x();
            Eigen::VectorXd duxdD = eX;
            Eigen::VectorXd duydD = eY;

            v->dXdD.resize(2);
            v->dXdD[0] = dreda*dadl*(dldux*duxdD+dlduy*duydD);
            v->dXdD[1] = dimda*dadl*(dldux*duxdD+dlduy*duydD);

            //std::cout<<"fuck: "<<v->dXdD[0].rows()<<" "<<v->dXdD[0].cols()<<std::endl;
        }
    }
    else {
        if (useSmoothestField) {
            std::cerr<<"Not Implemented!"<<std::endl;
        }
        else {
            std::cerr<<"Not Implemented!"<<std::endl;
        }
    }
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

int computeStripePatterns(
                          bool is_target,
                          Eigen::VectorXd& v_tar,
                          std::vector<Eigen::MatrixXd>& dpdD,
                          Eigen::VectorXd &param,
                          const Eigen::MatrixXd &V,
                          Eigen::MatrixXi &F,
                          Eigen::MatrixXd &D,
                          double theta,
                          double frequency,
                          int numCoords,
                          Eigen::VectorXi &branchIndex,
                          Eigen::MatrixXd &parameterization,
                          Eigen::MatrixXi &zeroIndex,
                          std::vector<bool> &isBorder,
                          bool useSmoothestField,
                          std::string *error)
{
    assert(V.cols() == 3);
    assert(F.cols() == 3);
    if (!(numCoords == 1 || numCoords == 2)) {
        throw std::runtime_error("Invalid number of coordinate functions.");
    }

    // Compute field-aligned parametrization
    MeshEigen mesh;
    if (auto ret = mesh.init(V, F)) {
        return ret;
    }
    if(!is_target){
        mesh.target = false;
        mesh.v_tar = v_tar;
    } 
    mesh.setDirectionField(D, useSmoothestField, theta);
    mesh.setDirectionFieldDerivative(D, useSmoothestField, theta);

    // double eps = 1e-5;
    // for(VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
    // {
    //     auto i = v->index;
    //     double X_o = v->directionField.im;

    //     D(i,1) += eps;
    //     mesh.setDirectionField(D, useSmoothestField, theta);
    //     D(i,1) -= eps;

    //     double X_p = v->directionField.im;
    //     mesh.setDirectionField(D, useSmoothestField, theta);

    //     double numerical = (X_p-X_o)/eps;
    //     std::cout<<i<<" numerical: "<<numerical<<" analytical: "<<v->dXdD[1](1)<<" "<<numerical/(v->dXdD[1](1))<<std::endl;
        
    // }


    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        auto z = v->directionField;
        if (!(std::isfinite(z.re) && std::isfinite(z.im))) {
            if (error) {
                *error =
                    "Direction field has NaN or infinite values. "
                    "Make sure the input mesh is clean of isolated and non-manifold vertices.";
            }
            return -1;
        }
    }
    mesh.nCoordinateFunctions = 1;
    mesh.lambda = frequency;
    mesh.parameterize();
    if (numCoords == 2) {
        mesh.glueParameterization();
    }

    // Convert result
    branchIndex.resize(F.rows());
    branchIndex.setZero();
    parameterization.resize(3 * F.rows(), numCoords);
    parameterization.setZero();
    zeroIndex.resize(F.rows(), numCoords);
    zeroIndex.setZero();
    isBorder.resize(F.rows());

    int faceIndex = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++, ++faceIndex) {
        isBorder[faceIndex] = f->isBoundary();
        if (f->isBoundary()) continue;

        double k = f->fieldIndex(2.);
        branchIndex(faceIndex) = (int)k;
        for (int p = 0; p < numCoords; ++p) {
            zeroIndex(faceIndex, p) = (int)f->paramIndex[p];
        }

        // We don't do anything special for faces with singularities for now
        int i = 0;
        HalfEdgeCIter he = f->he;
        do {
            Complex g = he->texcoord;
            F(faceIndex, i) = he->vertex->index;
            parameterization(3 * faceIndex + i, 0) = g.re;
            if (numCoords == 2) {
                parameterization(3 * faceIndex + i, 1) = g.im;
            }
            i++;
            he = he->next;
        } while (he != f->he);
    }
    param.setZero();
    dpdD.resize(V.rows());
    int index = 0;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        dpdD[index] = v->dcolorsdD;
        index++;
    }
    param = mesh.colors;

    if(is_target) v_tar = mesh.v_tar;

    return 0;
}

}  // namespace DDG
