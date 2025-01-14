#pragma once

#include <eigen3/Eigen/Dense>
#include <string>
#include <vector>

namespace DDG {

///
/// Compute a stripe pattern aligned to an input direction field.
///
/// @param[in]     V                        #V x 3 input mesh vertices.
/// @param[in,out] F                        #F x 3 input mesh faces. Facet may be reordered by the
///                                         algorithm.
/// @param[in]     D                        #V x 3 input direction field per vertex. If the matrix
///                                         is empty, then either the smoothest field or the
///                                         curvature field is used, based on the value of the input
///                                         flag `useSmoothestField`.
/// @param[in]     theta                    Rotate input field by a given angle (between 0 and 2*pi).
/// @param[in]     frequency                Frequency of output stripes.
/// @param[in]     numCoordinateFunctions   Number of parametrization to compute (D=1 or 2).
/// @param[out]    branchIndex              #F x 1 branch index of each parametrization. A zero
///                                         indicates a regular face.
/// @param[out]    parameterization         3*#F x D parametrization at each face corner.
/// @param[out]    zeroIndex                #F x D zero index of each parametrization.
/// @param[out]    isBorder                 #F x 1 array indicating if a face is on the border of
///                                         the mesh.
/// @param[in]     useSmoothestField        If no input direction field is provided, computes either
///                                         the smoothest field, or the field of principal
///                                         curvature, based on the value of this flag.
/// @param[out]    error                    Error message in case of failure. direction field.
///
/// @return        0 in case of success.
///
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
                          int numCoordinateFunctions,
                          Eigen::VectorXi &branchIndex,
                          Eigen::MatrixXd &parameterization,
                          Eigen::MatrixXi &zeroIndex,
                          std::vector<bool> &isBorder,
                          bool useSmoothestField = false,
                          std::string *error = nullptr);

}  // namespace DDG
