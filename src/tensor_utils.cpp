#include "tensor_utils.hpp"

TensorUtils::TensorUtils() : ad_helper(dim * dim) {
      P_flat = Vector<double>(dim * dim);
      C_flat = FullMatrix<double>(dim * dim, dim * dim);
    }

void TensorUtils::compute_tensors(const Tensor<2, dim> &F, Tensor<2, dim> &P,
                         Tensor<4, dim> &C) {
      const FEValuesExtractors::Tensor<2> F_linearized(0);
      const FEValuesExtractors::Tensor<4> C_linearized(0);
      const typename AutoDiff::Types<ADNumberType>::tape_index tape_index = 0;

      const bool is_recording =
          ad_helper.start_recording_operations(tape_index);
      if (is_recording) {
        ad_helper.register_independent_variable(F, F_linearized);
        ADNumberType W_ad = compute_W(F);
        ad_helper.register_dependent_variable(W_ad);
        ad_helper.stop_recording_operations(false);
      } else {
        ad_helper.activate_recorded_tape(tape_index);
        ad_helper.set_independent_variable(F, F_linearized);
      }

      ad_helper.compute_gradient(P_flat);
      ad_helper.compute_hessian(C_flat);

      P = ad_helper.extract_gradient_component(P_flat, F_linearized);
      C = ad_helper.extract_hessian_component(C_flat, F_linearized,
                                              F_linearized);
    }

typename TensorUtils::ADNumberType TensorUtils::compute_W(const Tensor<2, dim, ADNumberType> &F) const {
      const ADNumberType J_ad = determinant(F);
      const ADTensor2 C_ad = transpose(F) * F;
      const ADNumberType I1_ad = trace(C_ad);
      ADNumberType psi_ad =
          (mu_hook / 2.0) * (I1_ad * std::pow(J_ad, -2.0 / 3.0) - 3.0);
      psi_ad += (k_hook / 2.0) * std::pow(J_ad - 1.0, 2.0);

      return psi_ad;
    }
