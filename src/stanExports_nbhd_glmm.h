// Generated by rstantools.  Do not edit by hand.

#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_nbhd_glmm_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 57> locations_array__ =
  {" (found before start of program)",
  " (in 'nbhd_glmm', line 37, column 4 to column 18)",
  " (in 'nbhd_glmm', line 38, column 4 to column 19)",
  " (in 'nbhd_glmm', line 39, column 4 to column 23)",
  " (in 'nbhd_glmm', line 40, column 4 to column 21)",
  " (in 'nbhd_glmm', line 56, column 4 to column 31)",
  " (in 'nbhd_glmm', line 58, column 4 to column 59)",
  " (in 'nbhd_glmm', line 44, column 15 to column 16)",
  " (in 'nbhd_glmm', line 44, column 8 to column 63)",
  " (in 'nbhd_glmm', line 45, column 8 to column 29)",
  " (in 'nbhd_glmm', line 46, column 8 to column 26)",
  " (in 'nbhd_glmm', line 43, column 21 to line 47, column 5)",
  " (in 'nbhd_glmm', line 43, column 4 to line 47, column 5)",
  " (in 'nbhd_glmm', line 48, column 4 to column 37)",
  " (in 'nbhd_glmm', line 49, column 4 to column 51)",
  " (in 'nbhd_glmm', line 51, column 4 to column 27)",
  " (in 'nbhd_glmm', line 52, column 4 to column 23)",
  " (in 'nbhd_glmm', line 8, column 4 to column 19)",
  " (in 'nbhd_glmm', line 9, column 10 to column 11)",
  " (in 'nbhd_glmm', line 9, column 4 to column 13)",
  " (in 'nbhd_glmm', line 10, column 4 to column 19)",
  " (in 'nbhd_glmm', line 11, column 11 to column 12)",
  " (in 'nbhd_glmm', line 11, column 14 to column 15)",
  " (in 'nbhd_glmm', line 11, column 4 to column 19)",
  " (in 'nbhd_glmm', line 13, column 4 to column 22)",
  " (in 'nbhd_glmm', line 14, column 20 to column 21)",
  " (in 'nbhd_glmm', line 14, column 4 to column 23)",
  " (in 'nbhd_glmm', line 15, column 4 to column 24)",
  " (in 'nbhd_glmm', line 16, column 4 to column 29)",
  " (in 'nbhd_glmm', line 17, column 4 to column 19)",
  " (in 'nbhd_glmm', line 20, column 4 to column 19)",
  " (in 'nbhd_glmm', line 21, column 11 to column 12)",
  " (in 'nbhd_glmm', line 21, column 14 to column 16)",
  " (in 'nbhd_glmm', line 21, column 4 to column 21)",
  " (in 'nbhd_glmm', line 22, column 11 to column 13)",
  " (in 'nbhd_glmm', line 22, column 4 to column 23)",
  " (in 'nbhd_glmm', line 24, column 11 to column 12)",
  " (in 'nbhd_glmm', line 24, column 14 to column 16)",
  " (in 'nbhd_glmm', line 24, column 4 to column 21)",
  " (in 'nbhd_glmm', line 25, column 11 to column 13)",
  " (in 'nbhd_glmm', line 25, column 15 to column 17)",
  " (in 'nbhd_glmm', line 25, column 4 to column 22)",
  " (in 'nbhd_glmm', line 26, column 11 to column 13)",
  " (in 'nbhd_glmm', line 26, column 15 to column 17)",
  " (in 'nbhd_glmm', line 26, column 4 to column 26)",
  " (in 'nbhd_glmm', line 28, column 8 to column 38)",
  " (in 'nbhd_glmm', line 29, column 8 to column 46)",
  " (in 'nbhd_glmm', line 27, column 19 to line 30, column 5)",
  " (in 'nbhd_glmm', line 27, column 4 to line 30, column 5)",
  " (in 'nbhd_glmm', line 32, column 4 to column 37)",
  " (in 'nbhd_glmm', line 33, column 4 to column 37)",
  " (in 'nbhd_glmm', line 34, column 4 to column 25)",
  " (in 'nbhd_glmm', line 37, column 11 to column 13)",
  " (in 'nbhd_glmm', line 40, column 11 to column 15)",
  " (in 'nbhd_glmm', line 56, column 11 to column 13)",
  " (in 'nbhd_glmm', line 4, column 8 to column 32)",
  " (in 'nbhd_glmm', line 3, column 29 to line 5, column 5)"};
template <typename T0__,
          stan::require_all_t<stan::is_col_vector<T0__>,
                              stan::is_vt_not_complex<T0__>>* = nullptr>
Eigen::Matrix<stan::promote_args_t<stan::base_type_t<T0__>>,-1,1>
cloglog(const T0__& p_arg__, std::ostream* pstream__);
template <typename T0__,
          stan::require_all_t<stan::is_col_vector<T0__>,
                              stan::is_vt_not_complex<T0__>>*>
Eigen::Matrix<stan::promote_args_t<stan::base_type_t<T0__>>,-1,1>
cloglog(const T0__& p_arg__, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::base_type_t<T0__>>;
  int current_statement__ = 0;
  const auto& p = stan::math::to_ref(p_arg__);
  static constexpr bool propto__ = true;
  // suppress unused var warning
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  // suppress unused var warning
  (void) DUMMY_VAR__;
  try {
    current_statement__ = 55;
    return stan::math::log(
             stan::math::minus(stan::math::log(stan::math::subtract(1, p))));
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
  }
}
#include <stan_meta_header.hpp>
class model_nbhd_glmm final : public model_base_crtp<model_nbhd_glmm> {
private:
  int N;
  std::vector<int> Y;
  int K;
  Eigen::Matrix<double,-1,-1> X_data__;
  int N_id;
  std::vector<int> id;
  double bQ_prior_scale;
  int bQ_prior_df;
  int prior_only;
  int Kc;
  Eigen::Matrix<double,-1,-1> Xc_data__;
  Eigen::Matrix<double,-1,1> means_X_data__;
  Eigen::Matrix<double,-1,-1> XQ_data__;
  Eigen::Matrix<double,-1,-1> XR_data__;
  Eigen::Matrix<double,-1,-1> XR_inv_data__;
  Eigen::Map<Eigen::Matrix<double,-1,-1>> X{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> Xc{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,1>> means_X{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> XQ{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> XR{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> XR_inv{nullptr, 0, 0};
public:
  ~model_nbhd_glmm() {}
  model_nbhd_glmm(stan::io::var_context& context__, unsigned int
                  random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_nbhd_glmm_namespace::model_nbhd_glmm";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 17;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      current_statement__ = 17;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 17;
      stan::math::check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 18;
      stan::math::validate_non_negative_index("Y", "N", N);
      current_statement__ = 19;
      context__.validate_dims("data initialization", "Y", "int",
        std::vector<size_t>{static_cast<size_t>(N)});
      Y = std::vector<int>(N, std::numeric_limits<int>::min());
      current_statement__ = 19;
      Y = context__.vals_i("Y");
      current_statement__ = 20;
      context__.validate_dims("data initialization", "K", "int",
        std::vector<size_t>{});
      K = std::numeric_limits<int>::min();
      current_statement__ = 20;
      K = context__.vals_i("K")[(1 - 1)];
      current_statement__ = 20;
      stan::math::check_greater_or_equal(function__, "K", K, 1);
      current_statement__ = 21;
      stan::math::validate_non_negative_index("X", "N", N);
      current_statement__ = 22;
      stan::math::validate_non_negative_index("X", "K", K);
      current_statement__ = 23;
      context__.validate_dims("data initialization", "X", "double",
        std::vector<size_t>{static_cast<size_t>(N), static_cast<size_t>(K)});
      X_data__ = Eigen::Matrix<double,-1,-1>::Constant(N, K,
                   std::numeric_limits<double>::quiet_NaN());
      new (&X) Eigen::Map<Eigen::Matrix<double,-1,-1>>(X_data__.data(), N, K);
      {
        std::vector<local_scalar_t__> X_flat__;
        current_statement__ = 23;
        X_flat__ = context__.vals_r("X");
        current_statement__ = 23;
        pos__ = 1;
        current_statement__ = 23;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 23;
          for (int sym2__ = 1; sym2__ <= N; ++sym2__) {
            current_statement__ = 23;
            stan::model::assign(X, X_flat__[(pos__ - 1)],
              "assigning variable X", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 23;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 24;
      context__.validate_dims("data initialization", "N_id", "int",
        std::vector<size_t>{});
      N_id = std::numeric_limits<int>::min();
      current_statement__ = 24;
      N_id = context__.vals_i("N_id")[(1 - 1)];
      current_statement__ = 24;
      stan::math::check_greater_or_equal(function__, "N_id", N_id, 1);
      current_statement__ = 25;
      stan::math::validate_non_negative_index("id", "N", N);
      current_statement__ = 26;
      context__.validate_dims("data initialization", "id", "int",
        std::vector<size_t>{static_cast<size_t>(N)});
      id = std::vector<int>(N, std::numeric_limits<int>::min());
      current_statement__ = 26;
      id = context__.vals_i("id");
      current_statement__ = 26;
      stan::math::check_greater_or_equal(function__, "id", id, 1);
      current_statement__ = 27;
      context__.validate_dims("data initialization", "bQ_prior_scale",
        "double", std::vector<size_t>{});
      bQ_prior_scale = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 27;
      bQ_prior_scale = context__.vals_r("bQ_prior_scale")[(1 - 1)];
      current_statement__ = 28;
      context__.validate_dims("data initialization", "bQ_prior_df", "int",
        std::vector<size_t>{});
      bQ_prior_df = std::numeric_limits<int>::min();
      current_statement__ = 28;
      bQ_prior_df = context__.vals_i("bQ_prior_df")[(1 - 1)];
      current_statement__ = 28;
      stan::math::check_greater_or_equal(function__, "bQ_prior_df",
        bQ_prior_df, 1);
      current_statement__ = 29;
      context__.validate_dims("data initialization", "prior_only", "int",
        std::vector<size_t>{});
      prior_only = std::numeric_limits<int>::min();
      current_statement__ = 29;
      prior_only = context__.vals_i("prior_only")[(1 - 1)];
      current_statement__ = 30;
      Kc = std::numeric_limits<int>::min();
      current_statement__ = 30;
      Kc = (K - 1);
      current_statement__ = 31;
      stan::math::validate_non_negative_index("Xc", "N", N);
      current_statement__ = 32;
      stan::math::validate_non_negative_index("Xc", "Kc", Kc);
      current_statement__ = 33;
      Xc_data__ = Eigen::Matrix<double,-1,-1>::Constant(N, Kc,
                    std::numeric_limits<double>::quiet_NaN());
      new (&Xc) Eigen::Map<Eigen::Matrix<double,-1,-1>>(Xc_data__.data(), N,
        Kc);
      current_statement__ = 34;
      stan::math::validate_non_negative_index("means_X", "Kc", Kc);
      current_statement__ = 35;
      means_X_data__ = Eigen::Matrix<double,-1,1>::Constant(Kc,
                         std::numeric_limits<double>::quiet_NaN());
      new (&means_X)
        Eigen::Map<Eigen::Matrix<double,-1,1>>(means_X_data__.data(), Kc);
      current_statement__ = 36;
      stan::math::validate_non_negative_index("XQ", "N", N);
      current_statement__ = 37;
      stan::math::validate_non_negative_index("XQ", "Kc", Kc);
      current_statement__ = 38;
      XQ_data__ = Eigen::Matrix<double,-1,-1>::Constant(N, Kc,
                    std::numeric_limits<double>::quiet_NaN());
      new (&XQ) Eigen::Map<Eigen::Matrix<double,-1,-1>>(XQ_data__.data(), N,
        Kc);
      current_statement__ = 39;
      stan::math::validate_non_negative_index("XR", "Kc", Kc);
      current_statement__ = 40;
      stan::math::validate_non_negative_index("XR", "Kc", Kc);
      current_statement__ = 41;
      XR_data__ = Eigen::Matrix<double,-1,-1>::Constant(Kc, Kc,
                    std::numeric_limits<double>::quiet_NaN());
      new (&XR) Eigen::Map<Eigen::Matrix<double,-1,-1>>(XR_data__.data(), Kc,
        Kc);
      current_statement__ = 42;
      stan::math::validate_non_negative_index("XR_inv", "Kc", Kc);
      current_statement__ = 43;
      stan::math::validate_non_negative_index("XR_inv", "Kc", Kc);
      current_statement__ = 44;
      XR_inv_data__ = Eigen::Matrix<double,-1,-1>::Constant(Kc, Kc,
                        std::numeric_limits<double>::quiet_NaN());
      new (&XR_inv)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(XR_inv_data__.data(), Kc, Kc);
      current_statement__ = 48;
      for (int i = 2; i <= K; ++i) {
        current_statement__ = 45;
        stan::model::assign(means_X,
          stan::math::mean(
            stan::model::rvalue(X, "X", stan::model::index_omni(),
              stan::model::index_uni(i))), "assigning variable means_X",
          stan::model::index_uni((i - 1)));
        current_statement__ = 46;
        stan::model::assign(Xc,
          stan::math::subtract(
            stan::model::rvalue(X, "X", stan::model::index_omni(),
              stan::model::index_uni(i)),
            stan::model::rvalue(means_X, "means_X",
              stan::model::index_uni((i - 1)))), "assigning variable Xc",
          stan::model::index_omni(), stan::model::index_uni((i - 1)));
      }
      current_statement__ = 49;
      stan::model::assign(XQ,
        stan::math::multiply(stan::math::qr_thin_Q(Xc),
          stan::math::sqrt((N - 1))), "assigning variable XQ");
      current_statement__ = 50;
      stan::model::assign(XR,
        stan::math::divide(stan::math::qr_thin_R(Xc),
          stan::math::sqrt((N - 1))), "assigning variable XR");
      current_statement__ = 51;
      stan::model::assign(XR_inv, stan::math::inverse(XR),
        "assigning variable XR_inv");
      current_statement__ = 52;
      stan::math::validate_non_negative_index("bQ", "Kc", Kc);
      current_statement__ = 53;
      stan::math::validate_non_negative_index("z_1", "N_id", N_id);
      current_statement__ = 54;
      stan::math::validate_non_negative_index("b", "Kc", Kc);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = Kc + 1 + 1 + N_id;
  }
  inline std::string model_name() const final {
    return "model_nbhd_glmm";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_nbhd_glmm_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<local_scalar_t__,-1,1> bQ =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(Kc, DUMMY_VAR__);
      current_statement__ = 1;
      bQ = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(Kc);
      local_scalar_t__ Intercept = DUMMY_VAR__;
      current_statement__ = 2;
      Intercept = in__.template read<local_scalar_t__>();
      local_scalar_t__ sd_1 = DUMMY_VAR__;
      current_statement__ = 3;
      sd_1 = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
               lp__);
      Eigen::Matrix<local_scalar_t__,-1,1> z_1 =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N_id, DUMMY_VAR__);
      current_statement__ = 4;
      z_1 = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(N_id);
      {
        current_statement__ = 12;
        if (stan::math::logical_negation(prior_only)) {
          current_statement__ = 7;
          stan::math::validate_non_negative_index("mu", "N", N);
          Eigen::Matrix<local_scalar_t__,-1,1> mu =
            Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N, DUMMY_VAR__);
          current_statement__ = 8;
          stan::model::assign(mu,
            stan::math::add(
              stan::math::add(Intercept, stan::math::multiply(XQ, bQ)),
              stan::math::multiply((20.0 * sd_1),
                stan::model::rvalue(z_1, "z_1", stan::model::index_multi(id)))),
            "assigning variable mu");
          current_statement__ = 9;
          stan::model::assign(mu,
            stan::math::inv_cloglog(stan::model::deep_copy(mu)),
            "assigning variable mu");
          current_statement__ = 10;
          lp_accum__.add(stan::math::bernoulli_lpmf<propto__>(Y, mu));
        }
        current_statement__ = 13;
        lp_accum__.add(stan::math::student_t_lpdf<propto__>(Intercept, 3, 0,
                         2.5));
        current_statement__ = 14;
        lp_accum__.add(stan::math::student_t_lpdf<propto__>(bQ, bQ_prior_df,
                         0, bQ_prior_scale));
        current_statement__ = 15;
        lp_accum__.add(stan::math::gamma_lpdf<propto__>(sd_1, 2.0, 2.0));
        current_statement__ = 16;
        lp_accum__.add(stan::math::std_normal_lpdf<propto__>(z_1));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_nbhd_glmm_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<double,-1,1> bQ =
        Eigen::Matrix<double,-1,1>::Constant(Kc,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      bQ = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(Kc);
      double Intercept = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      Intercept = in__.template read<local_scalar_t__>();
      double sd_1 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 3;
      sd_1 = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
               lp__);
      Eigen::Matrix<double,-1,1> z_1 =
        Eigen::Matrix<double,-1,1>::Constant(N_id,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 4;
      z_1 = in__.template read<Eigen::Matrix<local_scalar_t__,-1,1>>(N_id);
      out__.write(bQ);
      out__.write(Intercept);
      out__.write(sd_1);
      out__.write(z_1);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      Eigen::Matrix<double,-1,1> b =
        Eigen::Matrix<double,-1,1>::Constant(Kc,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 5;
      stan::model::assign(b, stan::math::multiply(XR_inv, bQ),
        "assigning variable b");
      double b_Intercept = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 6;
      b_Intercept = (Intercept - stan::math::dot_product(means_X, b));
      out__.write(b);
      out__.write(b_Intercept);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,1> bQ =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(Kc, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(bQ,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,1>>(Kc),
        "assigning variable bQ");
      out__.write(bQ);
      local_scalar_t__ Intercept = DUMMY_VAR__;
      current_statement__ = 2;
      Intercept = in__.read<local_scalar_t__>();
      out__.write(Intercept);
      local_scalar_t__ sd_1 = DUMMY_VAR__;
      current_statement__ = 3;
      sd_1 = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, sd_1);
      Eigen::Matrix<local_scalar_t__,-1,1> z_1 =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N_id, DUMMY_VAR__);
      current_statement__ = 4;
      stan::model::assign(z_1,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,1>>(N_id),
        "assigning variable z_1");
      out__.write(z_1);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "bQ", "double",
        std::vector<size_t>{static_cast<size_t>(Kc)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "Intercept",
        "double", std::vector<size_t>{});
      current_statement__ = 3;
      context__.validate_dims("parameter initialization", "sd_1", "double",
        std::vector<size_t>{});
      current_statement__ = 4;
      context__.validate_dims("parameter initialization", "z_1", "double",
        std::vector<size_t>{static_cast<size_t>(N_id)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,1> bQ =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(Kc, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> bQ_flat__;
        current_statement__ = 1;
        bQ_flat__ = context__.vals_r("bQ");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= Kc; ++sym1__) {
          current_statement__ = 1;
          stan::model::assign(bQ, bQ_flat__[(pos__ - 1)],
            "assigning variable bQ", stan::model::index_uni(sym1__));
          current_statement__ = 1;
          pos__ = (pos__ + 1);
        }
      }
      out__.write(bQ);
      local_scalar_t__ Intercept = DUMMY_VAR__;
      current_statement__ = 2;
      Intercept = context__.vals_r("Intercept")[(1 - 1)];
      out__.write(Intercept);
      local_scalar_t__ sd_1 = DUMMY_VAR__;
      current_statement__ = 3;
      sd_1 = context__.vals_r("sd_1")[(1 - 1)];
      out__.write_free_lb(0, sd_1);
      Eigen::Matrix<local_scalar_t__,-1,1> z_1 =
        Eigen::Matrix<local_scalar_t__,-1,1>::Constant(N_id, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> z_1_flat__;
        current_statement__ = 4;
        z_1_flat__ = context__.vals_r("z_1");
        current_statement__ = 4;
        pos__ = 1;
        current_statement__ = 4;
        for (int sym1__ = 1; sym1__ <= N_id; ++sym1__) {
          current_statement__ = 4;
          stan::model::assign(z_1, z_1_flat__[(pos__ - 1)],
            "assigning variable z_1", stan::model::index_uni(sym1__));
          current_statement__ = 4;
          pos__ = (pos__ + 1);
        }
      }
      out__.write(z_1);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"bQ", "Intercept", "sd_1", "z_1"};
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"b", "b_Intercept"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(
                                                                    Kc)},
                std::vector<size_t>{}, std::vector<size_t>{},
                std::vector<size_t>{static_cast<size_t>(N_id)}};
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(Kc)},
             std::vector<size_t>{}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= Kc; ++sym1__) {
      param_names__.emplace_back(std::string() + "bQ" + '.' +
        std::to_string(sym1__));
    }
    param_names__.emplace_back(std::string() + "Intercept");
    param_names__.emplace_back(std::string() + "sd_1");
    for (int sym1__ = 1; sym1__ <= N_id; ++sym1__) {
      param_names__.emplace_back(std::string() + "z_1" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= Kc; ++sym1__) {
        param_names__.emplace_back(std::string() + "b" + '.' +
          std::to_string(sym1__));
      }
      param_names__.emplace_back(std::string() + "b_Intercept");
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= Kc; ++sym1__) {
      param_names__.emplace_back(std::string() + "bQ" + '.' +
        std::to_string(sym1__));
    }
    param_names__.emplace_back(std::string() + "Intercept");
    param_names__.emplace_back(std::string() + "sd_1");
    for (int sym1__ = 1; sym1__ <= N_id; ++sym1__) {
      param_names__.emplace_back(std::string() + "z_1" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= Kc; ++sym1__) {
        param_names__.emplace_back(std::string() + "b" + '.' +
          std::to_string(sym1__));
      }
      param_names__.emplace_back(std::string() + "b_Intercept");
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"bQ\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Kc) + "},\"block\":\"parameters\"},{\"name\":\"Intercept\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sd_1\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"z_1\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N_id) + "},\"block\":\"parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Kc) + "},\"block\":\"generated_quantities\"},{\"name\":\"b_Intercept\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"bQ\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Kc) + "},\"block\":\"parameters\"},{\"name\":\"Intercept\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sd_1\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"z_1\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N_id) + "},\"block\":\"parameters\"},{\"name\":\"b\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Kc) + "},\"block\":\"generated_quantities\"},{\"name\":\"b_Intercept\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (((Kc + 1) + 1) + N_id);
    const size_t num_transformed = emit_transformed_parameters * (0);
    const size_t num_gen_quantities = emit_generated_quantities * ((Kc + 1));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (((Kc + 1) + 1) + N_id);
    const size_t num_transformed = emit_transformed_parameters * (0);
    const size_t num_gen_quantities = emit_generated_quantities * ((Kc + 1));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_nbhd_glmm_namespace::model_nbhd_glmm;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_nbhd_glmm_namespace::profiles__;
}
#endif
#endif
