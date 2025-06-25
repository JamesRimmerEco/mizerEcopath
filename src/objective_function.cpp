#include <TMB.hpp>

// Helper: single‐ or double‐sigmoid fishing mortality
template<class Type>
vector<Type> calculate_F_mort(
        Type  s1_left,     // ascending intercept
        Type  s2_left,     // ascending slope
        Type  s1_right,    // descending intercept
        Type  s2_right,    // descending slope
        int   use_double,  // 1 = dome, 0 = single
        Type  catchability,
        const vector<Type> &l
) {
    // Ascending limb
    vector<Type> sel_left = Type(1.0) / (Type(1.0) + exp(s1_left - s2_left * l));
    // Descending limb
    vector<Type> sel_right(l.size());
    if (use_double == 1) {
        sel_right = Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * l));
    } else {
        sel_right.setConstant(Type(1.0));
    }
    vector<Type> F_mort = catchability * sel_left.cwiseProduct(sel_right);
    TMBAD_ASSERT((F_mort.array().isFinite() && (F_mort.array() >= 0)).all());
    return F_mort;
}

// Helper function: Steady-state numbers-at-size
template<class Type>
vector<Type> calculate_N(
        const vector<Type> &mort,
        const vector<Type> &growth,
        Type biomass,
        const vector<Type> &w,
        const vector<Type> &dw
) {
    int n = dw.size();
    vector<Type> N(n);
    N(0) = Type(1.0);
    for (int i = 1; i < n; ++i) {
        Type denom = growth(i) + mort(i)*dw(i);
        N(i) = N(i-1) * growth(i-1) / denom;
    }
    TMBAD_ASSERT((N.array().isFinite() && (N.array() >= 0)).all());
    // rescale to match total biomass
    vector<Type> biomass_bins = N * w * dw;
    Type total_biom = biomass_bins.sum();
    N *= (biomass / total_biom);
    return N;
}

// Main objective
template<class Type>
Type objective_function<Type>::operator() () {
    // --- Data ---
    DATA_INTEGER(use_counts);
    DATA_VECTOR(counts);
    DATA_IVECTOR(bin_index);
    DATA_IVECTOR(f_index);
    DATA_VECTOR(coeff_fj);
    DATA_VECTOR(coeff_fj1);

    DATA_INTEGER(use_double); // 1 = double, 0 = single

    DATA_VECTOR(dw);
    DATA_VECTOR(w);
    DATA_VECTOR(l);
    DATA_SCALAR(yield);
    DATA_SCALAR(production);
    DATA_SCALAR(biomass);
    DATA_VECTOR(growth);
    DATA_SCALAR(w_mat);
    DATA_SCALAR(d);
    DATA_SCALAR(yield_lambda);
    DATA_SCALAR(production_lambda);

    // --- Parameters ---
    PARAMETER(l50);
    PARAMETER(ratio);
    PARAMETER(d50);
    PARAMETER(mu_mat);
    PARAMETER(catchability);
    PARAMETER(r_right);

    // Reconstruct the right‐hand 50% point
    Type l50_right = l50 + CppAD::abs(d50);

    // Convert to logistic intercepts/slopes
    Type s1_left  = log(Type(3.0)) * (Type(2.0)*l50) /
        (l50 - l50*ratio);
    Type s2_left  = log(Type(9.0)) /
        (l50 - l50*ratio);
    Type s1_right = log(Type(3.0)) *
        (Type(2.0)*l50_right) /
            (l50_right - l50_right/r_right);
    Type s2_right = log(Type(9.0)) /
        (l50_right - l50_right/r_right);

    // --- Fishing mortality ---
    vector<Type> F_mort = calculate_F_mort(
        s1_left,  s2_left,
        s1_right, s2_right,
        use_double,
        catchability,
        l
    );

    // --- Total mortality & steady‐state N ---
    vector<Type> mort = mu_mat * pow(w / w_mat, d) + F_mort;
    vector<Type> N    = calculate_N(mort, growth, biomass, w, dw);

    // --- Catch density ---
    vector<Type> catch_dens = N * F_mort;

    // --- Negative Log-Likelihood ---
    Type nll = Type(0.0);
    if (use_counts == 1) {
        int nb = counts.size(), ns = bin_index.size();
        vector<Type> probs(nb);
        probs.fill(Type(1e-10));
        for (int k = 0; k < ns; ++k) {
            int i = bin_index(k), j = f_index(k);
            probs(i) += coeff_fj(k)*catch_dens(j)
                + coeff_fj1(k)*catch_dens(j+1);
        }
        probs /= probs.sum();
        nll -= dmultinom(counts, probs, true);
    }

    // Yield penalty
    if (yield_lambda > 0) {
        Type model_yield = (catch_dens * w * dw).sum();
        REPORT(model_yield);
        nll += yield_lambda * pow(log(model_yield/yield), Type(2));
    }
    // Production penalty
    if (production_lambda > 0) {
        Type model_prod = (N * growth * dw).sum();
        REPORT(model_prod);
        nll += production_lambda * pow(log(model_prod/production), Type(2));
    }

    TMBAD_ASSERT(nll >= 0);
    return nll;
}

// Required to register the template
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
