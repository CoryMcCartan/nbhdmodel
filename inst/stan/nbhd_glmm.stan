functions {
    // compute the cloglog link [0, 1] -> (-inf, inf)
    vector cloglog(vector p) {
        return log(-log(1 - p));
    }
}

data {
    int<lower=1> N;  // total number of observations
    int Y[N];  // response variable
    int<lower=1> K;  // number of population-level effects
    matrix[N, K] X;  // population-level design matrix
    // data for group-level effects of ID 1
    int<lower=1> N_id;  // number of grouping levels
    int<lower=1> id[N];  // grouping indicator per observation
    real bQ_prior_scale; // scale of QR prior
    int<lower=1> bQ_prior_df; // df of QR prior
    int prior_only;  // should the likelihood be ignored?
}

transformed data {
    int Kc = K - 1;
    matrix[N, Kc] Xc;  // centered version of X without an intercept
    vector[Kc] means_X;  // column means of X before centering
    // matrices for QR decomposition
    matrix[N, Kc] XQ;
    matrix[Kc, Kc] XR;
    matrix[Kc, Kc] XR_inv;
    for (i in 2:K) {
        means_X[i - 1] = mean(X[, i]);
        Xc[, i - 1] = X[, i] - means_X[i - 1];
    }
    // compute and scale QR decomposition
    XQ = qr_thin_Q(Xc) * sqrt(N - 1);
    XR = qr_thin_R(Xc) / sqrt(N - 1);
    XR_inv = inverse(XR);
}

parameters {
    vector[Kc] bQ;  // regression coefficients at QR scale
    real Intercept;  // temporary intercept for centered predictors
    real<lower=0> sd_1;  // group-level standard deviations
    vector[N_id] z_1;  // standardized group-level effects
}

model {
    if (!prior_only) {
        vector[N] mu = Intercept + XQ * bQ + 20.0*sd_1*z_1[id];
        mu = inv_cloglog(mu);
        Y ~ bernoulli(mu);
    }

    Intercept ~ student_t(3, 0, 2.5);
    bQ ~ student_t(bQ_prior_df, 0, bQ_prior_scale);
    //sd_1 ~ student_t(3, 0, 2.5);
    sd_1 ~ gamma(2.0, 2.0);
    z_1 ~ std_normal();
}

generated quantities {
    // obtain the actual coefficients
    vector[Kc] b = XR_inv * bQ;
    // actual population-level intercept
    real b_Intercept = Intercept - dot_product(means_X, b);
}
