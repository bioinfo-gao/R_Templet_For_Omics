Using brms
Peter Ralph
Wed Nov 9 20:59:03 2022
Running brms
The paper introducing brms, available through vignette("brms_overview"), provides a good technical introduction to how brms works.

Some other useful places for documentation are:
    
    Vignettes listed under vignette(package='brms')
help(brmsformula) - very detailed! Describes all (?) the modeling options.
help(brmsfamily) - lists defaults, etcetera for each family.
help(set_prior) for how to choose priors
methods(class = "brmsfit") for methods of a brmsfit object
Here, we’ll look under the hood, in detail, at one of its examples, the kidney dataset, whose analysis is described in the brms_overview vignette. This has data describing the first and (possibly) second recurrence time of infection in thirty-eight patients.

library(brms)
data(kidney)
kidney$patient <- factor(kidney$patient)
head(kidney)
##   time censored patient recur age    sex disease
## 1    8        0       1     1  28   male   other
## 2   23        0       2     1  48 female      GN
## 3   22        0       3     1  32   male   other
## 4  447        0       4     1  31 female   other
## 5   30        0       5     1  10   male   other
## 6   24        0       6     1  16 female   other
Here is the model described in the paper:
    
    fit1 <- brm(formula = time | cens(censored) ~ age * sex + disease + (1 + age|patient),
                data = kidney, family = lognormal(),
                prior = c(set_prior("normal(0,5)", class = "b"),
                          set_prior("cauchy(0,2)", class = "sd"),
                          set_prior("lkj(2)", class = "cor")),
                warmup = 1000, iter = 2000, chains = 4,
                control = list(adapt_delta = 0.95))
## Compiling Stan program...
## Start sampling
fit1
##  Family: lognormal 
##   Links: mu = identity; sigma = identity 
## Formula: time | cens(censored) ~ age * sex + disease + (1 + age | patient) 
##    Data: kidney (Number of observations: 76) 
##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
##          total post-warmup draws = 4000
## 
## Group-Level Effects: 
## ~patient (Number of levels: 38) 
##                    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sd(Intercept)          0.41      0.29     0.02     1.09 1.00     1179     2016
## sd(age)                0.01      0.01     0.00     0.02 1.00      794     1249
## cor(Intercept,age)    -0.17      0.46    -0.89     0.74 1.00     1985     1817
## 
## Population-Level Effects: 
##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## Intercept         2.74      0.98     0.84     4.70 1.00     1475     1999
## age               0.01      0.03    -0.03     0.06 1.00     1070     1629
## sexfemale         2.42      1.16     0.13     4.74 1.00     1229     1867
## diseaseGN        -0.37      0.53    -1.42     0.66 1.00     1930     2045
## diseaseAN        -0.50      0.51    -1.54     0.51 1.00     2360     2636
## diseasePKD        0.62      0.73    -0.78     2.01 1.00     2053     2693
## age:sexfemale    -0.02      0.03    -0.08     0.03 1.00     1031     1628
## 
## Family Specific Parameters: 
##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## sigma     1.14      0.13     0.91     1.43 1.00     2271     2467
## 
## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).
Note that the cens(censored) call is documented in help("addition-terms") and help(brmsformula), which says that

With the exception of categorical, ordinal, and mixture families,
left, right, and interval censoring can be modeled through ‘y |
    cens(censored) ~ predictors’. The censoring variable (named
                                                          ‘censored’ in this example) should contain the values ‘'left'’,
‘'none'’, ‘'right'’, and ‘'interval'’ (or equivalently ‘-1’, ‘0’,
                                       ‘1’, and ‘2’) to indicate that the corresponding observation is
left censored, not censored, right censored, or interval censored.
For interval censored data, a second variable (let's call it ‘y2’)
     has to be passed to ‘cens’. In this case, the formula has the
     structure ‘y | cens(censored, y2) ~ predictors’.  While the lower
     bounds are given in ‘y’, the upper bounds are given in ‘y2’ for
     interval censored data. Intervals are assumed to be open on the
     left and closed on the right: ‘(y, y2]’.
The model
Translating the formula, we’re trying to fit the following Bayesian hierarchical model. Let timei
 be the i
th time, with censoredi=0
 if this is a recurrence, and censoredi=1
 otherwise (i.e., if it is censored). Let’s use Si
 to denote the actual recurrence time. Then the basic model is
timeitimeiSi=Siif censoredi=0<Siif censoredi=1∼logNormal(μi,σ).
How does μ
 relate to the linear predictor? Well, since a link is not specified, checking help(brmsfamily) tells us that the defaults are

     lognormal(link = "identity", link_sigma = "log")
and so the parameterization will be equivalent to:
μi=α+βage, sexiagei+βdiseasediseasei+νpatienti+νage, patientiagei
Here α
 is the intercept, and the various β
 and ν
 are other parameters. This is overparameterized; we’ll verify what parameters are actually used below.

For priors, we have for the “fixed effects”:
βage, s∼Normal(0,5).
There is no prior specified in the call for the (“random” or “group”) patient-specific effects, ν
, and in fact there doesn’t seem to be a way to specify their priors directly: the prior on (νp,νa,p)
 will always be Normal, (but by default correlated between patients p
). The prior standard deviations of each of these two effects are specified, though; this is the class = "sd" argument to brm( ) above:
(νpatienti,νpatienti, agei)σ1σ2L∼Normal(0,diag(σ)LLTdiag(σ))∼Cauchy(0,2)∼Cauchy(0,2)∼LKJ(2)
The last term is an LKJ prior on the correlation matrix.

The only remaining aspect of the model (not specified in the brm() call) is the intercept, α
. For discussion of the prior on this (which is somewhat involved because of mean-centering), see help(brmsformula).

The Stan model
Let’s make sure we understand exactly what’s going on here. Use the source, Rey: let’s look at the underlying Stan code. The call stancode(fit1) produces this, which we’ll look at in pieces.

The data block
First, comes the data block:

// generated with brms 2.10.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=-1,upper=2> cens[N];  // indicates censoring
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
It is pretty generic, since presumably the same structure gets used for many types of model, but it is commented – nice! We can check that the data that brms actually passes to Stan matches this: standata(fit1) returns this:

str(standata(fit1))
## List of 12
##  $ N         : int 76
##  $ Y         : num [1:76(1d)] 8 23 22 447 30 24 7 511 53 15 ...
##  $ cens      : num [1:76(1d)] 0 0 0 0 0 0 0 0 0 0 ...
##  $ K         : int 7
##  $ X         : num [1:76, 1:7] 1 1 1 1 1 1 1 1 1 1 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:76] "1" "2" "3" "4" ...
##   .. ..$ : chr [1:7] "Intercept" "age" "sexfemale" "diseaseGN" ...
##   ..- attr(*, "assign")= int [1:7] 0 1 2 3 3 3 4
##   ..- attr(*, "contrasts")=List of 2
##   .. ..$ sex    : num [1:2, 1] 0 1
##   .. .. ..- attr(*, "dimnames")=List of 2
##   .. .. .. ..$ : chr [1:2] "male" "female"
##   .. .. .. ..$ : chr "female"
##   .. ..$ disease: num [1:4, 1:3] 0 1 0 0 0 0 1 0 0 0 ...
##   .. .. ..- attr(*, "dimnames")=List of 2
##   .. .. .. ..$ : chr [1:4] "other" "GN" "AN" "PKD"
##   .. .. .. ..$ : chr [1:3] "GN" "AN" "PKD"
##  $ Z_1_1     : num [1:76(1d)] 1 1 1 1 1 1 1 1 1 1 ...
##   ..- attr(*, "dimnames")=List of 1
##   .. ..$ : chr [1:76] "1" "2" "3" "4" ...
##  $ Z_1_2     : num [1:76(1d)] 28 48 32 31 10 16 51 55 69 51 ...
##   ..- attr(*, "dimnames")=List of 1
##   .. ..$ : chr [1:76] "1" "2" "3" "4" ...
##  $ J_1       : int [1:76(1d)] 1 2 3 4 5 6 7 8 9 10 ...
##  $ N_1       : int 38
##  $ M_1       : int 2
##  $ NC_1      : int 1
##  $ prior_only: int 0
##  - attr(*, "class")= chr [1:2] "standata" "list"
A key object here is X, the design matrix for the covariates. In this case, this is

head(standata(fit1)$X)
##   Intercept age sexfemale diseaseGN diseaseAN diseasePKD age:sexfemale
## 1         1  28         0         0         0          0             0
## 2         1  48         1         1         0          0            48
## 3         1  32         0         0         0          0             0
## 4         1  31         1         0         0          0            31
## 5         1  10         0         0         0          0             0
## 6         1  16         1         0         0          0            16
The transformed data block just centers the columns of X, producing Xc:

transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
Note also that the term (1 + age | patient) has produced standata(fit1)$M_1 = 2 coefficients for each of the standata(fit1)$N_1 = 38 levels of the patient factor.

The parameters block
There are five (vectors or matrices) of coefficents and one (Cholesky factor of a) correlation matrix in the parameters block:

parameters {
  vector[Kc] b;  // population-level effects
  // temporary intercept for centered predictors
  real Intercept;
  real<lower=0> sigma;  // residual SD
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_1] L_1;
}
Here, there is:

one b for each column of the design matrix, X;
an Intercept;
a standard deviation, sigma, for the logNormal response
two sd_1s for the patient effects: one for the patient-specific intercept, and one for the patient-specific age effect
a (2 x 38) matrix z_1 encoding each of the patient effects, decorrelated
a (Cholesky factor for the 2 x 2) correlation matrix for the two patient effects
In transformed parameters we get

r_1, the (38 x 2) matrix of the actual patent effects
and some stuff for optimization.

transformed parameters {
  // actual group-level effects
  matrix[N_1, M_1] r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
// using vectors speeds up indexing in loops
vector[N_1] r_1_1 = r_1[, 1];
vector[N_1] r_1_2 = r_1[, 2];
}
Model block
First, we construct the linear predictor, mu:
    
    model {
        // initialize linear predictor term
        vector[N] mu = Intercept + Xc * b;
        for (n in 1:N) {
            // add more terms to the linear predictor
            mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
        }
        Then, comes the sampling statements. Recall that target += normal_lpdf(b | 0, 5); is the same as b ~ normal(0, 5);, so these are priors:
            
            // priors including all constants
        target += normal_lpdf(b | 0,5);
        target += student_t_lpdf(Intercept | 3, 4, 10);
        The strange terms here, e.g., subtracting off student_t_lccdf(0 | ...) are because we’ve got a half Student’s t
        prior (i.e., we condition on it being above 0), so brms is (admirably) keeping the log posterior density correct, even up to constants. (Here, _lcdf is “log cumulative distribution function”, and _lccdf is “log complementary cumulative distribution function.)
        
        target += student_t_lpdf(sigma | 3, 0, 10)
        - 1 * student_t_lccdf(0 | 3, 0, 10);
        target += cauchy_lpdf(sd_1 | 0,2)
        - 2 * cauchy_lccdf(0 | 0,2);
        target += normal_lpdf(to_vector(z_1) | 0, 1);
        target += lkj_corr_cholesky_lpdf(L_1 | 2);
        And, here’s where censoring comes in: if the observation is not censored (if cens[n] == 0) then the sampling statement is as usual, adding lognormal_lpdf(Y[n] | mu[n], sigma); but if it is right-censored (cens[n] == 1), then we need to add the log probability that the logNormal is greater than Y[n], and so add lognormmal_lccdf(Y[n] | mu[n], sigma).
        
        // likelihood including all constants
        if (!prior_only) {
            for (n in 1:N) {
                // special treatment of censored data
                if (cens[n] == 0) {
                    target += lognormal_lpdf(Y[n] | mu[n], sigma);
                } else if (cens[n] == 1) {
                    target += lognormal_lccdf(Y[n] | mu[n], sigma);
                } else if (cens[n] == -1) {
                    target += lognormal_lcdf(Y[n] | mu[n], sigma);
                }
            }
        }
    }
And, that’s it for the model!
    
    Generated quantities
Here’s where two normalizations for optimization get undone: the mean-centering of X, and the decorrelation of patient effects, extracting the non-mean-centered intercept and the actual correlation matrix of patient effects.

generated quantities {
    // actual population-level intercept
    real b_Intercept = Intercept - dot_product(means_X, b);
    // group-level correlations
    corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
    vector<lower=-1,upper=1>[NC_1] cor_1;
    // extract upper diagonal of correlation matrix
    for (k in 1:M_1) {
        for (j in 1:(k - 1)) {
            cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
        }
    }
}
Extracting information
Ok, now that we know precisely what’s being estimated, what information can we get out of the fitted model?
    
    One thing useful to know about is the concept of distributional parameters, or “dpars”. These are the parameters to the response distribution: for instance, those of the logNormal are mu and sigma:
    
    fit1$family$dpars
## [1] "mu"    "sigma"
Since the (transformed) linear predictor is turned into mu, this means that each observation has it’s own posterior distribution on mu (and on sigma, too).

posterior_samples( ) : All the samples
A straightforward way to get the samples out of the brmsfit object is to use posterior_samples( ) (or as.matrix( )) on it; this will give you one column per parameter and one row per sample. Specific parameters can be obtained by the pars argument (which takes regular expressions). For instance, here are posterior samples of the “slope” parameters:
    
    head(posterior_samples(fit1, pars="b_"))
##   b_Intercept        b_age b_sexfemale b_diseaseGN b_diseaseAN b_diseasePKD b_age:sexfemale
## 1    2.041359  0.039193905   4.4811974  -0.8015959 -0.65459314    0.3064934     -0.07510717
## 2    4.927252 -0.032136290  -0.2090837  -0.4479103 -0.09809829    0.8209213      0.02781259
## 3    4.323357 -0.027484442   0.5044778  -0.4336619  0.10034811    0.6383120      0.02536837
## 4    4.169794 -0.022590576   0.1477550  -0.4698482 -0.86476422    0.9181340      0.03824760
## 5    5.125690 -0.042474171  -0.6314809  -0.5944212 -0.39030076    1.2497668      0.04456000
## 6    2.800304  0.009680656   2.1796588  -0.2117511 -0.75716502    0.2059602     -0.01259277
and here are the SD and correlation parameters:
    
    head(posterior_samples(fit1, pars=c("sd", "cor", "sigma")))
##   sd_patient__Intercept sd_patient__age cor_patient__Intercept__age    sigma
## 1           0.243324980    1.431458e-04                  -0.6145436 1.293972
## 2           0.009494062    1.995648e-05                   0.5904942 1.139450
## 3           0.025387962    6.603650e-05                   0.8473875 1.308999
## 4           0.268901508    3.632726e-05                   0.9405463 1.302209
## 5           0.103122401    1.234581e-05                   0.7234974 1.345990
## 6           0.044433133    1.663014e-05                  -0.5477178 1.022988
and here are some of the patient-specific effects:
    
    head(posterior_samples(fit1, pars="^r_"))[,1:5]
##   r_patient[1,Intercept] r_patient[2,Intercept] r_patient[3,Intercept] r_patient[4,Intercept] r_patient[5,Intercept]
## 1             0.50150450           -0.086635015            -0.38823154           -0.146653443           -0.161207530
## 2            -0.01450695            0.001081755             0.01151147            0.013138762           -0.002226891
## 3            -0.03131236           -0.026826527            -0.00705684            0.006915832           -0.046216653
## 4            -0.11871762           -0.595192980             0.09752397            0.293331071           -0.229389286
## 5            -0.15855766           -0.172572301             0.08654596           -0.012225816           -0.136896996
## 6            -0.02309118            0.032195989            -0.03540636            0.030838992           -0.033534607
fixef( ) : fixed effects
The “fixed effects” (i.e., the b_ parameters) can be extracted with fixef( ), either as a summary table (if summary=TRUE) or as a matrix of samples (otherwise). Here it is for this model:
    
    fixef(fit1)
##                  Estimate  Est.Error        Q2.5      Q97.5
## Intercept      2.74297839 0.98179304  0.84025703 4.69512966
## age            0.01402661 0.02530538 -0.03450567 0.06385842
## sexfemale      2.41837062 1.16469139  0.12983696 4.73713428
## diseaseGN     -0.36751978 0.52518004 -1.41939034 0.66052565
## diseaseAN     -0.50043002 0.50891622 -1.53753157 0.50840965
## diseasePKD     0.62207929 0.72807200 -0.78246718 2.00589764
## age:sexfemale -0.02159992 0.02681184 -0.07532626 0.03058815
The Estimate and Est.Error columns give the mean and SD, respectively, of the posterior distributions of the listed parameters.

ranef( ) : “random” effects
The “group-level” effects - here, the patient-specific intercepts and slopes - are returned, as a 3D array, by ranef( ), in a similar format to fixef( ). For instance, here are the patient-specific age effects:
    
    head(ranef(fit1, pars='age'))
## $patient
## , , age
## 
##         Estimate   Est.Error        Q2.5      Q97.5
## 1  -1.228074e-03 0.010599731 -0.02394011 0.01987465
## 2  -2.279677e-03 0.010314875 -0.02698522 0.01726670
## 3  -5.515474e-06 0.010235396 -0.02128074 0.02219155
## 4   2.070935e-03 0.010917735 -0.01841183 0.02733859
## 5  -6.172472e-05 0.010461073 -0.02199139 0.02238146
## 6  -6.709983e-04 0.010372109 -0.02295300 0.02056866
## 7  -3.703092e-03 0.010176718 -0.02817039 0.01383671
## 8   1.320567e-03 0.009127154 -0.01692709 0.02216013
## 9   1.837007e-03 0.008716498 -0.01465174 0.02185969
## 10  2.953904e-03 0.009477916 -0.01458713 0.02587324
## 11 -1.739521e-03 0.009600864 -0.02522923 0.01581569
## 12  2.007554e-04 0.010316832 -0.02179273 0.02228575
## 13 -6.964507e-04 0.010144372 -0.02347726 0.02073956
## 14  3.989352e-03 0.011400428 -0.01563742 0.03179674
## 15  5.205242e-04 0.011170096 -0.02387479 0.02428181
## 16  1.859519e-04 0.009563485 -0.01937539 0.02218193
## 17  1.799486e-03 0.009351846 -0.01543432 0.02339253
## 18  1.141962e-03 0.009790292 -0.01789273 0.02292219
## 19  3.862817e-03 0.010877725 -0.01392571 0.03059823
## 20 -1.732127e-03 0.010033630 -0.02532850 0.01719203
## 21  5.398422e-03 0.011829486 -0.01389832 0.03390025
## 22  1.512313e-03 0.011292346 -0.02129536 0.02654629
## 23 -3.082049e-03 0.008824515 -0.02366211 0.01260636
## 24  2.506236e-04 0.010066061 -0.02010505 0.02286438
## 25  5.556263e-04 0.009810147 -0.01948684 0.02167248
## 26  4.950458e-03 0.010561234 -0.01166078 0.03069769
## 27 -1.562869e-04 0.010528780 -0.02442324 0.02182849
## 28 -2.905325e-03 0.009618501 -0.02562252 0.01363890
## 29 -4.347550e-03 0.010633456 -0.03013918 0.01343027
## 30 -1.673425e-03 0.009216802 -0.02476945 0.01513619
## 31 -2.265434e-03 0.009125489 -0.02299769 0.01405937
## 32 -9.540620e-04 0.009491849 -0.02314714 0.01714449
## 33 -4.592717e-03 0.009872965 -0.03007115 0.01085338
## 34  1.612935e-03 0.009846022 -0.01776378 0.02412524
## 35 -1.761414e-03 0.011630702 -0.02765050 0.02292380
## 36  1.276532e-03 0.010091617 -0.01966985 0.02498874
## 37 -2.320132e-03 0.010660002 -0.02752366 0.01645066
## 38  1.779828e-04 0.009750630 -0.02080964 0.02130386
predict( ) : the responses
The predict.brmsfit( ) method will produce the posterior distribution of the responses, which is in this case the times. This is either the recurrence time, for uncensored observations, or a uniform time between 0 and the recurrence time, for censored ones. With summary=FALSE, this gets a (num steps) x (num data points matrix of samples from the posterior distribution, and with summary=TRUE it returns a summary table, e.g.:
                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                     head(predict(fit1))
                                                                                                                                                                                                                                                                                                                 ##       Estimate Est.Error      Q2.5     Q97.5
                                                                                                                                                                                                                                                                                                                 ## [1,]  44.50189  97.24697  1.655912  239.8330
                                                                                                                                                                                                                                                                                                                 ## [2,] 156.46379 359.46260  5.234299  806.1940
                                                                                                                                                                                                                                                                                                                 ## [3,]  55.43995 105.98715  2.233496  303.8072
                                                                                                                                                                                                                                                                                                                 ## [4,] 384.20582 699.12042 14.321316 2215.3332
                                                                                                                                                                                                                                                                                                                 ## [5,]  45.25001 120.27611  1.346303  267.4825
                                                                                                                                                                                                                                                                                                                 ## [6,] 288.55267 544.61076 11.654704 1525.1934
                                                                                                                                                                                                                                                                                                                 Under the hood, this just calls rlnorm( ) with posterior samples from the dpars,
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 fitted( ) : the “mean” response
                                                                                                                                                                                                                                                                                                                 The method fitted.brmsfit( ) can produce the posterior distribution of either the linear predictor (if scale="linear") or the mean response (if scale="response"). Again, summary can be used to return either a summary table or the full set of samples.
                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                 head(fitted(fit1, summary=TRUE))
                                                                                                                                                                                                                                                                                                                 ##       Estimate Est.Error       Q2.5    Q97.5
                                                                                                                                                                                                                                                                                                                 ## [1,]  44.15714  26.78380  12.235186 111.6373
                                                                                                                                                                                                                                                                                                                 ## [2,] 160.10570 114.80782  36.559616 453.5086
                                                                                                                                                                                                                                                                                                                 ## [3,]  55.18891  33.55914  17.170841 141.8991
                                                                                                                                                                                                                                                                                                                 ## [4,] 392.70267 231.08303 139.904747 971.5443
                                                                                                                                                                                                                                                                                                                 ## [5,]  45.91564  41.22531   9.151228 154.0889
                                                                                                                                                                                                                                                                                                                 ## [6,] 297.94963 177.68905  89.744473 748.5797
                                                                                                                                                                                                                                                                                                                 You can also pass fitted( ) the name of a distributional parameter to have it return the posterior distribution of those parameters (one per data point). For instance, here are summaries of the posterior of mu the top few data points:
                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                     head(fitted(fit1, dpar='mu'))
                                                                                                                                                                                                                                                                                                                 ##      Estimate Est.Error     Q2.5    Q97.5
                                                                                                                                                                                                                                                                                                                 ## [1,] 2.969156 0.5252858 1.889021 3.958614
                                                                                                                                                                                                                                                                                                                 ## [2,] 4.221122 0.5823402 2.995580 5.352074
                                                                                                                                                                                                                                                                                                                 ## [3,] 3.204536 0.5058786 2.216919 4.237261
                                                                                                                                                                                                                                                                                                                 ## [4,] 5.180230 0.4991074 4.258862 6.233864
                                                                                                                                                                                                                                                                                                                 ## [5,] 2.895303 0.7027871 1.521475 4.315845
                                                                                                                                                                                                                                                                                                                 ## [6,] 4.887733 0.5035245 3.902700 5.848140
                                                                                                                                                                                                                                                                                                                 Under the hood, the standard fitted( ) function uses the samples of the dpars and the analytical relationship between those and the mean. For instance, in this case, here’s what ends up being called:
                                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                     > brms:::fitted_lognormal
                                                                                                                                                                                                                                                                                                                 function (draws)
                                                                                                                                                                                                                                                                                                                 {
                                                                                                                                                                                                                                                                                                                     with(draws$dpars, exp(mu + sigma^2/2))
                                                                                                                                                                                                                                                                                                                 }