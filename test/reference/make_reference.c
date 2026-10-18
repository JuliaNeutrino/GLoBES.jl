/*
 * Generates test/reference_values.jl from the C GLoBES library.
 *
 * Build and run (from test/data, so that the AEDL files are found):
 *
 *   gcc make_reference.c $(globes-config --include) $(globes-config --libs) -o make_reference
 *   cd ../data && ../reference/make_reference > ../reference_values.jl
 *
 * Every number the Julia test suite compares against comes from here, so the
 * setups below must stay in sync with test/runtests.jl. Indices in the keys
 * are 1-based to match the Julia API.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <globes/globes.h>

static FILE *out;

static void put_vec(const char *key, const double *v, int n)
{
  fprintf(out, "R[\"%s\"] = [", key);
  for (int i = 0; i < n; i++)
    fprintf(out, "%s%.17g", i ? ", " : "", v[i]);
  fprintf(out, "]\n");
}

static void put_num(const char *key, double x)
{
  fprintf(out, "R[\"%s\"] = %.17g\n", key, x);
}

static void put_params(const char *key, glb_params p)
{
  double v[6 + 32];
  int n = 0;
  for (int i = 0; i < 6; i++) v[n++] = glbGetOscParams(p, i);
  for (int i = 0; i < glb_num_of_exps; i++) v[n++] = glbGetDensityParams(p, i);
  put_vec(key, v, n);
}

/* Same as example5.c: flux, far/near fiducial mass and energy calibration */
#define EXP_FAR  0
#define EXP_NEAR 1
static double sq(double x) { return x * x; }
static double chiDCNorm(int exp, int rule, int n_params, double *x, double *errors,
                        void *user_data)
{
  int n_bins = glbGetNumberOfBins(EXP_FAR);
  double *true_N = glbGetRuleRatePtr(EXP_NEAR, 0);
  double *true_F = glbGetRuleRatePtr(EXP_FAR, 0);
  double fit_N[n_bins], fit_F[n_bins];
  double emin, emax, chi2 = 0.0;
  int ew_low, ew_high;
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);
  glbShiftEnergyScale(x[3], glbGetSignalFitRatePtr(EXP_FAR, 0), fit_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_NEAR, 0), fit_N, n_bins, emin, emax);
  double nF = 1.0 + x[0] + x[1], nN = 1.0 + x[0] + x[2];
  for (int i = ew_low; i <= ew_high; i++) {
    if (true_F[i] > 0) chi2 += sq(true_F[i] - nF * fit_F[i]) / true_F[i];
    if (true_N[i] > 0) chi2 += sq(true_N[i] - nN * fit_N[i]) / true_N[i];
  }
  for (int i = 0; i < n_params; i++)
    chi2 += sq(x[i] / errors[i]);
  return chi2;
}

struct var { const char *name; double value; };
struct setup {
  const char *name;
  const char *files[2];
  struct var vars[4];
};

static const struct setup setups[] = {
  { "NOvA",          { "NOvA.glb" } },
  { "NOvA_prepost",  { "NOvA_pre_post_smearing.glb" } },
  { "T2K",           { "T2K.glb" } },
  { "T2HK",          { "T2HK.glb" } },
  { "NFstandard",    { "NFstandard.glb" } },
  { "NFvar",         { "NFvar.glb" },    { { "BASELINE", 3000.0 }, { "emax", 50.0 } } },
  { "NF_hR_lT",      { "NF_hR_lT.glb" }, { { "BASELINE", 4000.0 }, { "emax", 50.0 } } },
  { "BB_350",        { "BB_350.glb" } },
  { "BBvar_WC",      { "BBvar_WC.glb" },
      { { "gammafactor", 350.0 }, { "baselinefactor", 730.0/350.0 }, { "EXP_FACTOR", 1.0 } } },
  { "Reactor1",      { "Reactor1.glb" } },
  { "Reactor2",      { "Reactor2.glb" } },
  { "DChooz",        { "D-Chooz_far.glb", "D-Chooz_near.glb" } },
};
static const int n_setups = sizeof(setups) / sizeof(setups[0]);

/* Built-in chi^2 functions tried on every rule, with their error vectors */
static const char *alt_chi[] = { "chiSpectrumTilt", "chiNoSysSpectrum", "chiSpectrumOnly",
  "chiTotalRatesTilt", "chiNoSysTotalRates", "chiSpectrumCalib", "chiZero" };
static double alt_err[][4] = { { 0.025, 0.01, 0.05, 0.01 }, { 0 }, { 0 },
  { 0.025, 0.01, 0.05, 0.01 }, { 0 }, { 0.025, 0.01, 0.05, 0.01 }, { 0 } };

static void define_params(glb_params p, double th12, double th13, double th23,
                          double delta, double dm21, double dm31, double rho)
{
  glbDefineParams(p, th12, th13, th23, delta, dm21, dm31);
  glbSetDensityParams(p, rho, GLB_ALL);
}

static void true_point(glb_params p)
{ define_params(p, asin(sqrt(0.307)), 0.149, asin(sqrt(0.55)), 1.2*M_PI, 7.4e-5, 2.5e-3, 1.0); }

static void test_point(glb_params p)
{ define_params(p, asin(sqrt(0.307)), 0.149, asin(sqrt(0.45)), 0.5*M_PI, 7.4e-5, 2.45e-3, 1.0); }

static void run_setup(const struct setup *s)
{
  char key[256];
  const char *nm = s->name;

  glbClearAEDLVariables();
  for (int i = 0; i < 4 && s->vars[i].name; i++)
    glbDefineAEDLVariable(s->vars[i].name, s->vars[i].value);
  glbClearExperimentList();
  for (int i = 0; i < 2 && s->files[i]; i++)
    glbInitExperiment((char *)s->files[i], &glb_experiment_list[0], &glb_num_of_exps);

  glb_params p0 = glbAllocParams(), pt = glbAllocParams(), fit = glbAllocParams();
  glb_params err = glbAllocParams();
  true_point(p0);
  test_point(pt);
  define_params(err, 0.1*asin(sqrt(0.307)), 0, 0.1*asin(sqrt(0.55)), 0, 0.05*7.4e-5, 0.05*2.5e-3, 0.05);
  glbSetOscillationParameters(p0);
  glbSetRates();
  glbSetCentralValues(p0);
  glbSetInputErrors(err);
  glbSetStartingValues(p0);

  /* Event rates */
  for (int e = 0; e < glb_num_of_exps; e++) {
    int nb = glbGetNumberOfBins(e);
    snprintf(key, sizeof key, "%s/%d/bincenters", nm, e+1);
    put_vec(key, glbGetBinCentersListPtr(e), nb);
    for (int r = 0; r < glbGetNumberOfRules(e); r++) {
      snprintf(key, sizeof key, "%s/%d/rule%d/rate", nm, e+1, r+1);
      put_vec(key, glbGetRuleRatePtr(e, r), nb);
      snprintf(key, sizeof key, "%s/%d/rule%d/signal", nm, e+1, r+1);
      put_vec(key, glbGetSignalRatePtr(e, r), nb);
      snprintf(key, sizeof key, "%s/%d/rule%d/bg", nm, e+1, r+1);
      put_vec(key, glbGetBGRatePtr(e, r), nb);
      snprintf(key, sizeof key, "%s/%d/rule%d/total", nm, e+1, r+1);
      put_num(key, glbTotalRuleRate(e, r, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG));
    }
    for (int c = 0; c < glbGetNumberOfChannels(e); c++) {
      snprintf(key, sizeof key, "%s/%d/channel%d/post", nm, e+1, c+1);
      put_vec(key, glbGetChannelRatePtr(e, c, GLB_POST), nb);
    }
  }

  /* chi^2 without minimisation over oscillation parameters */
  snprintf(key, sizeof key, "%s/chisys", nm);
  put_num(key, glbChiSys(pt, GLB_ALL, GLB_ALL));
  for (int e = 0; e < glb_num_of_exps; e++)
    for (int r = 0; r < glbGetNumberOfRules(e); r++) {
      snprintf(key, sizeof key, "%s/%d/rule%d/chisys", nm, e+1, r+1);
      put_num(key, glbChiSys(pt, e, r));
    }
  glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
  snprintf(key, sizeof key, "%s/chisys_off", nm);
  put_num(key, glbChiSys(pt, GLB_ALL, GLB_ALL));
  glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);

  /* chi^2 with minimisation */
  snprintf(key, sizeof key, "%s/chitheta13", nm);
  put_num(key, glbChiTheta13(pt, fit, GLB_ALL));
  snprintf(key, sizeof key, "%s/chidelta", nm);
  put_num(key, glbChiDelta(pt, fit, GLB_ALL));
  snprintf(key, sizeof key, "%s/chidelta_fit", nm);
  put_params(key, fit);
  snprintf(key, sizeof key, "%s/chitheta13delta", nm);
  put_num(key, glbChiTheta13Delta(pt, NULL, GLB_ALL));

  glb_projection proj = glbAllocProjection();
  glbDefineProjection(proj, GLB_FIXED, GLB_FIXED, GLB_FREE, GLB_FIXED, GLB_FIXED, GLB_FREE);
  glbSetDensityProjectionFlag(proj, GLB_FREE, GLB_ALL);
  glbSetProjection(proj);
  snprintf(key, sizeof key, "%s/chinp", nm);
  put_num(key, glbChiNP(pt, NULL, GLB_ALL));
  glbFreeProjection(proj);

  /* Built-in chi^2 functions (DChooz uses its own, defined in the AEDL file) */
  if (strcmp(nm, "DChooz") != 0) {
    for (int k = 0; k < (int)(sizeof alt_chi / sizeof alt_chi[0]); k++) {
      glbSetChiFunction(GLB_ALL, GLB_ALL, GLB_ON, alt_chi[k], alt_err[k]);
      snprintf(key, sizeof key, "%s/alt/%s", nm, alt_chi[k]);
      put_num(key, glbChiSys(pt, GLB_ALL, GLB_ALL));
    }

    /* The same functions at fixed nuisance parameters, without minimisation.
     * The fit rates are still those of the test point from the last glbChiSys. */
    static double x[4] = { 0.01, -0.005, 0.02, 0.003 };
    static const int dims[] = { 4, 0, 1, 4, 0, 4, 0 };
    for (int e = 0; e < glb_num_of_exps; e++)
      for (int r = 0; r < glbGetNumberOfRules(e); r++)
        for (int k = 0; k < (int)(sizeof alt_chi / sizeof alt_chi[0]); k++) {
          snprintf(key, sizeof key, "%s/%d/rule%d/fixed/%s", nm, e+1, r+1, alt_chi[k]);
          put_num(key, glbGetChiFunctionPtr(alt_chi[k])(e, r, dims[k], x, alt_err[k], NULL));
        }
  }

  /* Scaling the exposure */
  glbSetTargetMass(0, 2.0 * glbGetTargetMass(0));
  glbSetRunningTime(0, 0, 0.5 * glbGetRunningTime(0, 0));
  glbSetRates();
  for (int r = 0; r < glbGetNumberOfRules(0); r++) {
    snprintf(key, sizeof key, "%s/1/rule%d/total_scaled", nm, r+1);
    put_num(key, glbTotalRuleRate(0, r, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG));
  }

  glbFreeParams(p0); glbFreeParams(pt); glbFreeParams(fit); glbFreeParams(err);
}

static void probabilities(void)
{
  static const double energies[] = { 0.3, 0.5, 1.0, 2.0, 3.5, 7.0, 15.0, 40.0 };
  const int nE = sizeof energies / sizeof energies[0];
  double v[9 * 2 * nE];
  char key[64];
  glb_params p0 = glbAllocParams();
  true_point(p0);

  glbClearAEDLVariables();
  glbClearExperimentList();
  glbInitExperiment("NFstandard.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glbSetOscillationParameters(p0);

  for (int kind = 0; kind < 3; kind++) {
    int n = 0;
    for (int s = 0; s < 2; s++)
      for (int k = 0; k < nE; k++)
        for (int a = 1; a <= 3; a++)
          for (int b = 1; b <= 3; b++) {
            int cp = s ? -1 : +1;
            if (kind == 0)      v[n++] = glbVacuumProbability(a, b, cp, energies[k], 1300.0);
            else if (kind == 1) v[n++] = glbConstantDensityProbability(a, b, cp, energies[k], 1300.0, 2.8);
            else                v[n++] = glbProfileProbability(0, a, b, cp, energies[k]);
          }
    snprintf(key, sizeof key, "prob/%s", kind == 0 ? "vacuum" : kind == 1 ? "constant" : "profile");
    put_vec(key, v, n);
  }
  put_vec("prob/energies", energies, nE);

  /* Energy scale shift of a fixed spectrum */
  double in[20], sh[20];
  for (int i = 0; i < 20; i++) in[i] = 100.0 * exp(-0.5 * sq((i - 8.0) / 3.0)) + i;
  static const double shifts[] = { -0.3, -0.04, 0.01, 0.07, 0.3 };
  for (int k = 0; k < (int)(sizeof shifts / sizeof shifts[0]); k++) {
    glbShiftEnergyScale(shifts[k], in, sh, 20, 0.5, 10.5);
    snprintf(key, sizeof key, "shift/%d", k+1);
    put_vec(key, sh, 20);
  }
  glbFreeParams(p0);
}

int main(int argc, char *argv[])
{
  out = stdout;
  glbInit(argv[0]);
  glbSetVerbosityLevel(0);
  glbDefineChiFunction(&chiDCNorm, 5, "chiDCNorm", NULL);

  fprintf(out, "# Generated by test/reference/make_reference.c from C GLoBES -- do not edit.\n");
  fprintf(out, "const R = Dict{String,Any}()\n");
  probabilities();
  for (int i = 0; i < n_setups; i++)
    run_setup(&setups[i]);
  return 0;
}
