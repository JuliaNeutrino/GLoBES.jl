/* Reference numbers of C GLoBES 3.2.18 for the DUNE configuration (chiMultiExp) */
#include <stdio.h>
#include <math.h>
#include <globes/globes.h>
int main(int argc, char *argv[])
{
  glbInit(argv[0]);
  glbInitExperiment("DUNE_GLoBES.glb", &glb_experiment_list[0], &glb_num_of_exps);
  glb_params t = glbAllocParams(), s = glbAllocParams();
  glbDefineParams(t, asin(sqrt(0.307)), asin(sqrt(0.02203)), asin(sqrt(0.572)), 198.0*M_PI/180, 7.41e-5, 2.511e-3);
  glbSetDensityParams(t, 1.0, GLB_ALL);
  glbSetOscillationParameters(t); glbSetRates();
  int nr = glbGetNumberOfRules(0);
  for (int r = 0; r < nr; r++)
    printf("rule %d signal %.10g bg %.10g\n", r, glbTotalRuleRate(0, r, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG),
           glbTotalRuleRate(0, r, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_BG));
  double dcp[] = {0.0, 90.0, 180.0, 270.0};
  for (int i = 0; i < 4; i++) {
    glbCopyParams(t, s);
    glbSetOscParams(s, dcp[i]*M_PI/180, GLB_DELTA_CP);
    printf("dcp %g chisys_all %.10g", dcp[i], glbChiSys(s, GLB_ALL, GLB_ALL));
    for (int r = 0; r < nr; r++) printf(" r%d %.10g", r, glbChiSys(s, 0, r));
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
    printf(" nosys %.10g\n", glbChiSys(s, GLB_ALL, GLB_ALL));
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON);
  }
  glbCopyParams(t, s);
  glbSetOscParams(s, -2.45e-3, GLB_DM_31);
  printf("IO chisys %.10g\n", glbChiSys(s, GLB_ALL, GLB_ALL));
  return 0;
}
