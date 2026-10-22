#include "engine.c"
/* P(mu->e) and P(mu->mu) of the SNSI + decoherence engine, one layer */
int main(int argc, char *argv[])
{
  glbInit(argv[0]);
  glbRegisterProbabilityEngine(18, &my_probability_matrix, &my_set_oscillation_parameters,
                               &my_get_oscillation_parameters, NULL);
  double th12 = asin(sqrt(atof(argv[1]))), th13 = asin(sqrt(atof(argv[2]))), th23 = asin(sqrt(atof(argv[3])));
  double dcp = atof(argv[4]), sdm = atof(argv[5]), ldm = atof(argv[6]);
  double gamma = atof(argv[7]), eemu = atof(argv[8]), phemu = atof(argv[9]), eetau = atof(argv[10]), phetau = atof(argv[11]);
  glb_params p = glbAllocParams();
  glbDefineParams(p, th12, th13, th23, dcp, sdm, ldm);
  for (int i = 6; i < 18; i++) glbSetOscParams(p, 0.0, i);
  glbSetOscParams(p, gamma, 6); glbSetOscParams(p, gamma, 7); glbSetOscParams(p, gamma, 8);
  glbSetOscParams(p, eemu, 10); glbSetOscParams(p, phemu, 11);
  glbSetOscParams(p, eetau, 12); glbSetOscParams(p, phetau, 13);
  my_set_oscillation_parameters(p, NULL);
  double L = 1300.0, rho = 2.848, P[3][3];
  double Es[] = {0.5, 1.0, 2.5, 5.0, 11.0};
  for (int cp = 1; cp >= -1; cp -= 2)
    for (int k = 0; k < 5; k++) {
      my_probability_matrix(P, cp, Es[k], 1, &L, &rho, 0.0, NULL);
      printf("%d %g %.12e %.12e\n", cp, Es[k], P[1][0], P[1][1]);
    }
  return 0;
}
