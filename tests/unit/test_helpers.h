// test_helpers.h - small fixtures shared by the unit tests.
#ifndef UNIT_TEST_HELPERS_H_
#define UNIT_TEST_HELPERS_H_

#include "test_access.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace ut {

// Run fn() in a forked child and return its exit status (or 128+signal).
// Used for code paths that call exit(EXIT_FAILURE) ("death tests").
inline int exit_status_of(const std::function<void()>& fn) {
  fflush(stdout); fflush(stderr);
  pid_t pid = fork();
  if (pid == 0) {
    // silence the child's diagnostics
    if (!freopen("/dev/null", "w", stderr)) {}
    if (!freopen("/dev/null", "w", stdout)) {}
    fn();
    _exit(0);
  }
  int st = 0;
  waitpid(pid, &st, 0);
  if (WIFEXITED(st)) return WEXITSTATUS(st);
  if (WIFSIGNALED(st)) return 128 + WTERMSIG(st);
  return -1;
}

// A plain soil layer with explicit properties (no CohortLookup needed).
inline SoilLayer* make_soil_layer(double z, double dz, double poro,
                                  double liq, double ice, double tem,
                                  int frozen, double frozenfrac) {
  SoilLayer* l = new SoilLayer();
  l->isMineral = true;
  l->z = z; l->dz = dz; l->poro = poro;
  l->liq = liq; l->ice = ice; l->tem = tem;
  l->frozen = frozen; l->frozenfrac = frozenfrac;
  l->tcsolid = 2.0; l->bulkden = 2700.0 * (1.0 - poro);
  l->vhcsolid = 2.0e6;
  l->derivePhysicalProperty();   // minliq, maxliq, maxice, tcsat*, tcdry
  l->hksat = 0.005; l->psisat = -100.0; l->bsw = 5.0;
  l->pce_f = 0.; l->pce_t = 0.;
  l->solind = 0; l->indl = 0;
  return l;
}

// Ground with a list of soil layers (top to bottom); Ground owns the layers.
struct Column {
  Ground g;
  std::vector<Layer*> L;
  explicit Column(const std::vector<SoilLayer*>& layers) {
    for (size_t i = 0; i < layers.size(); ++i) {
      g.insertBack(layers[i]);
      L.push_back(layers[i]);
    }
    reindex();
    g.fstsoill = L.front(); g.lstsoill = L.back();
    g.fstminel = L.front(); g.lstminel = L.back();
    g.frontsz.clear(); g.frontstype.clear();
  }
  void reindex() {
    int k = 1;
    for (Layer* l = g.toplayer; l != NULL; l = l->nextl, ++k) {
      l->indl = k; if (l->isSoil) l->solind = k;
    }
  }
  double total_water() const {
    double s = 0.;
    for (size_t i = 0; i < L.size(); ++i) s += L[i]->liq + L[i]->ice;
    return s;
  }
};

} // namespace ut
#endif
