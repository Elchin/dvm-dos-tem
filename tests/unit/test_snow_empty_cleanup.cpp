// model_settings.snow_empty_cleanup (fix/snow-empty-layers): root cause of the r1 1 Jun 2024 melt-out closure miss
// (fix/meltout-closure, design doc 15g). When the thermal step melts a lower snow layer completely while the layer
// above keeps ice, the daily snow update left the emptied layer in the column: constructSnowLayers() drops empty
// layers from the top only, and combineSnowLayers() tests the pre-compaction dz (dz -> 0 only in updateThick()).
// The layer (ice = liq = 0, dz = 0, zero heat capacity) then entered the next day's thermal step.
// snow_empty_cleanup = 1 removes such layers at the end of Snow_Env::updateDailyM (both solvers); 0 = legacy.
#include "doctest.h"
#include "test_helpers.h"
#include <cmath>

using ut::make_soil_layer;

namespace {
SnowLayer* snl(double dz, double ice, double liq, double tem, double rho, int frozen = 1) {
  SnowLayer* s = new SnowLayer(); s->dz = dz; s->ice = ice; s->liq = liq; s->tem = tem; s->rho = rho; s->z = 0.0;
  s->frozen = frozen; s->age = 10.0; return s;
}
std::vector<SoilLayer*> soil3() {
  return {make_soil_layer(0.0, 0.05, 0.6, 15, 0, -2.0, 1, 1), make_soil_layer(0.05, 0.1, 0.5, 20, 0, -3.0, 1, 1),
          make_soil_layer(0.15, 0.2, 0.4, 30, 0, -4.0, 1, 1)};
}
struct SnowCol {
  ut::Column col; EnvData ed; CohortData cd; Snow_Env sn; Soil_Env se;
  explicit SnowCol(const std::vector<SnowLayer*>& snow, int cleanup = 1, int solver = 0) : col(soil3()) {
    for (size_t k = snow.size(); k-- > 0;) { col.g.insertFront(snow[k]); }
    col.reindex(); col.g.fstsoill = col.L.front();
    col.g.snowdimpar.newden = 50.0; col.g.snowdimpar.denmax = 250.0;
    col.g.snow_empty_cleanup = cleanup;
    col.g.resortGroundLayers();
    sn.setGround(&col.g); sn.setEnvData(&ed); sn.setCohortData(&cd);
    sn.snowenvpar.albmax = 0.8; sn.snowenvpar.albmin = 0.4; sn.thermal_solver = solver;
    se.setGround(&col.g); se.setEnvData(&ed); se.setCohortData(&cd); se.snowenvp = &sn;
    se.thermal_solver = solver; se.thermal_sfcc = 1; se.thermal_substeps = 4;
    cd.m_vegd.fpcsum = 0.0; ed.d_v2g.sthfl = 0.0; ed.d_v2g.sdrip = 0.0; ed.d_a2l.snfl = 0.0;
    ed.d_v2g.swthfl = 0.0; ed.d_a2l.nirr = 0.0; ed.d_soid.nfactor = 1.0; ed.d_snw2a.sublim = 0.0; ed.d_snw2soi.melt = 0.0;
  }
  int nsnow() const { int n = 0; for (Layer* l = col.g.toplayer; l && l->isSnow; l = l->nextl) n++; return n; }
  double snowW() const { double w = 0; for (Layer* l = col.g.toplayer; l && l->isSnow; l = l->nextl) w += l->ice + l->liq; return w; }
  double allW() const { double w = 0; for (Layer* l = col.g.toplayer; l; l = l->nextl) w += l->ice + l->liq; return w; }
};
// list integrity: forward/backward walks agree, prev/next consistent, indices 1..n, snow above soil, fstsoill first soil
void checkLinkage(Ground& g) {
  std::vector<Layer*> fw, bw;
  for (Layer* l = g.toplayer; l; l = l->nextl) fw.push_back(l);
  for (Layer* l = g.botlayer; l; l = l->prevl) bw.push_back(l);
  REQUIRE(fw.size() == bw.size());
  for (size_t i = 0; i < fw.size(); ++i) {
    CHECK(fw[i] == bw[fw.size() - 1 - i]);
    CHECK(fw[i]->indl == (int)i + 1);
    CHECK(fw[i]->prevl == (i ? fw[i - 1] : NULL));
    CHECK(fw[i]->nextl == (i + 1 < fw.size() ? fw[i + 1] : NULL));
  }
  CHECK(g.toplayer->prevl == NULL); CHECK(g.botlayer->nextl == NULL);
  bool soil = false;
  for (Layer* l : fw) { if (l->isSoil) soil = true; else if (l->isSnow) CHECK(!soil); }
  Layer* f = g.toplayer; while (f && !f->isSoil) f = f->nextl;
  CHECK(g.fstsoill == f);
}
double colH(Ground& g) { double h = 0; for (Layer* l = g.toplayer; l; l = l->nextl) h += Soil_Env::enthLayerH(l); return h; }
}

TEST_SUITE("snow_empty_cleanup") {

TEST_CASE("switch defaults to 0 (legacy) in ModelData and Ground; 1 is read; other values are FATAL") {
  ModelData m0; CHECK(m0.snow_empty_cleanup == 0);
  Json::Value c0; c0["model_settings"]["svp_ice_fix"] = 0; ModelData mj0(c0); CHECK(mj0.snow_empty_cleanup == 0);
  Json::Value c1; c1["model_settings"]["snow_empty_cleanup"] = 1; ModelData mj1(c1); CHECK(mj1.snow_empty_cleanup == 1);
  Json::Value c2; c2["model_settings"]["snow_empty_cleanup"] = 2;
  CHECK(ut::exit_status_of([&] { ModelData m(c2); }) == EXIT_FAILURE);
  Ground g; CHECK(g.snow_empty_cleanup == 0);
}

TEST_CASE("an empty snow layer between snow layers is removed; linkage, indices and mass are kept") {
  SnowLayer* a = snl(0.05, 10.0, 0.0, -3.0, 200.0);
  SnowLayer* e = snl(0.0, 0.0, 0.0, -1.0, 210.0);
  SnowLayer* d = snl(0.08, 16.0, 0.0, -5.0, 200.0);
  SnowCol c({a, e, d});
  const double W0 = c.allW();
  REQUIRE(c.nsnow() == 3);
  CHECK(c.col.g.removeEmptySnowLayers());
  c.col.g.resortGroundLayers();
  CHECK(c.nsnow() == 2);
  CHECK(c.col.g.toplayer == a); CHECK(a->nextl == d); CHECK(d->prevl == a);
  CHECK(d->nextl == c.col.L.front()); CHECK(c.col.L.front()->prevl == d);
  CHECK(c.allW() == W0);                                   // exact: nothing moved
  CHECK(a->ice == 10.0); CHECK(d->ice == 16.0); CHECK(a->dz == 0.05); CHECK(d->dz == 0.08);
  CHECK(c.col.g.snow_empty_removed == 1); CHECK(c.col.g.snow_empty_days == 1);
  checkLinkage(c.col.g);
  c.col.g.updateSnowHorizon();
  CHECK(c.col.g.snow.numl == 2); CHECK(c.col.g.snow.thick == doctest::Approx(0.13));
  CHECK_FALSE(c.col.g.removeEmptySnowLayers());            // idempotent
}

TEST_CASE("empty top layer, empty layer on the soil and a cascade of empties are all removed") {
  SnowLayer* e0 = snl(0.0, 0.0, 0.0, -1.0, 200.0);
  SnowLayer* a = snl(0.05, 10.0, 0.0, -3.0, 200.0);
  SnowLayer* e1 = snl(0.0, 0.0, 0.0, 0.0, 200.0);
  SnowLayer* e2 = snl(0.0, 0.0, 0.0, -2.0, 200.0);
  SnowCol c({e0, a, e1, e2});
  const double W0 = c.allW();
  CHECK(c.col.g.removeEmptySnowLayers());
  c.col.g.resortGroundLayers();
  CHECK(c.nsnow() == 1); CHECK(c.col.g.toplayer == a); CHECK(a->prevl == NULL);
  CHECK(a->nextl == c.col.L.front());
  CHECK(c.allW() == W0); CHECK(c.col.g.snow_empty_removed == 3); CHECK(c.col.g.snow_empty_days == 1);
  checkLinkage(c.col.g);
}

TEST_CASE("every snow layer empty: the soil becomes the top layer, snow counts go to zero") {
  SnowCol c({snl(0.0, 0.0, 0.0, -1.0, 200.0), snl(0.0, 0.0, 0.0, -1.0, 200.0)});
  const double W0 = c.allW();
  CHECK(c.col.g.removeEmptySnowLayers());
  c.col.g.resortGroundLayers();
  CHECK(c.nsnow() == 0); CHECK(c.col.g.toplayer == c.col.L.front()); CHECK(c.col.g.toplayer->isSoil);
  CHECK(c.allW() == W0); CHECK(c.col.g.snow.extramass == 0.0);
  checkLinkage(c.col.g);
  c.col.g.updateSnowHorizon();
  CHECK(c.col.g.snow.numl == 0); CHECK(c.col.g.snow.thick == 0.0); CHECK(c.col.g.snow.swe == 0.0);
}

TEST_CASE("residual mass below the threshold is conserved: to the snow layer above, else below, else extramass") {
  SnowLayer* a = snl(0.05, 10.0, 0.0, -3.0, 200.0);
  SnowLayer* e = snl(2.0e-12, 4.0e-10, 3.0e-10, -1.0, 200.0);   // < 1e-9 kg m-2 each, dz <= 1e-9 m
  SnowCol c({a, e});
  CHECK(c.col.g.removeEmptySnowLayers());
  CHECK(c.nsnow() == 1);
  CHECK(a->ice == 10.0 + 4.0e-10); CHECK(a->liq == 3.0e-10);
  CHECK(c.col.g.snow_empty_mass == doctest::Approx(7.0e-10));

  SnowLayer* t = snl(1.0e-12, 5.0e-10, 0.0, -1.0, 200.0);         // top, snow below
  SnowLayer* b = snl(0.08, 16.0, 0.0, -5.0, 200.0);
  SnowCol c2({t, b});
  CHECK(c2.col.g.removeEmptySnowLayers());
  CHECK(c2.col.g.toplayer == b); CHECK(b->ice == 16.0 + 5.0e-10);

  SnowCol c3({snl(1.0e-12, 5.0e-10, 0.0, -1.0, 200.0)});         // only snow layer
  CHECK(c3.col.g.removeEmptySnowLayers());
  CHECK(c3.nsnow() == 0); CHECK(c3.col.g.snow.extramass == 5.0e-10);
}

TEST_CASE("layers that are not empty are kept: thin snow with mass, or zero mass with thickness") {
  SnowCol c({snl(1.0e-8, 2.0e-6, 0.0, -1.0, 200.0), snl(0.01, 0.0, 0.0, -1.0, 200.0), snl(0.08, 16.0, 0.0, -5.0, 200.0)});
  CHECK_FALSE(c.col.g.removeEmptySnowLayers());
  CHECK(c.nsnow() == 3); CHECK(c.col.g.snow_empty_removed == 0);
}

TEST_CASE("daily snow update after a lower layer melted out: off keeps a zero-thickness empty layer, on removes it") {
  for (int solver = 0; solver <= 1; ++solver) {
    CAPTURE(solver);
    double Woff = 0, Won = 0, Hoff = 0, Hon = 0, moff = 0, mon = 0;
    for (int on = 0; on <= 1; ++on) {
      CAPTURE(on);
      // top keeps ice; the bottom layer melted completely in the thermal step (all its ice is now liq, frozen 0)
      SnowLayer* a = snl(0.05, 10.0, 0.0, -0.5, 200.0);
      SnowLayer* b = snl(0.08, 0.0, 3.0, 0.0, 200.0, 0);
      SnowCol c({a, b}, on, solver);
      const double W0 = c.allW();
      Soil_Env::EnthSnap es;
      if (solver == 1) c.se.enthSnap(es);
      const double H0 = colH(c.col.g);
      c.sn.updateDailyM(-5.0);
      if (solver == 1) c.se.enthRemapSnow(es, -5.0);
      const double melt = c.ed.d_snw2soi.melt;
      CHECK(melt == 3.0);                                   // the bottom layer's water left as melt
      CHECK(c.allW() + melt == doctest::Approx(W0).epsilon(1e-15));
      checkLinkage(c.col.g);
      if (on == 0) {
        REQUIRE(c.nsnow() == 2);                            // legacy: empty layer survives into the next day
        CHECK(b->ice == 0.0); CHECK(b->liq == 0.0); CHECK(b->dz == 0.0);
        CHECK(c.col.g.snow.numl == 2);
        Woff = c.allW(); Hoff = colH(c.col.g); moff = melt;
      } else {
        CHECK(c.nsnow() == 1); CHECK(c.col.g.toplayer == a);
        CHECK(a->nextl == c.col.L.front());
        CHECK(c.col.g.snow.numl == 1); CHECK(c.col.g.snow_empty_removed == 1);
        CHECK(c.ed.d_snws.swesum == doctest::Approx(a->ice));
        Won = c.allW(); Hon = colH(c.col.g); mon = melt;
      }
      if (solver == 1) {                                    // the remap books the snow-window change exactly
        CHECK(colH(c.col.g) - H0 == doctest::Approx(c.se.enth_bk_snow).epsilon(1e-12).scale(1.0));
      }
    }
    CHECK(Won == Woff); CHECK(mon == moff);                // same water and melt
    CHECK(Hon == doctest::Approx(Hoff).epsilon(1e-14));    // the empty layer held no heat
  }
}

TEST_CASE("new solver: next day's thermal step on the cleaned r1-like column has no zero-capacity snow cell and closes") {
  // r1 30-31 May 2024 analogue: a remnant top layer over a bottom layer that melted out today
  SnowLayer* top = snl(1.5749e-4, 9.020455e-3, 0.0, -0.009999366, 57.28);
  SnowLayer* bot = snl(0.08, 0.0, 2.0, 0.0, 209.87, 0);
  SnowCol c({top, bot}, 1, 1);
  Soil_Env::EnthSnap es; c.se.enthSnap(es);
  c.sn.updateDailyM(-0.5);
  c.se.enthRemapSnow(es, -0.5);
  checkLinkage(c.col.g);
  for (Layer* l = c.col.g.toplayer; l && l->isSnow; l = l->nextl) CHECK(l->ice + l->liq > 0.0);
  const double H0 = colH(c.col.g); double W0 = 0; for (Layer* l = c.col.g.toplayer; l; l = l->nextl) W0 += l->liq + l->ice;
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, 136.59, q);
  CHECK(q.fallbacks == 0);
  CHECK(fabs(q.res) / 86400.0 <= 1e-6);
  double W1 = 0; for (Layer* l = c.col.g.toplayer; l; l = l->nextl) W1 += l->liq + l->ice;
  CHECK(W1 == doctest::Approx(W0).epsilon(1e-14));
  CHECK(colH(c.col.g) - H0 == doctest::Approx(q.init + q.equil + q.gtop + q.gbot + q.src + q.res + q.clip).epsilon(1e-12).scale(1.0));
}

}
