#include "../include/Organic.h"
#include <cmath>

Organic::Organic() : humic_max_thick(0.0) {
};

Organic::~Organic() {
};

void Organic::ShlwThickScheme(const double & totthickness) {
  shlwthick = totthickness;

  if (totthickness<0.) {
    shlwthick = 0.;
  }

  shlwnum   = 0;

  //These values are from Shuhua Yi's paper, 2010
  if(totthickness<=0.00) {
    shlwdz[0]=MISSING_D;
    shlwdz[1]=MISSING_D;
    shlwdz[2]=MISSING_D;
    shlwnum =0;
  } else if (totthickness<0.04) {
    shlwdz[0]=totthickness;
    shlwdz[1]=MISSING_D;
    shlwdz[2]=MISSING_D;
    shlwnum =1;
  } else if (totthickness<0.06) {
    shlwdz[0]=0.02;
    shlwdz[1]=totthickness-0.02;
    shlwdz[2]=MISSING_D;
    shlwnum =2;
  } else if (totthickness<0.10) {
    shlwdz[0]=0.03;
    shlwdz[1]=totthickness-0.03;
    shlwdz[2]=MISSING_D;
    shlwnum =2;
  } else if (totthickness<0.15) {
    shlwdz[0]=0.02;
    shlwdz[1]=0.04;
    shlwdz[2]=totthickness-0.06;
    shlwnum =3;
  } else if (totthickness<0.20) {
    shlwdz[0]=0.03;
    shlwdz[1]=0.06;
    shlwdz[2]=totthickness -0.09;
    shlwnum =3;
  } else if (totthickness<0.28) {
    shlwdz[0]=0.04;
    shlwdz[1]=0.08;
    shlwdz[2]=totthickness -0.12;
    shlwnum =3;
  } else if (totthickness<0.40) {
    shlwdz[0]=0.05;
    shlwdz[1]=0.11;
    shlwdz[2]=totthickness -0.16;
    shlwnum =3;
  } else if (totthickness >=0.4) {
    shlwdz[0]=0.1;
    shlwdz[1]=0.2;
    shlwdz[2]=totthickness -0.3;
    shlwnum =3;
  }

  if (shlwnum>0) {
    lstshlwdz = shlwdz[shlwnum-1];
  }
};



void Organic::DeepThickScheme(const double & totthickness) {
  deepthick = totthickness;

  if(totthickness<=0.0) {
    deepthick =0.0;
    deepdz[0]=MISSING_D;
    deepdz[1]=MISSING_D;
    deepdz[2]=MISSING_D;
    deepnum =0;
    return;
  }

  deepnum = 0;

  if(lstshlwdz>0.0) {
    if(deepthick < 3* lstshlwdz) {
      deepdz[0]=deepthick;
      deepdz[1]=MISSING_D;
      deepdz[2]=MISSING_D;
      deepnum =1;
    } else if(deepthick >= 3* lstshlwdz && deepthick<6*lstshlwdz) {
      deepdz[0]=1./3. * deepthick;
      deepdz[1]=2./3. * deepthick;
      deepdz[2]=MISSING_D;
      deepnum =2;
    } else {
      deepdz[0]=1./6. * deepthick;
      deepdz[1]=2./6. * deepthick;
      deepdz[2]=3./6. * deepthick;
      deepnum =3;
    }
  } else {
    if(deepthick <= 0.02) {
      deepdz[0]=deepthick;
      deepdz[1]=MISSING_D;
      deepdz[2]=MISSING_D;
      deepnum =1;
    } else if(deepthick<=0.06) {
      deepdz[0]=1./3. * deepthick;
      deepdz[1]=2./3. * deepthick;
      deepdz[2]=MISSING_D;
      deepnum =2;
    } else {
      deepdz[0]=1./6. * deepthick;
      deepdz[1]=2./6. * deepthick;
      deepdz[2]=3./6. * deepthick;
      deepnum =3;
    }
  }

  // v5q-humic-stefan: optional max-thickness cap on the deep (humic) sublayers (0 = off, legacy)
  if (humic_max_thick > 0.0) {
    applyHumicMaxThick(humic_max_thick, deepthick, deepdz, deepnum);
  }
};

// v5q-humic-stefan: equal re-split of the deep horizon so no sublayer exceeds 'cap' (if feasible
// within MAX_DEP_LAY). No-op when cap <= 0, deepthick <= 0, or every current sublayer is <= cap.
// Refs: humic layering and the organic C<->thickness relation of the dynamic organic soil module:
//   Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010, doi:10.1029/2010JG001302.
//   The cap itself is a vertical-resolution (discretization) choice: no published source.
bool Organic::applyHumicMaxThick(double cap, double deepthick, double deepdz[], int & deepnum) {
  if (!(cap > 0.0) || !(deepthick > 0.0) || deepnum < 1) { return false; }
  double maxdz = 0.0;
  for (int k = 0; k < deepnum && k < MAX_DEP_LAY; k++) { maxdz = fmax(maxdz, deepdz[k]); }
  if (maxdz <= cap * (1.0 + 1e-12)) { return false; }
  int n = (int)std::ceil(deepthick / cap - 1e-12);
  if (n < deepnum) { n = deepnum; }
  if (n > MAX_DEP_LAY) { n = MAX_DEP_LAY; }
  for (int k = 0; k < MAX_DEP_LAY; k++) { deepdz[k] = (k < n) ? deepthick / n : MISSING_D; }
  deepnum = n;
  return true;
}

// if shlw peat thickness from input soil profile
void Organic::assignShlwThicknesses(int soiltype[], double soildz[],
                                    const int & soilmaxnum) {
  shlwnum   = 0;
  shlwthick = 0;

  for(int i=0; i<soilmaxnum; i++) {
    if(soiltype[i] ==1) {
      shlwdz[shlwnum] = soildz[i];
      shlwnum ++;
      shlwthick += soildz[i];
    } else {
      if(soiltype[i]>2) {
        break;
      }
    }
  }
};

// if deep peat thickness from input soil profile
void Organic::assignDeepThicknesses(int soiltype[], double soildz[],
                                    const int & soilmaxnum) {
  deepnum =0;
  deepthick =0;

  for(int i=0; i<soilmaxnum; i++) {
    if(soiltype[i] ==2) {
      deepdz[deepnum] = soildz[i];
      deepnum ++;
      deepthick += soildz[i];
    } else {
      if(soiltype[i]>2) {
        break;
      }
    }
  }
};

