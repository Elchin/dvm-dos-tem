import re
p='src/Soil_Env_enthalpy.cpp'; s=open(p).read()
s=s.replace('#include <iostream>\n','#include <iostream>\n#include <cstdlib>\n#include <fstream>\n',1)
anchor='  const double phys = I * h_in - out * hp;\n'
dbg='''  { static const char* dbgf = getenv("ENTH_DBG_FILE");   // DIAG ONLY (not for commit)
    static long dcount = 0; dcount++;
    if (dbgf != NULL && fabs(columnEnthalpyE() - s.hcol) > 5.0e5) {
      std::ofstream f(dbgf, std::ios::app);
      f << "WATERWIN call=" << dcount << " dH=" << (columnEnthalpyE() - s.hcol) << " I=" << I << " P=" << P << " G=" << G << " hp=" << hp
        << " drainl=" << (ground->drainl ? ground->drainl->indl : -1) << " shlwfrz=" << (ground->fstshlwl ? ground->fstshlwl->frozen : -9)
        << " mossfrz=" << (ground->fstmossl ? ground->fstmossl->frozen : -9) << " ststate=" << ground->ststate << "\\n";
      for (int i = 0; i < n; ++i) { Layer* l = now[i];
        f << "  " << i << " t=" << (l->isSnow?"S":(l->isRock?"R":(l->isMineral?"M":"O"))) << " z=" << l->z << " dz=" << l->dz << " T0=" << s.T[i] << " T1=" << l->tem
          << " liq0=" << s.liq[i] << " liq1=" << l->liq << " ice0=" << s.ice[i] << " ice1=" << l->ice << " frz=" << l->frozen << " ff=" << l->frozenfrac
          << " poro=" << l->poro << " maxice=" << l->maxice << " maxliq=" << l->maxliq << "\\n"; }
    } }
'''
assert anchor in s; s=s.replace(anchor,dbg+anchor,1)
open(p,'w').write(s)
