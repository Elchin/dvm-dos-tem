#!/bin/bash
# mkvar.sh NAME SUBSTEPS [extra "key": val lines...]  -> cs/NAME from cs/HPCW_s24src (config edits only)
W=~/tem-dev/ssf-work; n=$1; s=$2; shift 2
rm -rf $W/cs/$n; cp -a $W/cs/HPCW_s24src $W/cs/$n
for st in main spinup; do sed -i "s/\"thermal_substeps\": 24/\"thermal_substeps\": $s/" $W/cs/$n/$st/config/config.js
  for kv in "$@"; do sed -i "s/\"thermal_substeps\": $s,/\"thermal_substeps\": $s, $kv,/" $W/cs/$n/$st/config/config.js; done; done
grep -h -o '"thermal_substeps": [0-9]*[^"]*"[a-z_]*' $W/cs/$n/*/config/config.js
