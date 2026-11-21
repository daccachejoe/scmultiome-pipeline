#!/bin/bash
# MACS2 wrapper for the ArchR peaks branch (02_call_peaks.R passes this as
# ArchR's pathToMacs2). ArchR calls it as `<wrapper> callpeak ...`.
#
# Why a wrapper: on Minerva, loading the macs/2.1.0 module (python 2.7 and an
# old libstdc++) inside the R job breaks Rcpp. Loading it here keeps it in the
# MACS2 subprocess only.
#   archr_macs2_path set   -> exec that macs2 binary directly (e.g. a conda env)
#   otherwise              -> module purge + module load $archr_macs2_module
# Both keys are in config/pipeline.config and reach this script as env vars.

if [ -n "${archr_macs2_path:-}" ]; then
    exec "$archr_macs2_path" "$@"
fi

if ! type module >/dev/null 2>&1; then
    for init in "${LMOD_PKG:-/usr/share/lmod/lmod}/init/bash" /usr/share/lmod/lmod/init/bash /etc/profile.d/modules.sh; do
        [ -f "$init" ] && { source "$init"; break; }
    done
fi
module purge >/dev/null 2>&1
module load "${archr_macs2_module:-macs/2.1.0}" >/dev/null 2>&1
exec macs2 "$@"
