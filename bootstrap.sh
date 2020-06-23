#!/bin/bash

# stop on any error
set -e

# utilities
download()
{
	# download a file if it doesn't exist.
	# download <destination> <URL>

	local F="$1"
	local URL="$2"
	if ! test -f "$F"; then
		echo "$F: "
		curl -L# -o "$F" "$URL"
	else
		echo "$F exists; skipping download."
	fi
}

#
# Download required SPICE utilities
#
echo "## Downloading SPICE utilities:"
mkdir -p spice
for F in msopck bspidmod pinpoint; do
	download "spice/$F" "http://naif.jpl.nasa.gov/pub/naif/utilities/MacIntel_OSX_64bit/$F"
	chmod +x "spice/$F"
done
echo

#
# Download required SPICE kernels
#
echo "## Downloading SPICE kernels:"
mkdir -p kernels
download kernels/latest_leapseconds.tls https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls
download kernels/earth_070425_370426_predict.bpc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/earth_070425_370426_predict.bpc
download kernels/earth_latest_high_prec.bpc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
download kernels/de430.bsp https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
download kernels/earth_topo_050714.tf https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/stations/earth_topo_050714.tf
echo

#
# Locate where OpenOrb data may be, suggest what to use for OORB_DATA variable
#
OORB=$(command -v oorb)
if [[ ! -z $OORB ]]; then
	OORB_DATA="$(dirname $(dirname $OORB))/share/openorb"
	if [[ -f "$OORB_DATA/de430.dat" ]]; then
		echo "Found OpenOrb data files. Set OORB_DATA by running:"
		echo '   export OORB_DATA="'"$OORB_DATA"'"'
		echo "(or the csh equivalent, if that's what you're using)".
	else
		OORB_DATA=;
	fi
fi
if [[ -z "OORB_DATA" ]]; then
	echo "Make sure to set OORB_DATA to wherever your openorb instance is."
fi
