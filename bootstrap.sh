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
# Download required SPICE kernels
#
echo "## Downloading SPICE kernels:"
download oif/data/latest_leapseconds.tls https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls
download oif/data/earth_070425_370426_predict.bpc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/earth_070425_370426_predict.bpc
download oif/data/earth_latest_high_prec.bpc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
download oif/data/de430.bsp https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
download oif/data/earth_topo_201023.tf https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/stations/earth_topo_201023.tf
echo "done."
