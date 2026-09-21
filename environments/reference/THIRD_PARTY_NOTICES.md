# Third-party source used by the environment recipe

The recipe downloads dependencies from their original distributions. Their licenses and copyright notices remain with those distributions. Conda archives and third-party compiled binaries are not included in this directory.

## R dbscan 1.2.6

Upstream project: https://github.com/mhahsler/dbscan

Source archive: https://cran.r-project.org/src/contrib/Archive/dbscan/dbscan_1.2.6.tar.gz

SHA256: `b4eab5a7ec4bdc2b65a1e89ac1cdd5860fe64f407953f19fc8e0af2922f0c9a6`

The upstream `DESCRIPTION` states `GPL (>= 2)`. The affected `src/mst.cpp` carries the notice “Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.” and identifies its license as GNU GPL version 3. The package also contains ANN code attributed to the University of Maryland, Sunil Arya and David Mount. The source archive's full notices must be preserved.

DG-scRNA's `install_dbscan_index64.R` applies one code change in `src/mst.cpp`: both lower-triangle indices passed to `LT_POS0` are cast to `R_xlen_t` before integer products. It also changes the package version from `1.2.6` to `1.2.6.9001` to identify the modification. The algorithm and clustering parameters are unchanged. The installer records the original archive checksum, modified source checksum, destination library and build version.

The installer retains the downloaded source archive and modified source tree in the user's runtime prefix. This repository distributes the patch recipe rather than a modified binary. Anyone redistributing the resulting environment must retain the upstream notices and meet the corresponding source obligations of the included licenses. DG-scRNA's GPLv3 text is in the repository root `LICENSE`.
