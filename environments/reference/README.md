# Pinned reference runtimes

This recipe installs separate R and Python prefixes for the reference backend. It reuses the exact 244 R-prefix and 87 Python-prefix Conda package records validated by the isolated GBM fixture. These are runtime dependency closures, not a minimal requirements list. The locks target Linux x86-64; they do not claim portability to macOS, Windows or other CPU architectures.

```bash
bash environments/reference/install.sh /conda-base/bin/python /env/dgscrna-reference
/env/dgscrna-reference/python/bin/python -m pip install --no-deps \
  /downloads/dgscrna-2.0.0rc1-py3-none-any.whl
```

Run in a compute allocation on HPC. `examples/reference/install.sbatch` provides a scheduler-neutral template. The root must be fresh: existing `r`, `python` or `package_cache` directories cause an error. Installation uses copies in an isolated package cache, verifies exact package archives/builds, and records the actual loaded R/Python package paths. User/site Conda configuration is excluded without editing it.

The install script downloads the original `dbscan_1.2.6.tar.gz` source from CRAN and checks SHA256 `b4eab5a7ec4bdc2b65a1e89ac1cdd5860fe64f407953f19fc8e0af2922f0c9a6`. Set `DGSCRNA_DBSCAN_SOURCE` to an existing copy if required. It then applies the recorded index-width correction, builds version `1.2.6.9001`, and saves its source and installation manifest under the new prefix. No patched binary is bundled. See `THIRD_PARTY_NOTICES.md`.

Files:

- `*-linux-64.explicit.txt`: exact Conda URLs and MD5 identifiers.
- `*.provenance.json`: version/build records and SHA256 where supplied by Conda; site-specific prefix names removed.
- `create_isolated.py`: Conda API installation with an empty configuration search path.
- `fetch_dbscan.py`, `install_dbscan_index64.R`: checked source acquisition and the explicit patch/build.
- `check_isolation.{py,R}`: loaded-package path verification and a small CPU Torch smoke check.
- `verify_downloads.py`: pristine archive checksums and installed build verification.

The successful historical fixture verified all 331 archive/build records and matched fresh fitted outputs, including DL probabilities and terminal labels. That result certifies the fixture and its recorded environment. Run the package's reference tests in the target environment before extending this conclusion to other inputs.
