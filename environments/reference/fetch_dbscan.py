"""Obtain and checksum the unmodified CRAN dbscan 1.2.6 source archive."""
from pathlib import Path
import argparse
import hashlib
import shutil
import urllib.error
import urllib.request

SHA256 = "b4eab5a7ec4bdc2b65a1e89ac1cdd5860fe64f407953f19fc8e0af2922f0c9a6"
URLS = (
    "https://cran.r-project.org/src/contrib/dbscan_1.2.6.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/dbscan/dbscan_1.2.6.tar.gz",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--source", type=Path, help="Use an already downloaded upstream archive")
    args = parser.parse_args()
    dest = args.destination.resolve()
    if dest.exists():
        parser.error("Destination exists; use a fresh file")
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".part")
    if args.source:
        shutil.copyfile(args.source, partial)
    else:
        for index, url in enumerate(URLS):
            try:
                with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as stream:
                    shutil.copyfileobj(response, stream)
                break
            except (urllib.error.URLError, TimeoutError):
                partial.unlink(missing_ok=True)
                if index == len(URLS) - 1:
                    raise
    with partial.open("rb") as stream:
        found = hashlib.file_digest(stream, "sha256").hexdigest()
    if found != SHA256:
        partial.unlink()
        raise RuntimeError(f"Unexpected dbscan source checksum: {found}")
    partial.replace(dest)
    print(f"VERIFIED_DBSCAN_SOURCE {dest} {SHA256}")


if __name__ == "__main__":
    main()
