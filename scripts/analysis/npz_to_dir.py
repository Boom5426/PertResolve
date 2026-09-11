"""Convert an extract written as .npz into the memory-mappable directory form.

One array at a time, so the peak is the largest single array rather than the archive.
Written because three extracts already exist as .npz and re-extracting them costs about
fifty minutes of streaming that the conversion does not.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        raise SystemExit("usage: npz_to_dir.py <in.npz> <out_dir>")
    source, target = Path(args[0]), Path(args[1])
    target.mkdir(parents=True, exist_ok=True)
    with np.load(source, allow_pickle=False) as blob:
        for key in ("data", "indices", "indptr", "shape", "labels"):
            array = blob[key]
            np.save(target / f"{key}.npy", array)
            print(f"  {key:10} {array.dtype}  {array.shape}  "
                  f"{array.nbytes / 1e9:.2f} GB", flush=True)
            del array
        (target / "provenance.json").write_text(str(blob["provenance"][0]))
    json.loads((target / "provenance.json").read_text())      # fail loudly if malformed
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
