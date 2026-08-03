"""Regenerate every independent refined Fig. 5 panel."""
from __future__ import annotations

import importlib

MODULES = [
    "fig5a_replicate_ceiling",
    "fig5b_ceiling_vs_models",
    "fig5c_phase_map",
    "fig5d_predictor_family",
    "fig5e_rank_recovery_concept",
    "fig5f_external_benchmarks",
    "fig5g_recovery_support",
]


def main() -> None:
    for name in MODULES:
        print(f"rendering {name}")
        importlib.import_module(name).main()


if __name__ == "__main__":
    main()
