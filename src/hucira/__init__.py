import logging
from importlib.metadata import version

from .data import load_CIP_signatures, load_cytokine_info, load_human_cytokine_dict, load_Lupus_data, load_MS_CSF_data
from .pl import plot_communication, plot_significant_results
from .tl import (
    check_robustness,
    create_celltype_combos,
    get_all_senders_and_receivers,
    get_one_senders_and_receivers,
    get_robust_significant_results,
    run_all_enrichment_test,
    run_one_enrichment_test,
)

# Configure library-level logger so INFO messages are visible by default.
logging.getLogger(__name__).addHandler(logging.StreamHandler())
logging.getLogger(__name__).setLevel(logging.INFO)

__all__ = [
    "load_cytokine_dict_data",
    "load_MS_data",
    "load_cytokine_info",
    "load_CIP_signatures",
    "check_robustness",
    "create_celltype_combos",
    "get_robust_significant_results",
    "run_one_enrichment_test",
    "run_all_enrichment_test",
    "plot_significant_results",
    "plot_communication",
    "get_one_senders_and_receivers",
    "get_all_senders_and_receivers",
]


__version__ = version("hucira")
