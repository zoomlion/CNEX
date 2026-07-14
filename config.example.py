"""Default CNEX configuration — copy to config.py and edit."""
import os

ROOT = os.path.dirname(os.path.abspath(__file__))

# Internal binaries (relative to ROOT)
VALIDATE = os.path.join(ROOT, "bin", "validate")
ASSEMBLE = os.path.join(ROOT, "bin", "assemble")
MERTABLE = os.path.join(ROOT, "bin", "mertable")
PIGZ     = os.path.join(ROOT, "src", "pigz", "pigz")

# External tools (customize per machine — use name in PATH or full path)
FAMSA    = "famsa"
FastTree = "FastTree"
ASTRAL   = "astral"
IQTREE3  = "iqtree3"

# Pipeline thresholds
MIN_CNE_PER_SPECIES = 100
ALIGNMENT_JOBS = 4           # parallel FAMSA processes
DEFAULT_THREADS = 20         # default threads for ASTRAL
IQTREE_THREADS = 20          # IQ-TREE 3 thread count

# Phylogeny defaults
DEFAULT_METHOD = "concat"                    # "concat" or "astral"
CONCAT_IDENTITY_THRESHOLDS = "0.4,0.6,0.8"  # comma-separated identity thresholds; empty = single run
ASTRAL_IDENTITY_THRESHOLDS = "0.4,0.6,0.8"  # comma-separated identity thresholds; empty = single run
ELEMENT_TAGS_FILE = ""                       # path to element_tags.tsv (ele_id\\ttag with header), empty = no tags
DRY_RUN = False                              # True = only generate scripts, don't run
