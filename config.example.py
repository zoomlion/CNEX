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
THREADS = 20                 # global worker count: FAMSA/FastTree parallelism, IQ-TREE/ASTRAL threads

# Phylogeny defaults
DEFAULT_METHOD = "concat"                    # "concat" or "astral"
CONCAT_LENGTH_QUANTILES = "25,50,75"  # comma-separated length quantiles; empty = single run
ASTRAL_BLOCK_GAPS = "1000,2000"       # kb thresholds for astral block clustering; comma-separated
GFF_DIR = ""                           # path to GFF3/*.gff* directory (for classify_elements.py / astral)
ELEMENT_TAGS_FILE = ""                       # path to element_tags.tsv (ele_id\\ttype, etc), empty = no tags
PARTITION = False                            # True = output partition file for IQ-TREE
DRY_RUN = False                              # True = only generate scripts, don't run
