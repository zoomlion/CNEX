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
IQTREE_THREADS = 20
