#!/usr/bin/env python3
import subprocess,sys
raise SystemExit(subprocess.call([sys.executable,'/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification_v2/scripts/filter2_v2_pipeline.py','selftest',*sys.argv[1:]]))
