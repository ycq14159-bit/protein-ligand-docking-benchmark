#!/usr/bin/env python3
import subprocess,sys
raise SystemExit(subprocess.call([sys.executable,'/root/autodl-tmp/benchmark_1.0/filter_2_ligand_qualification/scripts/filter2_pipeline.py','preflight',*sys.argv[1:]]))
