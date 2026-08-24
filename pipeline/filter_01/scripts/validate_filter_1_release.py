#!/usr/bin/env python3
import subprocess,sys
raise SystemExit(subprocess.call([sys.executable, '/root/autodl-tmp/benchmark_1.0/filter_1_protein_receptor_qualification/scripts/filter1_pipeline.py', 'validate', *sys.argv[1:]]))
