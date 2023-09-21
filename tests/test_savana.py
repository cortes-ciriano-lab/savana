"""
Testing module for SAVANA
Created: 10/01/2023
Python 3.9.6
Hillary Elrick
"""
#!/usr/bin/env python3

import subprocess

from pathlib import Path

ROOTDIR = Path(__file__).parent.parent

def test_install():
	""" test SAVANA is installed successfully """
	cmd = f"python3 {ROOTDIR}/savana/savana.py --help"
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	_, err = p.communicate()
	if p.returncode != 0:
		raise RuntimeError(f"FAILED: {cmd}\n{err}")
