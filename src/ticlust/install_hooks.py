from setuptools.command.install import install
from ticlust_helper.get_vsearch import get_vsearch_bin_path
from pathlib import Path
import sys

class CustomInstall(install):
    def run(self):
        self.get_vsearch_path()
        install.run(self)
    
    def get_vsearch_path(self):
        try:
            vsearch_bin_file = get_vsearch_bin_path()
            Path(vsearch_bin_file).chmod(0o777)
        except Exception as e:
            print(f"WARNING: Could not download the binary! {e}", file=sys.stderr)