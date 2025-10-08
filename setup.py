from pathlib import Path
from setuptools import setup, find_packages

ROOT = Path(__file__).resolve().parent
README_PATH = ROOT / "README.md"
LONG_DESCRIPTION = README_PATH.read_text(encoding="utf-8") if README_PATH.exists() else ""

setup(
    name="artifactsRemoval",
    version="0.2.0",
    description="Artifact removal workflows for Visium and Visium HD spatial transcriptomics data",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author="Kummerfeld Lab",
    author_email="kummerfeldlab@umich.edu",
    url="https://github.com/KummerfeldLab/BLADE",
    packages=find_packages(exclude=("old_ver", "tests", "tests.*")),
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "h5py",
        "pyarrow",
    ],
    keywords=["spatial transcriptomics", "artifact removal", "Visium"],
    include_package_data=True,
)
