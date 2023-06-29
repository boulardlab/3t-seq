import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="utilities-tabaro",
    version="0.0.1",
    author="Francesco Tabaro",
    author_email="francesco.tabaro@embl.it",
    description="Utilities and convenience functions for Snakemake rules and random stuff - mainly personal",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
