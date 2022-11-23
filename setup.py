#!/usr/bin/env python

from setuptools import find_packages, setup

with open("requirements.txt") as fh:
    install_requires = fh.readlines()

with open("README.md") as fh:
    long_description = fh.read()

# following src dir layout according to
# https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
version = "0.1.0"
setup(
    name="cell-annotation-service-client",
    version=version,
    description="Cell Annotation Service Client",
    author="Variants Team",
    author_email="variants@broadinstitute.org",
    long_description=long_description,
    install_requires=install_requires,
    tests_require=["coverage", "pytest"],
    python_requires=">=3.7",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    include_package_data=True,
)
