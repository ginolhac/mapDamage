#!/usr/bin/env python3
from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            "mapdamage.seqtk",
            sources=["mapdamage/seqtk/seqtk.c"],
            libraries=["z"],
        )
    ],
)
