"""Python setup.py for ecg_hrv_analysis"""

import io
import os

from setuptools import find_packages, setup


def read(*paths, **kwargs):
    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="ecg_hrv_analysis",
    description="project_description",
    url="https://github.com/Maria-Iosif/ECG-HRV-Analysis",
    long_description=read("README.md"),
    author="Maria Iosif",
    packages=find_packages(exclude=[".github"]),
    install_requires=read_requirements("requirements.txt"),
)
