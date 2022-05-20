# python setup.py bdist_wheel
# pip install

from setuptools import setup, find_packages
import pathlib, os

with open("README.md", "r") as fh:
    long_description = fh.read()




setup(

    name = 'rovar',
    author="Anders Lindanger Roevaer",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/andersbhm/useful_functions_phd",

    version='0.1',
    classifiers=[
             "Programming Language :: Python :: 3",
             "License :: OSI Approved :: MIT License"],
    packages=find_packages(), #include/exclude arguments take * aswildcard, . for any sub-package names
)
