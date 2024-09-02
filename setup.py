from setuptools import setup, find_packages
import os

def readme():
    with open('README.md', 'r') as fh:
        return fh.read()

this_directory = os.path.abspath(os.path.dirname(__file__))

if os.path.exists(os.path.join(this_directory, 'requirements.txt')):
    with open(os.path.join(this_directory, 'requirements.txt'), 'r') as f:
        requirements = [line.split('#')[0].strip() for line in f.readlines()]
        requirements = [line for line in requirements if line]
else:
    requirements = []

setup(
    name="Syndirella",
    version="2.0.1-alpha",
    author="Kate Fieseler",
    author_email="kate.fieseler@stats.ox.ac.uk",
    description="Synthetically directed elaborations",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/kate-fie/syndirella",
    packages=find_packages(),
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'syndirella=syndirella.cli:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    python_requires='>=3.10',
    include_package_data=True,
    zip_safe=False,
)
