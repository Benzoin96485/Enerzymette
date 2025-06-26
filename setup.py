from setuptools import setup, find_packages

setup(
    name='enerzymette',
    version='0.0.1',
    install_requires=['ase'],
    entry_points={'console_scripts': ['enerzymette=enerzymette.cli:main']},
    packages=find_packages(include=["enerzymette", "enerzymette.*"]),
    auth='Benzoin96485',
    author_email='luowl7@mit.edu',
    description='A collection of small scripts that are useful for Enerzyme',
    url='https://github.com/Benzoin96485/Enerzymette',
    zip_safe=True,
)
