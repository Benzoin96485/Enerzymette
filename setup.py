from setuptools import setup, find_packages

setup(
    name='enerzymette',
    version='0.0.1',
    install_requires=['ase', 'mendeleev'],
    entry_points={'console_scripts': ['enerzymette=enerzymette.cli:main']},
    packages=find_packages(include=["enerzymette", "enerzymette.*"]),
    auth='Benzoin96485',
    author_email='luowl7@mit.edu',
    description='Enerzyme workflow manager - Around Neural Biochemistry',
    url='https://github.com/Benzoin96485/Enerzymette',
    zip_safe=True,
)
