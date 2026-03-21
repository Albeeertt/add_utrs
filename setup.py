from setuptools import setup, find_packages

with open("requeriments.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="ss_utr",
    version="0.1.2",
    packages=find_packages(include=["ss_utr", "ss_utr.*"]),
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'ss_utr=ss_utr.cli:main',
        ],
    },
)
