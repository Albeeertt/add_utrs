from setuptools import setup, find_packages

with open("requeriments.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="add_utrs",
    version="1.0.1",
    packages=find_packages(include=["add_utrs", "add_utrs.*"]),
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'add_utrs=add_utrs.cli:main',
        ],
    },
)
