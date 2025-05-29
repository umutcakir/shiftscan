from setuptools import setup, find_packages

setup(
    name="shiftscan",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
        "biotite>=1.3.0",
        "pandas>=2.2.3",
        "joblib>=1.5.1",
        "tqdm>=4.67.1"
    ],
    entry_points={
        "console_scripts": [
            "shiftscan=shiftscan.cli:run",
        ],
    },
    author="Umut Cakir & Ali Yurtseven",
    author_email="hpumut@gmail.com",
    description="ShiftScan CLI: shiftscan",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url= "https://github.com/umutcakir/shiftscan",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)