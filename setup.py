import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="symcirc",
    version="0.0.3",
    author="Matyas Vasek",
    author_email="matyas.vasek@gmail.com",
    description="A python package for symbolic electronic circuit analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT license",
    url="https://github.com/MatyasVasek/SymCirc",
    project_urls={
        "Bug Tracker": "https://github.com/MatyasVasek/SymCirc/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    install_requires=['sympy', ],
)