import setuptools

#with open("README.rst", "r") as fh:
#    long_description = fh.read()
with open("requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh]

setuptools.setup(
    name="telo-permutations",
    version="1.0",
    author="Kristina Buss",
    author_email="kristina.buss@asu.edu",
    description="A Python tool to find length-limited permutations of a provided telomere repeat sequence within a given genome",
    #long_description=long_description,
    #long_description_content_type="text/x-rst",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=requirements,
)
