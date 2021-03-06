from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(

    name = "hydes",
    version = "0.0.4",
    author = "Andres Vicente Arevalo",
    author_email = "andres.vicente.arevalo@gmail.com",
    description = "Hydrodinamical equation solver",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/andreuva/HyDES",
    project_urls={
        'Bug Reports': 'https://github.com/andreuva/HyDES/issues',
        'Funding': 'https://donate.pypi.org',
        'Source': 'https://github.com/andreuva/HyDES',
    },
    package_dir={'': 'src'},
    classifiers = [
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent"
                ],
    py_modules=['hydes', 'calc', 'domain', 'display', 'initCond', 'state', 'animation'],
    python_requires = ">=3.6",
    install_requires=['numpy', 'matplotlib', "opencv-python"],
    extras_require={'dev': ["setuptools>=54", "wheel", "twine"]},
    package_data={'sample_params': ['params.json']},
)