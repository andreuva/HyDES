from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(

    name = "hydes",
    version = "0.0.2",
    author = "Andres Vicente Arevalo",
    author_email = "andres.vicente.arevalo@gmail.com",
    description = "Hydrodinamical equation solver",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/andreuva/HyDES",
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/andreuva/HyDES/issues',
        'Funding': 'https://donate.pypi.org',
        'Source': 'https://github.com/andreuva/HyDES',
    },
    package_dir={'': 'src'},  # Optional
    classifiers = [
                "Programming Language :: Python :: 3",
                "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent"
                ],
    packages=find_packages(where='src'),
    python_requires = ">=3.6",
    install_requires=['numpy', 'matplotlib'],
    extras_require={'movie': ['opencv-python']},
    package_data={'sample_params': ['params.json']},

)