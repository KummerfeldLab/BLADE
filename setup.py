from setuptools import setup, find_packages

setup(
    name='artifactsRemoval',  # Replace 'mypackage' with the name of your package
    version='0.1',  # Replace '0.1.0' with the current version of your package
    author='Yinzhao Wang',  # Replace 'Your Name' with your name
    author_email='wyz2016828@example.com',  # Replace with your email
    description='A package for detect and removal of artifacts in Spatial Transcriptomics Tissue samples',  # A short description
    long_description=open('README.md').read(),  # Long description read from the README.md
    long_description_content_type='text/markdown',  # This is important for formatting the long description
    url='https://github.com/KummerfeldLab/artifactDetection',  # Link to your package's GitHub repo
    packages=find_packages(),  # Finds all python modules in the directory automatically
    install_requires=['numpy', 
                      'pandas',
                      'scipy',
                      'matplotlib',
                      ],
    keywords='artifacts removal',  # Keywords for your project
    python_requires='>=3.6'  # Minimum version requirement of the Python runtime
)
