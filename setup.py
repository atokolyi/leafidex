from setuptools import setup

setup(name='TheSpliceGirls',
      version='0.1',
      description='An annotation pipeline for transcript splicing events',
      url='https://github.com/atokolyi/TheSpliceGirls',
      author='Alex Tokolyi',
      author_email='alex.tokolyi@gmail.com',
      license='gplv3',
      packages=['TheSpliceGirls'],
      install_requires=[
          'pandas','pyranges1','requests','numpy'
      ],
      zip_safe=False)
