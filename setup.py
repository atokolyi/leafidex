from setuptools import setup

setup(name='TheSpliceGirls',
      version='0.1',
      description='An annotation pipeline for transcript splicing events',
      url='https://github.com/atokolyi/TheSpliceGirls',
      author='Alex Tokolyi',
      author_email='alex.tokolyi@gmail.com',
      license='gplv3',
      packages=['TheSpliceGirls'],
      package_data={'TheSpliceGirls' :['TheSpliceGirls/data/*']},
      install_requires=[
          'pandas','pyranges1','requests'
      ],
      zip_safe=False)
