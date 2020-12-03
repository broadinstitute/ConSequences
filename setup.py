from setuptools import setup

setup(name='ConSequences',
      version='0.1',
      description='Suite to delineate contiguous and conserved sequences from assemblies and search for their presence in raw sequencing data.',
      url='http://github.com/broadinstitute/ConSequences/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['ConSequences'],
      scripts=['bin/delineateSegmentsOnReference.py', 'bin/generateReferenceMSA.py', 'bin/querySegmentInRawReads.py'],
      zip_safe=False)