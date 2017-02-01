from setuptools import setup

setup(name='meshmonk',
      version='0.0.1',
      description='Algorithms base on scipy',
      url='https://bitbucket.org/webmonks/kuleuven-algorithms',
      author='Jonatan Snyders',
      author_email='jonatan@webmonks.io',
      license='MIT',
      packages=['meshmonk'],
      install_requires=[
            "numpy==1.7.1",
            "scipy==0.18.1"
      ],
      zip_safe=False)
