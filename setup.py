from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'SANTON'
LONG_DESCRIPTION = 'Sequencing Analysis Toolkits for Off-target Nomination'

setup(
        name="santon", 
        version=VERSION,
        author="Huan Qiu, Ruijia Wang",
        author_email="help.qbio@vorbio.com",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
		license='MIT',
        packages=find_packages(),
        keywords=['python', 'santon'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Programming Language :: Python :: 3",
        ],
        zip_safe=False
)

