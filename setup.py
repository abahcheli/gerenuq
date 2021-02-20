from setuptools import setup, find_packages
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# # automatically captured required modules for install_requires in requirements.txt
# with open(path.join(HERE, 'requirements.txt'), encoding='utf-8') as f:
#     all_reqs = f.read().split('\n')

# install_requires = [x.strip() for x in all_reqs if ('git+' not in x) and (not x.startswith('#')) and (not x.startswith('-'))]
# dependency_links = [x.strip().replace('git+', '') for x in all_reqs \
#                     if 'git+' not in x]

setup(
    name = "gerenuq",
    version = "0.2.0",
    author = "Alec Bahcheli",
    author_email = "abahchel@uwo.ca",
    description = "Samfile long-read filtering script.",
    long_description = README,
    long_description_content_type="text/markdown",
    url="https://github.com/abahcheli/gerenuq",
    packages = find_packages(),
    install_requires = ["pysam", "pandas"],
    python_requires = '>=3.0',
    # scripts = ['gerenuq/gerenuq_cl.py'],
    entry_points = {
        'console_scripts': ['gerenuq=gerenuq.gerenuq_cl:main']
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Academic Free License (AFL)",
        "Operating System :: OS Independent",
    ]
)

