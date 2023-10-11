from setuptools import setup

setup(
    name="Connect_X",
    version="0.1",
    description="Application that efficiently connects genes to associated symptoms based on user-provided search parameters",
    author="Natalie Lipieta",
    author_email="natalielipieta@gmail.com",
    url="https://github.com/nlipieta/connect-x",
    scripts=["connect_x.py"],
    install_requires=['pandas', 'requests'],
)