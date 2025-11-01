from setuptools import setup

setup(
    install_requires = [
        f"oxDNA_analysis_tools @ file://localhost/{os.getcwd()}/analysis/"
    ],
)
