from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            "spectrum.mydpss",
            [
                "src/cpp/mydpss.c",
            ],
            export_symbols=["multitap"],
        )
    ],
)
