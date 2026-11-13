# oxpy instructions

This file provides guidance to coding agents working on `src/oxpy/`, the pybind11 Python bindings of the oxDNA engine. The instructions in the top-level `AGENTS.md` also apply.

## What oxpy is

`oxpy` exposes the C++ engine to Python: users drive simulations from scripts (`OxpyManager`), read and modify the configuration while it runs (`ConfigInfo`, particles, forces), and write observables in Python. `oat` (in `analysis/`) depends on it. The user documentation is in `docs/source/oxpy/` and is built (by Sphinx, from the docstrings of the compiled module) at https://lorenzo-rovigatti.github.io/oxDNA/oxpy/index.html.

## Layout

- `bindings.cpp`: the `PYBIND11_MODULE(core, m)` definition. It calls one `export_<Class>(m)` function per exposed class, registers `oxDNAException` as the Python exception `OxDNAError`, and creates the `analysis` submodule.
- `bindings_includes/`: one header per exposed class (`BaseParticle.h`, `ConfigInfo.h`, `input_file.h`, ...), each defining its `export_<Class>` function and, for classes that can be subclassed from Python, a trampoline class. External forces are in `bindings_includes/Forces/` (the `forces` submodule) and observables in `bindings_includes/Observables/`.
- `OxpyContext.{h,cpp}`: the `oxpy.Context` context manager, which initialises and tears down the engine's global state (logger, `CONFIG_INFO`, ...). Every simulation must be created inside a `with oxpy.Context():` block.
- `OxpyManager.{h,cpp}`: `OxpyManager`, a `SimManager` subclass with the methods that Python needs (`run`, `system_energy`, `update_temperature`, `add_output`, ...).
- `python_defs.h`: pybind11 includes and the `py` namespace alias; include it first in every binding header.
- `vector_matrix_casters.h`: type casters that convert `LR_vector`/`LR_matrix` to and from numpy arrays.
- `__init__.py`, `utils.py`: the pure-Python part of the package (`oxpy.utils` helpers such as `generate_default_input`, `Kelvin_to_oxDNA`). `__init__.py` loads the extension with `RTLD_GLOBAL` so that plugins can see the engine's symbols: don't remove it.
- `pyproject_cmake.toml`: the `pyproject.toml` of the package installed by `make install` (the version is substituted from the git tag). The top-level `pyproject.toml` builds the same package with scikit-build-core (`pip install .` from the repository root).
- `pybind11/`: a vendored copy of pybind11. Do not modify it; update it only as a whole, in a dedicated pull request.

## Build and test

```bash
mkdir -p build && cd build
cmake .. -DPython=ON          # add -DOxpySystemInstall=ON to install into the active environment instead of with --user
make -j8 && make install      # make install runs pip on build/python/, which contains the staged package
python ../examples/OXPY/run.py   # smoke test (also run by CI)
make test_oxpy                # runs the oxpy level of test/TestSuite.py (e.g. test/OXPY/oxpy_input.py)
```

CMake refuses to build with a default Python older than 3.11 (`-DOxpyVersionOverride=ON` skips the check). It uses the `python`/`python3` found first in `PATH`, so activate the target environment before running CMake. The extension links against the same `oxdna_common` library as the executables, so any change to the engine requires rebuilding and reinstalling oxpy before testing it from Python. After installing, `oat`'s tests (`pytest analysis/tests/test_cli/`) also exercise oxpy.

## Adding or changing bindings

- New public C++ API meant for Python goes in a header in `bindings_includes/` (in `Forces/` or `Observables/` for those), with an `export_<Class>(py::module &m)` function that is called from `bindings.cpp` (or from the `export_` function of the parent submodule, as `export_BaseForce` does for the forces).
- Bind the base classes before the classes that derive from them, and pass the base class to `py::class_` (e.g. `py::class_<DNANucleotide, BaseParticle, std::shared_ptr<DNANucleotide>>`). All the bound classes use `std::shared_ptr<T>` as holder type: keep this convention, since a class bound with different holder types in different places compiles but crashes at runtime.
- Classes that users subclass in Python (observables, forces) need a trampoline class with `PYBIND11_OVERLOAD`/`PYBIND11_OVERLOAD_PURE` for every virtual method that Python may override (see `PyBaseObservable` in `bindings_includes/Observables/BaseObservable.h`).
- Return references to objects owned by the engine (particles, their neighbours, molecules) with `py::return_value_policy::reference` (as `ConfigInfo.particles` and `BaseParticle.n3` do), so that Python never deletes them; never transfer the ownership of engine objects to Python.
- Every bound class, method, property and argument needs a docstring in the existing format: a raw `R"pbdoc(...)pbdoc"` string in numpydoc style, with `Parameters` and `Returns` sections and Sphinx cross-references (`:class:`, `:meth:`, `:attr:`). Name the arguments with `py::arg(...)` so that they can be passed by keyword and appear in the documentation.
- New classes must also be listed in the documentation, in the `autosummary`/`autoclass` blocks of the relevant page in `docs/source/oxpy/modules/` (`core/core.md`, `core/forces.md`, `core/observables.md`, `analysis.md`). If the feature changes how users drive simulations, update the guides in `docs/source/oxpy/` too (e.g. `extending/forces.md` and `extending/observables.md`).
- Keep the Python API backwards-compatible: do not rename or remove existing classes, methods, properties or arguments. Deprecate them first.
- Engine errors must reach Python as exceptions: throw `oxDNAException` in the C++ code (it is translated to `oxpy.core.OxDNAError`); never call `exit()` or `abort()` in code that oxpy can reach.
- Exposed C++ code must not rely on state that only exists in the command-line executables, and must leave the global state (`CONFIG_INFO`, logger, singletons) consistent so that several `Context` blocks can run in one Python process. Multiple consecutive contexts are known to segfault unless the manager is deleted at the end of each block (`del manager`, see the docs); don't make this worse, and say so in the PR if a change touches the context or manager lifetime.

## Tests

- Add a test for every new binding, either as an oxpy-level test folder in `test/` (an `oxpy_input.py` script and an `oxpy_compare` file, listed in `test/test_folder_list.txt`; see `test/OXPY`) or in `analysis/tests/` if the binding is used by oat.
- Test the behaviour through the Python API (the values a method returns, the effect of a setter on the simulation), not only that the call does not raise.
- Run `examples/OXPY/run.py` after any change to `OxpyContext`, `OxpyManager` or the module definition.
