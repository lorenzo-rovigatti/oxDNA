* The `number` and `number4` templated stuff has been removed from the code
* On the CPU side, `number` is now set at compile time to either float or double
* On the CUDA side, `number` has become `c_number` and `number4` has become `c_number4`. Both are set at compile time to `float` and `float4` or `double` and `LR_double4`. If the CUDA backend is compiled in `float` precision (*i.e.* if `-DCUDA_DOUBLE=OFF`, default behaviour), the `mixed` precision is available and can be enabled by putting in the input file `backend_precision=mixed`. The `backend_precision` option is otherwise not used any more.
* For plugins: before this change, one had to write two entry points, one for each precision (`float` and `double`). Now only one is required, which should be called `make_MyPlugin`, where `MyPlugin` is the name of the plugin class (and file).
* `BaseParticle`: `int_centers` is now a `std::vector` and hence it does not have to be allocated but `resize`'d. The number of interaction centers can be accessed with the `N_int_centers()` helper method.
* SimBackend's `_N` and ConfigInfo's `N` members have been removed. Both classes now possess a `N()` method that directly returns the number of particles.
* ConfigInfo's `particles` is now a method that returns a reference to the particles' `std::vector`
* ConfigInfo's `instance` returns a `std::shared_ptr` rather than a bare pointer
* The `RELEASE` macro defined in `src/defs.h` **must** match the current git tag or it won't be possible to commit changes to the repository. This is to make sure that the version printed by the executables is (at least decently) up to date.
