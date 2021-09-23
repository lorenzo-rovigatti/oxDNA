* The `number` and `number4` templated stuff has been removed from the code.
* As a result, most (if not all) of the `this->` occurrences are not necessary any more.
* On the CPU side, `number` is now set at compile time to either float or double.
* On the CUDA side, `number` has become `c_number` and `number4` has become `c_number4`. Both are set at compile time to `float` and `float4` or `double` and `LR_double4`. If the CUDA backend is compiled in `float` precision (*i.e.* if `-DCUDA_DOUBLE=OFF`, default behaviour), the `mixed` precision is available and can be enabled by putting in the input file `backend_precision=mixed`. The `backend_precision` option is otherwise not used any more.
* For plugins: before this change, one had to write two entry points, one for each precision (`float` and `double`). Now only one is required, which should be called `make_MyPlugin`, where `MyPlugin` is the name of the plugin class (and file).
* `BaseParticle`: `int_centers` is now a `std::vector` and hence it does not have to be allocated but `resize`'d. The number of interaction centers can be accessed with the `N_int_centers()` helper method.
* SimBackend's `_N` and ConfigInfo's `N` members have been removed. Both classes now possess a `N()` method that directly returns the number of particles.

* ConfigInfo's `particles` is now a method that returns a reference to the particles' `std::vector`.
* ConfigInfo's `instance` returns a `std::shared_ptr` rather than a bare pointer.
* The `RELEASE` macro defined in `src/defs.h` **must** match the current git tag or it won't be possible to commit changes to the repository. This is to make sure that the version printed by the executables is (at least decently) up to date.
* The functions that used to initialise, clean and print the `input_file` structure have now been included as methods in the class itself. 
* The `init` methods of the `ObservableOutput` and `BaseObservable` classes do not require any parameter. This is in contrast with previous versions in which both methods required a `ConfigInfo` object. 

## Interactions

* The signature of all the `pair_interaction*` methods have been changed. The parameter list went from `(BaseParticle *, BaseParticle *, LR_vector *, bool)` to `(BaseParticle *, BaseParticle *, bool, bool)`, where the third parameter (`compute_r`) tells the method whether `r` should be computed from scratch (the default). If not, the method should use the new `_computed_r` member, which should be set from the methods called with `compute_r==true`, or can be set externally through the `set_computed_r()` method. **Nota Bene:** many calls to the `pair_interaction_*` methods used `NULL` as a parameter. Since the latter is (on most systems) equal to zero, the compiler will happily cast it to `false` without complaining. This might give (and have given) raise to hard-to-find bugs.
* The need for the curiously-recurring template pattern has been replaced by a vector of std::function. The `ADD_INTERACTION_TO_MAP(index, method)` macro can be used to fill the interaction map with a method of the same class, however other classes' methods can also be added by using a more complicated (lambda-based) syntax. Thanks for this change there is no need for each interaction to also redefine the `pair_interaction_term` method. As a result, the `_pair_interaction_term_wrapper` method has been removed. This change breaks compilation of older interactions.
