# Events

oxDNA supports a basic [observer pattern](https://en.wikipedia.org/wiki/Observer_pattern) to associate callbacks and specific events. The system can be used from both C++ and Python.

## List of supported events

* `T_updated`: triggered when the simulation temperature is changed (for instance by calling {meth}`~oxpy.core.OxpyManager.update_temperature`).
* `box_updated`: triggered when the simulation box is changed. When this event is fired off the simulation is **always** in a valid state.
* `box_initialised`: triggered when the simulation box is re-initialised. Note that this event is fired off even if the box is re-initialised with the same values or during a trial Monte Carlo volume move. In this latter case the move may be reverted, and therefore you cannot assume that every time the event is triggered the simulation is in a valid state.
