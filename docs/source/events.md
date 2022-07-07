# Events

oxDNA supports a basic [observer pattern](https://en.wikipedia.org/wiki/Observer_pattern) to associate callbacks to specific events. The system can be used from both C++ and Python. The basic idea is to call `ConfigInfo`'s `subscribe` and `notify` methods to register callbacks and to trigger events, respectively.

````{admonition} Python
```python
import oxpy

[...]

# define a callback
def on_T_update():
  print("Temperature was updated to", conf_info.temperature())

# associate the callback to the "T_updated" event
conf_info = manager.config_info()
conf_info.subscribe("T_updated", on_T_update)

[...]

# somewhere else we fire off the event, which will trigger 
# the invoking of all the associated callbacks
conf_info.notify("T_updated")

```
```` 

````{admonition} C++
```c++
#include "../Utilities/ConfigInfo.h"

[...]

// associate a lambda function that calls a class method to the 
// "box_updated" event
CONFIG_INFO->subscribe("box_updated", [this]() { this->_on_box_update(); });

[...]

// somewhere else we fire off the event, which will trigger 
// the invoking of all the associated callbacks
CONFIG_INFO->notify("box_updated");
```
````

## List of supported events

* `T_updated`: triggered when the simulation temperature is changed (for instance by calling {meth}`~oxpy.core.OxpyManager.update_temperature`).
* `box_updated`: triggered when the simulation box is changed. When this event is fired off the simulation is **always** in a valid state.
* `box_initialised`: triggered when the simulation box is re-initialised. Note that this event is fired off even if the box is re-initialised with the same values or during a trial Monte Carlo volume move. In this latter case the move may be reverted, and therefore you cannot assume that every time the event is triggered the simulation is in a valid state.
