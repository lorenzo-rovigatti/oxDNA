# Modifying external forces at runtime

While it is not yet possible to write external forces entirely in Python, with oxpy you can change the parameters of some external forces while the simulation is running. The list of the available forces (and the parameters that can be changed) can be found [here](../modules/core/forces.md). All external forces can be accessed by using {class}`~oxpy.core.ConfigInfo`'s {attr}`~oxpy.core.ConfigInfo.forces` attribute.

Here is an example that can be run in the `examples/TRAPS` folder which doubles the stiffness of all forces every 20000 simulation steps:

```python
import oxpy

with oxpy.Context():
    manager = oxpy.OxpyManager("inputMD")

    for i in range(10):
        manager.run(20000)

        # double the stiffness of all forces after every iteration
        for force in manager.config_info().forces:
            force.stiff *= 2
```

Note that you can modify only some specific forces by using the `id` and `group_name` options in the external forces file.

For instance, imagine we have the following forces defined in the external forces file (similar to the `external.conf` file 
of the `examples/TRAPS` folder):

```text
{
type = trap
id = first_force
group_name = pulling_forces
particle = 0
stiff = 1.00
pos0 = 46.3113780977, 11.1604626391, 26.8730311801
rate = 0.
dir = 0.,0.,-1.
}

{
type = trap
id = last_force
group_name = pulling_forces
particle = 99
pos0 = 83.1532046751, 15.950789638, 37.3071701142
stiff = 1.0
rate = 0.
dir = 0.,0.,1.
}
```

We can now access only these two forces by using {class}`~oxpy.core.forces.BaseForce`'s {attr}`~oxpy.core.forces.BaseForce.group_name` attribute:

```python
for force in manager.config_info().forces:
	if force.group_name == "pulling_forces":
		force.stiff *= 2
```

or we can access one specific force by using {class}`~oxpy.core.ConfigInfo`'s {attr}`~oxpy.core.ConfigInfo.get_force_by_id`:

```python
my_force = manager.config_info().get_force_by_id("first_force")
my_force.stiff *= 2
```

```{note}
The {attr}`~oxpy.core.forces.BaseForce.id` and {attr}`~oxpy.core.forces.BaseForce.group_name` attributes can also be changed from the Python's side!
```
