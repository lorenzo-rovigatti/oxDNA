/*
 * ConfigInfo.h
 *
 *  Created on: Mar 21, 2020
 *      Author: lorenzo
 */

#ifndef OXPY_BINDINGS_INCLUDES_FLATTENEDCONFIGINFO_H_
#define OXPY_BINDINGS_INCLUDES_FLATTENEDCONFIGINFO_H_

#include "../python_defs.h"

#include <Utilities/FlattenedConfigInfo.h>

void export_FlattenedConfigInfo(py::module &m) {
	py::class_<FlattenedVectorArray> vector_array(m, "FlattenedVectorArray", py::buffer_protocol(), R"pbdoc(
        This is a view of a Nx3 matrix (where N is the number of particles in the simulation).

        Instances can be transparently converted to numpy arrays::

            # this assumes that va is a FlattenedVectorArray
            np_va = np.array(va)             # this will copy the data over to Python
            np_va = np.array(va, copy=False) # this will *not* copy the data

    )pbdoc");

	vector_array.def_buffer([](FlattenedVectorArray &va) -> py::buffer_info {
		return py::buffer_info(
				va.data.data(), 							/* Pointer to buffer */
				sizeof(number), 							/* Size of one scalar */
				py::format_descriptor<number>::format(), 	/* Python struct-style format descriptor */
				2, 											/* Number of dimensions */
				{	va.rows(), va.cols()}, 					/* Buffer dimensions */
				{	sizeof(number) * va.cols(), 			/* Strides (in bytes) for each index */
					sizeof(number)}
		);
	});

	py::class_<FlattenedConfigInfo, std::shared_ptr<FlattenedConfigInfo>> conf_info(m, "FlattenedConfigInfo", R"pbdoc(
         An object storing particle types, positions, a1 and a3 orientation vectors in flat arrays that can be used as or converted to numpy arrays.

         Note that the data stored in this object is updated maximum once per time step. The update is performed when the object
         is accessed from the :py:class:`ConfigInfo` object, meaning that it is better not to store pointers to this object if you
         do not want to risk accessing outdated data. The vector types are stored in :py:class:`FlattenedVectorArray` objects, which 
         can be used like this::

             flat_positions = manager.config_info().flattened_conf.positions
             np_poss = np.array(flat_positions, copy=False)
             print(np.average(flat_positions, axis=0) == np.average(np_poss, axis=0))

	)pbdoc");

	conf_info.def_readonly("positions", &FlattenedConfigInfo::positions, R"pbdoc(Particle positions)pbdoc");
	conf_info.def_readonly("a1s", &FlattenedConfigInfo::a1s, R"pbdoc(Particle a1 orientation vectors)pbdoc");
	conf_info.def_readonly("a3s", &FlattenedConfigInfo::a3s, R"pbdoc(Particle a3 orientation vectors)pbdoc");
	conf_info.def_readonly("types", &FlattenedConfigInfo::types, R"pbdoc(Particle types)pbdoc");
}

#endif /* OXPY_BINDINGS_INCLUDES_FLATTENEDCONFIGINFO_H_ */
