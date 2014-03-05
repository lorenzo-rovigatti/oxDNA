/*
 * DNA interaction class that does not use meshes
 * types
 */

#ifndef DNA_INTERACTION_NOMESH_H
#define DNA_INTERACTION_NOMESH_H

#include "DNAInteraction.h"

/**
 * @brief Handles interactions between DNA nucleotides without using meshes.
 *
 * This interaction is selected with
 * interaction_type = DNA_nomesh
 */
template <typename number>
class DNAInteraction_nomesh : public DNAInteraction<number> {
protected:

	/**
	 * @brief Custom function that returns f4.
	 *
	 * @param cost  argument of f4
	 * @param i     type of the interaction
	 */
	//virtual number _custom_f4 (number cost, int i) { return this->_query_mesh (cost, this->_mesh_f4[i]); }
	virtual number _custom_f4 (number cost, int i) {
		if (i != CXST_F4_THETA1) return this->_fakef4 (cost, (void *)&i);
		else return this->_fakef4_cxst_t1 (cost, (void *)&i);
	}

	/**
	 * @brief Custom function that returns the derivative of f4. See _custom_f4
	 *
	 * @param cost  argument of f4D
	 * @param i     type of the interaction
	 */
	//virtual number _custom_f4D (number cost, int i) { return this->_query_meshD (cost, this->_mesh_f4[i]); }
	virtual number _custom_f4D (number cost, int i) {
		if (i != CXST_F4_THETA1) return this->_fakef4D (cost, (void *)&i);
		else return this->_fakef4D_cxst_t1 (cost, (void *)&i);
	}

public:
	DNAInteraction_nomesh();
	virtual ~DNAInteraction_nomesh();

};

#endif /* DNA_INTERACTION_NOMESH_H */

