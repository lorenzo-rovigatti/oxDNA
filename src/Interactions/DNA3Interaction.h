#include "DNA2Interaction.h"


class DNA3Interaction: public DNA2Interaction {

protected:

public:
	DNA3Interaction(){
		std::cout << "\n \n \n  using DNA3Interaction \n \n \n" ;
		_grooving = 2; 
}
};

/**
 * @brief Handles interactions between DNA nucleotides without using meshes.
 *
 * This interaction is selected with
 * interaction_type = DNA2_nomesh
 */

class DNA3Interaction_nomesh : public DNA3Interaction {
protected:


public:

};

