#include "DNA3Interaction.h"

#include "../Particles/DNANucleotide.h"
#include <iostream>

DNA3Interaction::DNA3Interaction() :
                                DNA2Interaction() {
                std::cout << "\n \n \n  using DNA3Interaction \n \n \n" ;
                _grooving = 2; 
}

DNA3Interaction::~DNA3Interaction() {

}
