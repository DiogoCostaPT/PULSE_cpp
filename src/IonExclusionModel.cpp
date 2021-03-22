// Copyright 2021: Diogo Costa

// This program, PULSE_cpp, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "IonExclusionModel.h"

/* *****
 * Ion Exclusion model : Exchange with immobile phase (just exchange)
 * **** */

void IonExclusionModel(globalpar& gp,globalvar& gv,double* deltt)
{

    int il, ih;

    (*gv.exchange_is)  = (*deltt) * gp.alphaIE * (*gv.c_s) % (*gv.v_swe) ; 
    float exchange_is_l = 0.0f;

    // limit the flux to the available material
    for(il=0;il<gv.nl ;il++){
        for(ih=0;ih<gv.wetfront_cell ;ih++){
            if ((*gv.exchange_is).at(il,ih) > 0 && (*gv.v_liq).at(il,ih) > 0.0001){

                //std::cout << std::to_string((*gv.v_liq).at(il,ih)) << std::endl;
                exchange_is_l = fmin((*gv.exchange_is).at(il,ih),
                                    (*gv.c_s).at(il,ih) * (*gv.v_swe).at(il,ih));   

                (*gv.c_s).at(il,ih) = ((*gv.c_s).at(il,ih) * (*gv.v_swe).at(il,ih) 
                        - exchange_is_l) / (*gv.v_swe).at(il,ih);

                (*gv.c_m).at(il,ih) = ((*gv.c_m).at(il,ih) * (*gv.v_liq).at(il,ih) 
                        + exchange_is_l) / (*gv.v_liq).at(il,ih);

                (*gv.c_s).at(il,ih) = std::fmax((*gv.c_s).at(il,ih),0.0f);
                (*gv.c_m).at(il,ih) = std::fmax((*gv.c_m).at(il,ih),0.0f);  

            }
        }
    }
    
}