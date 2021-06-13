#include "write_substance_density_data.h"
#include <fstream>
#include <string>
#include <sstream>
#include "physical_const.h"
#include "utils.h"


namespace smoke_simulation{

void write_substance_density_data(Grid& all_grid, int file_number){
    std::ostringstream filename;
    filename<<"../substance_density_data/substance_density_"<<file_number<<".dat"<<std::flush;
//    std::string filename = "test.txt";
    std::ofstream writing_file;
    writing_file.open(filename.str(), std::ios::out);
    writing_file<<smoke_simulation::physical_const::kGrid_num_x<<" "
                <<smoke_simulation::physical_const::kGrid_num_y<<" "
                <<smoke_simulation::physical_const::kGrid_num_z<<std::endl;
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                writing_file<<all_grid.substance_density[get_voxel_center_index_3D(ix,iy,iz)]<<std::endl;
            }
        }
    }
}

}
