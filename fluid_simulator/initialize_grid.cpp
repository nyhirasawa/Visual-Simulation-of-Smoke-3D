#include "initialize_grid.h"
#include "physical_const.h"
#include "utils.h"

#include <math.h>


namespace smoke_simulation{
//グリッドを初期化する関数
void initialize_grid(Grid& all_grid){

    //velocity_in_voxel_face_xの初期化
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x+1;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix, iy, iz)]=0.0;
            }
        }
    }

    //velocity_in_voxel_face_yの初期化
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y+1;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                if(ix>(6.0*(smoke_simulation::physical_const::kGrid_num_x+1.0)/15.0)
                 &&ix<(9.0*(smoke_simulation::physical_const::kGrid_num_x+1.0)/15.0)
                 &&iy>(6.0*smoke_simulation::physical_const::kGrid_num_y/15.0)
                 &&iy<(9.0*smoke_simulation::physical_const::kGrid_num_y/15.0)
                 &&iz>(6.0*smoke_simulation::physical_const::kGrid_num_z/15.0)
                 &&iz<(9.0*smoke_simulation::physical_const::kGrid_num_z/15.0)){
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz)]=-10.0;
                }
                else{
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, iy, iz)]=0.0;
                }
            }
        }
    }

    //velocity_in_voxel_face_zの初期化
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z+1;iz++){
                all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, iz)]=0.0;
            }
        }
    }

    //voxel_centerの初期化
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                if(ix>(6.0*smoke_simulation::physical_const::kGrid_num_x/15.0)
                 &&ix<(9.0*smoke_simulation::physical_const::kGrid_num_x/15.0)
                 &&iy>(6.0*smoke_simulation::physical_const::kGrid_num_y/15.0)
                 &&iy<(9.0*smoke_simulation::physical_const::kGrid_num_y/15.0)
                 &&iz>(6.0*smoke_simulation::physical_const::kGrid_num_z/15.0)
                 &&iz<(9.0*smoke_simulation::physical_const::kGrid_num_z/15.0)){
                    all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz)]=5.0;
                }
                else{
                    all_grid.substance_density[get_voxel_center_index_3D(ix, iy, iz)]=0.0;
                }
            }
        }
    }

}
}//namespace smoke_simulation
