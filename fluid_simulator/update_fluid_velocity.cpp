#include "update_fluid_velocity.h"

#include <math.h>
#include<iostream>
#include <chrono>//時間計測用
#include <fstream>//ファイル書き出し用

#include "linear_solver.h"
#include "physical_const.h"
#include "sparse_matrix.h"
#include "utils.h"

namespace smoke_simulation{

void add_source(Grid& all_grid){
/*
    for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;++iy){
        if(iy>(7*smoke_simulation::physical_const::kGrid_num_y/15)
         &&iy<(8*smoke_simulation::physical_const::kGrid_num_y/15)){
            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0, iy)]=10.0;
            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(0, iy)]=0.0;
            all_grid.substance_density[get_voxel_center_index(0, iy)]=5.0;
//            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(1, iy)]=10.0;
//            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(1, iy)]=0.0;
//            all_grid.substance_density[get_voxel_center_index(1, iy)]=5.0;
        }
    }
*/
}

//壁における速度場をセットする関数
void set_boundary_velocity(Grid& all_grid){
    //y-z平面に平行な壁の速度(つまり速度のx成分)をセット
    for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
        for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0,iy,iz)]
                =-all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1,iy,iz)];
            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, iy, iz)]
                =-all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, iz)];
        }
    }
    for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, iy, 0)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1, iy, 0)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, iy, 1)])/2.0;
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, iy, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x,iy,0)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1,iy,0)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x,iy,1)])/2.0;
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, iy, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
    }
    for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0,0,iz)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1,0,iz)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0,1,iz)])/2.0;
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, iz)])/2.0;
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x,0,iz)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1,0,iz)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x,1,iz)])/2.0;
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
            =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
             +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-2, iz)])/2.0;
    }
    //x-z平面に平行な壁の速度(つまり速度のy成分)をセット
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 0, iz)]
                =-all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 1, iz)];
            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y, iz)]
                =-all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, iz)];
        }
    }
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 0, 0)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 1, 0)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 0, 1)])/2.0;
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y, 0)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y, 1)])/2.0;
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
    }
    for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 0, iz)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 1, iz)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(1, 0, iz)])/2.0;
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, iz)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 1, iz)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-2, 0, iz)])/2.0;
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y, iz)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(1, smoke_simulation::physical_const::kGrid_num_y, iz)])/2.0;
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y, iz)]
            =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
             +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y, iz)])/2.0;
    }
    //x-y平面に平行な壁の速度(つまり速度のz成分)をセット
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
            all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, 0)]
                =-all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, 1)];
            all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, smoke_simulation::physical_const::kGrid_num_z)]
                =-all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, iy, smoke_simulation::physical_const::kGrid_num_z-1)];
        }
    }
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, 0, 0)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, 0, 1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, 1, 0)])/2.0;
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, 1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, smoke_simulation::physical_const::kGrid_num_y-2, 0)])/2.0;
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, 0,smoke_simulation::physical_const::kGrid_num_z)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, 1, smoke_simulation::physical_const::kGrid_num_z)])/2.0;
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z)])/2.0;
    }
    for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, iy, 0)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, iy, 1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(1, iy, 0)])/2.0;
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, 0)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, 1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-2, iy, 0)])/2.0;
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, iy, smoke_simulation::physical_const::kGrid_num_z)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(1, iy, smoke_simulation::physical_const::kGrid_num_z)])/2.0;
        all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, smoke_simulation::physical_const::kGrid_num_z)]
            =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-2, iy, smoke_simulation::physical_const::kGrid_num_z)])/2.0;
    }

    //8個の隅の速度場は周辺から線形補間
    //x成分の速度
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, 0, 0)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1, 0, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, 1, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, 0, 1)])/3.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, 1)])/3.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, 0, 0)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, 1, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, 0, 1)])/3.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-2, 0)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-1, 1)])/3.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    //y成分の速度
    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 0, 0)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(1, 0, 0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 1, 0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 0, 1)])/3.0;
    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(1, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, 0)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-2, 0, 0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 1, 0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, 1)])/3.0;
    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-2, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
     all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y, 0)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(1, smoke_simulation::physical_const::kGrid_num_y, 0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y, 1)])/3.0;
     all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(1, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(0, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
     all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y, 0)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y,0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1,0)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y,1)])/3.0;
     all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    //z成分の速度
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, 0, 0)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(1, 0, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, 1, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, 0, 1)])/3.0;
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, 0)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-2, 0, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, 1, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, 1)])/3.0;
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, 1)])/3.0;
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-2, 0)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, 1)])/3.0;
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(1, 0, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, 1, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z-1)])/3.0;
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, smoke_simulation::physical_const::kGrid_num_z)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-2, 0, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, 1, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, smoke_simulation::physical_const::kGrid_num_z-1)])/3.0;
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)])/3.0;
    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z)]
        =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z)]
         +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)])/3.0;
}

//壁における圧力場をセットする関数
void set_boundary_pressure(Grid& all_grid){
    //y-z平面に平行な壁の速度(つまり速度のx成分)をセット
    for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
        for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
            all_grid.pressure[get_voxel_center_index_3D(0,iy,iz)]
                =all_grid.pressure[get_voxel_center_index_3D(1,iy,iz)];
            all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, iz)]
                =all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2, iy, iz)];
        }
    }
    //x-z平面に平行な壁の速度(つまり速度のy成分)をセット
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
            all_grid.pressure[get_voxel_center_index_3D(ix, 0, iz)]
                =all_grid.pressure[get_voxel_center_index_3D(ix, 1, iz)];
            all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
                =all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-2, iz)];
        }
    }
    //x-y平面に平行な壁の速度(つまり速度のz成分)をセット
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
            all_grid.pressure[get_voxel_center_index_3D(ix, iy, 0)]
                =all_grid.pressure[get_voxel_center_index_3D(ix, iy, 1)];
            all_grid.pressure[get_voxel_center_index_3D(ix, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
                =all_grid.pressure[get_voxel_center_index_3D(ix, iy, smoke_simulation::physical_const::kGrid_num_z-2)];
        }
    }
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        all_grid.pressure[get_voxel_center_index_3D(ix, 0, 0)]
            =(all_grid.pressure[get_voxel_center_index_3D(ix, 1, 0)]
             +all_grid.pressure[get_voxel_center_index_3D(ix, 0, 1)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(ix, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.pressure[get_voxel_center_index_3D(ix, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.pressure[get_voxel_center_index_3D(ix, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
            =(all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-2, 0)]
             +all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, 1)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.pressure[get_voxel_center_index_3D(ix, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
    }
    for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
        all_grid.pressure[get_voxel_center_index_3D(0, iy, 0)]
            =(all_grid.pressure[get_voxel_center_index_3D(1, iy, 0)]
             +all_grid.pressure[get_voxel_center_index_3D(0, iy, 1)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(0, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.pressure[get_voxel_center_index_3D(1, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.pressure[get_voxel_center_index_3D(0, iy, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1,iy,0)]
            =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2,iy,0)]
             +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1,iy,1)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
            =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2, iy, smoke_simulation::physical_const::kGrid_num_z-1)]
             +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, iy, smoke_simulation::physical_const::kGrid_num_z-2)])/2.0;
    }
    for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
        all_grid.pressure[get_voxel_center_index_3D(0,0,iz)]
            =(all_grid.pressure[get_voxel_center_index_3D(1,0,iz)]
             +all_grid.pressure[get_voxel_center_index_3D(0,1,iz)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
            =(all_grid.pressure[get_voxel_center_index_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
             +all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, iz)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1,0,iz)]
            =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2,0,iz)]
             +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1,1,iz)])/2.0;
        all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
            =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y-1, iz)]
             +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-2, iz)])/2.0;
    }
    //8個の隅の速度場は周辺から線形補間
    //x成分の速度
    all_grid.pressure[get_voxel_center_index_3D(0, 0, 0)]
        =(all_grid.pressure[get_voxel_center_index_3D(1, 0, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(0, 1, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(0, 0, 1)])/3.0;
    all_grid.pressure[get_voxel_center_index_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.pressure[get_voxel_center_index_3D(1, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(0, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(0, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
        =(all_grid.pressure[get_voxel_center_index_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, 1)])/3.0;
    all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.pressure[get_voxel_center_index_3D(1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(0, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, 0)]
        =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2, 0, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, 1, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, 1)])/3.0;
    all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2, 0, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, 1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, 0, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
    all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
        =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y-1, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-2, 0)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, 1)])/3.0;
    all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
        =(all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-2, smoke_simulation::physical_const::kGrid_num_z-1)]
         +all_grid.pressure[get_voxel_center_index_3D(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1, smoke_simulation::physical_const::kGrid_num_z-2)])/3.0;
}

//(vorticityの計算に使う)cell centered velocityの計算
void calc_cell_centered_velocity(Grid& all_grid){
    //cell_centered_velocityの計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy,iz)][0]
                    =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]+all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)])/2.0;
                all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy,iz)][1]
                    =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]+all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)])/2.0;
                all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy,iz)][2]
                    =(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]+all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)])/2.0;
            }
        }
    }
}

//vorticity confinement termの計算
void calc_vorticity_confinement(Grid& all_grid){
    //vorticityの計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                //範囲外を参照しないようにする処置
                if(ix==0||ix==smoke_simulation::physical_const::kGrid_num_x-1
                 ||iy==0||iy==smoke_simulation::physical_const::kGrid_num_y-1
                 ||iz==0||iz==smoke_simulation::physical_const::kGrid_num_z-1){
                    all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][0]=0.0;
                    all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][1]=0.0;
                    all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][2]=0.0;
                }
                else{
                    //3次元バージョン
                    all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][0]
                        =(all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy+1,iz)][2]
                         -all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy-1,iz)][2]
                         -all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy,iz+1)][1]
                         +all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy,iz-1)][1])
                         /(2.0*smoke_simulation::physical_const::kCell_length);
                    all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][1]
                        =(all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy,iz+1)][0]
                         -all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy,iz-1)][0]
                         -all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix+1,iy,iz)][2]
                         +all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix-1,iy,iz)][2])
                         /(2.0*smoke_simulation::physical_const::kCell_length);
                    all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][2]
                        =(all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix+1,iy,iz)][1]
                         -all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix-1,iy,iz)][1]
                         -all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy+1,iz)][0]
                         +all_grid.velocity_in_cell_center[get_voxel_center_index_3D(ix,iy-1,iz)][0])
                         /(2.0*smoke_simulation::physical_const::kCell_length);
                }
            }
        }
    }
    //ノルムを1に正規化したvorticity amplitude gradient(論文のN)の計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                //範囲外を参照しないようにする処置
                if(ix==0||ix==smoke_simulation::physical_const::kGrid_num_x-1
                 ||iy==0||iy==smoke_simulation::physical_const::kGrid_num_y-1
                 ||iz==0||iz==smoke_simulation::physical_const::kGrid_num_z-1){
                     all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][0]=0.0;
                     all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][1]=0.0;
                     all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][2]=0.0;
                }
                else{
                    double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index_3D(ix-1,iy,iz)][0]*all_grid.vorticity[get_voxel_center_index_3D(ix-1,iy,iz)][0]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix-1,iy,iz)][1]*all_grid.vorticity[get_voxel_center_index_3D(ix-1,iy,iz)][1]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix-1,iy,iz)][2]*all_grid.vorticity[get_voxel_center_index_3D(ix-1,iy,iz)][2]);
                    double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index_3D(ix+1,iy,iz)][0]*all_grid.vorticity[get_voxel_center_index_3D(ix+1,iy,iz)][0]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix+1,iy,iz)][1]*all_grid.vorticity[get_voxel_center_index_3D(ix+1,iy,iz)][1]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix+1,iy,iz)][2]*all_grid.vorticity[get_voxel_center_index_3D(ix+1,iy,iz)][2]);
                    double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index_3D(ix,iy-1,iz)][0]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy-1,iz)][0]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy-1,iz)][1]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy-1,iz)][1]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy-1,iz)][2]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy-1,iz)][2]);
                    double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index_3D(ix,iy+1,iz)][0]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy+1,iz)][0]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy+1,iz)][1]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy+1,iz)][1]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy+1,iz)][2]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy+1,iz)][2]);
                    double vorticity_amplitude_z0=sqrt(all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz-1)][0]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz-1)][0]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz-1)][1]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz-1)][1]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz-1)][2]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz-1)][2]);
                    double vorticity_amplitude_z1=sqrt(all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz+1)][0]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz+1)][0]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz+1)][1]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz+1)][1]
                                                      +all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz+1)][2]*all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz+1)][2]);
                    all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/(2.0*smoke_simulation::physical_const::kCell_length);
                    all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/(2.0*smoke_simulation::physical_const::kCell_length);
                    all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][2]=(vorticity_amplitude_z1-vorticity_amplitude_z0)/(2.0*smoke_simulation::physical_const::kCell_length);
                }
                double normalize_factor=sqrt(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][0]*all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][0]
                                            +all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][1]*all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][1]
                                            +all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][2]*all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][2]);
                //0除算を回避するための処理
                if(normalize_factor>0.0001){
                    all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][0]/=normalize_factor;
                    all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][1]/=normalize_factor;
                    all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix, iy, iz)][2]/=normalize_factor;
                }
            }
        }
    }
    //vorticity confinment term の計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][0]
                    +=smoke_simulation::physical_const::kConfinement_amplitude
                     *smoke_simulation::physical_const::kCell_length
                     *(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][1]
                      *all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][2]
                      -all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][2]
                      *all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][1]);
                all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][1]
                    +=smoke_simulation::physical_const::kConfinement_amplitude
                     *smoke_simulation::physical_const::kCell_length
                     *(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][2]
                      *all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][0]
                      -all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][0]
                      *all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][2]);
                all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][2]
                    +=smoke_simulation::physical_const::kConfinement_amplitude
                     *smoke_simulation::physical_const::kCell_length
                     *(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][0]
                      *all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][1]
                      -all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index_3D(ix,iy,iz)][1]
                      *all_grid.vorticity[get_voxel_center_index_3D(ix,iy,iz)][0]);
            }
        }
    }
}

//外力場による速度場の更新
void update_fluid_velocity_by_external_force(Grid& all_grid){
    //計算した外力場により速度場を更新
    //x成分の更新
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x+1;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                if(ix-1<0){
                    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                        +=all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][0]
                         *smoke_simulation::physical_const::kDt;
                }
                else if(ix>=smoke_simulation::physical_const::kGrid_num_x){
                    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                        +=all_grid.external_force_field[get_voxel_center_index_3D(ix-1,iy,iz)][0]
                         *smoke_simulation::physical_const::kDt;
                }
                else{
                    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                        +=((all_grid.external_force_field[get_voxel_center_index_3D(ix-1,iy,iz)][0]
                         +all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][0])/2.0)
                         *smoke_simulation::physical_const::kDt;
                }
            }
        }
    }
    //y成分の更新
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y+1;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                if(iy-1<0){
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                        +=all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][1]
                         *smoke_simulation::physical_const::kDt;
                }
                else if(iy>=smoke_simulation::physical_const::kGrid_num_y){
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                        +=all_grid.external_force_field[get_voxel_center_index_3D(ix,iy-1,iz)][1]
                         *smoke_simulation::physical_const::kDt;
                }
                else{
                    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                        +=((all_grid.external_force_field[get_voxel_center_index_3D(ix,iy-1,iz)][1]
                           +all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][1])/2.0)
                          *smoke_simulation::physical_const::kDt;
                }
            }
        }
    }
    //z成分の更新
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z+1;iz++){
                if(iz-1<0){
                    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]
                        +=all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][2]
                         *smoke_simulation::physical_const::kDt;
                }
                else if(iz>=smoke_simulation::physical_const::kGrid_num_z){
                    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]
                        +=all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz-1)][2]
                         *smoke_simulation::physical_const::kDt;
                }
                else{
                    all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]
                        +=((all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz-1)][2]
                           +all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][2])/2.0)
                          *smoke_simulation::physical_const::kDt;
                }
            }
        }
    }

}


//外力項の計算
void add_force_fluid(Grid& all_grid){
    //外力場をリセット
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][0]=0.0;
                all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][1]=0.0;
                all_grid.external_force_field[get_voxel_center_index_3D(ix,iy,iz)][2]=0.0;
            }
        }
    }
    calc_cell_centered_velocity(all_grid);
/*
    //重力の効果
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][1]-=9.8*all_grid.substance_density[get_voxel_center_index(ix, iy)]*smoke_simulation::physical_const::kDt;
        }
    }
*/
    calc_vorticity_confinement(all_grid);
    update_fluid_velocity_by_external_force(all_grid);
}

//x成分のvelocityのlinear interpolation
double linear_interpolation_x_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    //バックトレース先の座標のindex
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(double)(advected_index_z);
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-1){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-2;
        advected_z=(double)(advected_index_z);
    }
    //バックトレース先の速度を線形補間する
    double a0, a1;
    double b0, b1;
    double c0, c1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
    c0=advected_z-advected_index_z;
    c1=1.0-c0;
    return a1*b1*c1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x,advected_index_y,advected_index_z)]
          +a1*b1*c0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x,advected_index_y,advected_index_z+1)]
          +a1*b0*c1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x,advected_index_y+1,advected_index_z)]
          +a1*b0*c0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x,advected_index_y+1,advected_index_z+1)]
          +a0*b1*c1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x+1,advected_index_y,advected_index_z)]
          +a0*b1*c0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x+1,advected_index_y,advected_index_z+1)]
          +a0*b0*c1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x+1,advected_index_y+1,advected_index_z)]
          +a0*b0*c0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x+1,advected_index_y+1,advected_index_z+1)];
}

//y成分のvelocityのlinear interpolation
double linear_interpolation_y_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    //バックトレース先の座標のindex
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)(advected_index_x);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(double)(advected_index_x);
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)advected_index_y;
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-1;
        advected_y=(double)advected_index_y;
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(double)advected_index_z;
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-1){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-2;
        advected_z=(double)advected_index_z;
    }
    //バックトレース先の速度を線形補間する
    double a0, a1;
    double b0, b1;
    double c0, c1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
    c0=advected_z-advected_index_z;
    c1=1.0-c0;
    return a1*b1*c1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x,advected_index_y,advected_index_z)]
          +a1*b1*c0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x,advected_index_y,advected_index_z+1)]
          +a1*b0*c1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x,advected_index_y+1,advected_index_z)]
          +a1*b0*c0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x,advected_index_y+1,advected_index_z+1)]
          +a0*b1*c1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x+1,advected_index_y,advected_index_z)]
          +a0*b1*c0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x+1,advected_index_y,advected_index_z+1)]
          +a0*b0*c1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x+1,advected_index_y+1,advected_index_z)]
          +a0*b0*c0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x+1,advected_index_y+1,advected_index_z+1)];
}

//z成分のvelocityのlinear interpolation
double linear_interpolation_z_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    //バックトレース先の座標のindex
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)(advected_index_x);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(double)(advected_index_x);
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)advected_index_y;
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(double)advected_index_y;
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(double)advected_index_z;
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-1;
        advected_z=(double)advected_index_z;
    }
    //バックトレース先の速度を線形補間する
    double a0, a1;
    double b0, b1;
    double c0, c1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
    c0=advected_z-advected_index_z;
    c1=1.0-c0;
    return a1*b1*c1*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x,advected_index_y,advected_index_z)]
          +a1*b1*c0*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x,advected_index_y,advected_index_z+1)]
          +a1*b0*c1*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x,advected_index_y+1,advected_index_z)]
          +a1*b0*c0*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x,advected_index_y+1,advected_index_z+1)]
          +a0*b1*c1*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x+1,advected_index_y,advected_index_z)]
          +a0*b1*c0*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x+1,advected_index_y,advected_index_z+1)]
          +a0*b0*c1*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x+1,advected_index_y+1,advected_index_z)]
          +a0*b0*c0*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x+1,advected_index_y+1,advected_index_z+1)];
}


//x成分のvelocityのmonotonic cubic interpolation
double monotonic_cubic_interpolation_x_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(double)(advected_index_z);
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-1){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-2;
        advected_z=(double)(advected_index_z);
    }
    //範囲外を参照しないようにする処置
    if(advected_index_x<=1||advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1
     ||advected_index_y<=1||advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-2
     ||advected_index_z<=1||advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-2){
         return linear_interpolation_x_3D(advected_x, advected_y, advected_z, all_grid);
    }
    //monotonic cubic interpolation の処理
    else{
        double dk_0;
        double dk_1;
        double delta_k;
        double interpolated_velocity_x0[4][4];
        double interpolated_velocity_x1[4];
        //z軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                dk_0  =  (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z-1)])/2.0;
                dk_1  =  (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+2)]
                         -all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)])/2.0;
                delta_k = all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                double c0 = advected_z-advected_index_z;
                //3つのスロープの符号が異なる場合
                if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                    interpolated_velocity_x0[i][j]=(1.0-c0)*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)]
                                                  +      c0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)];
    //                interpolated_velocity_x[i]=all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)];
                }

                else{
                    //論文の式は間違っている。正しくはこっち
                    interpolated_velocity_x0[i][j]=(dk_0+dk_1-2*delta_k)*(c0*c0*c0)+(3*delta_k-2*dk_0-dk_1)*(c0*c0)+dk_0*c0
                                              +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                }
            }
        }

        //y軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            dk_0  =  (interpolated_velocity_x0[i][2]-interpolated_velocity_x0[i][0])/2.0;
            dk_1  =  (interpolated_velocity_x0[i][3]-interpolated_velocity_x0[i][1])/2.0;
            delta_k = interpolated_velocity_x0[i][2]-interpolated_velocity_x0[i][1];
            double b0 = advected_y-advected_index_y;
            //3つのスロープの符号が異なる場合
            if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                interpolated_velocity_x1[i]=(1.0-b0)*interpolated_velocity_x0[i][1]
                                          +b0*interpolated_velocity_x0[i][2];
            }

            else{
                //論文の式は間違っている。正しくはこっち
                interpolated_velocity_x1[i]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                          +interpolated_velocity_x0[i][1];
            }

        }
        //x軸方向に関するinterpolation
        dk_0=(interpolated_velocity_x1[2]-interpolated_velocity_x1[0])/2.0;
        dk_1=(interpolated_velocity_x1[3]-interpolated_velocity_x1[1])/2.0;
        delta_k=(interpolated_velocity_x1[2]-interpolated_velocity_x1[1]);
        double a0=advected_x-advected_index_x;
        //3つのスロープの符号が異なる場合
        if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
            return (1.0-a0)*interpolated_velocity_x1[1]+a0*interpolated_velocity_x1[2];
        }
        else{
            //論文の式は間違っている。正しくはこっち
            return (dk_0+dk_1-2*delta_k)*(a0*a0*a0)+(3*delta_k-2*dk_0-dk_1)*(a0*a0)+dk_0*a0
                    +interpolated_velocity_x1[1];
        }
    }
}

//y成分のvelocityのmonotonic cubic interpolation
double monotonic_cubic_interpolation_y_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-1;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(double)(advected_index_z);
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-1){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-2;
        advected_z=(double)(advected_index_z);
    }
    //範囲外を参照しないようにする処置
    if(advected_index_x<=1||advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-2
     ||advected_index_y<=1||advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1
     ||advected_index_z<=1||advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-2){
         return linear_interpolation_y_3D(advected_x, advected_y, advected_z, all_grid);
    }
    //monotonic cubic interpolation の処理
    else{
        double dk_0;
        double dk_1;
        double delta_k;
        double interpolated_velocity_y0[4][4];
        double interpolated_velocity_y1[4];
        //z軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                dk_0  =  (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z-1)])/2.0;
                dk_1  =  (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+2)]
                         -all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)])/2.0;
                delta_k = all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                double c0 = advected_z-advected_index_z;
                //3つのスロープの符号が異なる場合
                if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                    interpolated_velocity_y0[i][j]=(1.0-c0)*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)]
                                                  +      c0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)];
    //                interpolated_velocity_x[i]=all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)];
                }

                else{
                    //論文の式は間違っている。正しくはこっち
                    interpolated_velocity_y0[i][j]=(dk_0+dk_1-2*delta_k)*(c0*c0*c0)+(3*delta_k-2*dk_0-dk_1)*(c0*c0)+dk_0*c0
                                              +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                }
            }
        }

        //y軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            dk_0  =  (interpolated_velocity_y0[i][2]-interpolated_velocity_y0[i][0])/2.0;
            dk_1  =  (interpolated_velocity_y0[i][3]-interpolated_velocity_y0[i][1])/2.0;
            delta_k = interpolated_velocity_y0[i][2]-interpolated_velocity_y0[i][1];
            double b0 = advected_y-advected_index_y;
            //3つのスロープの符号が異なる場合
            if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                interpolated_velocity_y1[i]=(1.0-b0)*interpolated_velocity_y0[i][1]
                                          +b0*interpolated_velocity_y0[i][2];
            }

            else{
                //論文の式は間違っている。正しくはこっち
                interpolated_velocity_y1[i]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                          +interpolated_velocity_y0[i][1];
            }

        }
        //x軸方向に関するinterpolation
        dk_0=(interpolated_velocity_y1[2]-interpolated_velocity_y1[0])/2.0;
        dk_1=(interpolated_velocity_y1[3]-interpolated_velocity_y1[1])/2.0;
        delta_k=(interpolated_velocity_y1[2]-interpolated_velocity_y1[1]);
        double a0=advected_x-advected_index_x;
        //3つのスロープの符号が異なる場合
        if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
            return (1.0-a0)*interpolated_velocity_y1[1]+a0*interpolated_velocity_y1[2];
        }
        else{
            //論文の式は間違っている。正しくはこっち
            return (dk_0+dk_1-2*delta_k)*(a0*a0*a0)+(3*delta_k-2*dk_0-dk_1)*(a0*a0)+dk_0*a0
                    +interpolated_velocity_y1[1];
        }
    }
}

//z成分のvelocityのmonotonic cubic interpolation
double monotonic_cubic_interpolation_z_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(double)(advected_index_z);
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-1;
        advected_z=(double)(advected_index_z);
    }
    //範囲外を参照しないようにする処置
    if(advected_index_x<=1||advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-2
     ||advected_index_y<=1||advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-2
     ||advected_index_z<=1||advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-1){
         return linear_interpolation_z_3D(advected_x, advected_y, advected_z, all_grid);
    }
    //monotonic cubic interpolation の処理
    else{
        double dk_0;
        double dk_1;
        double delta_k;
        double interpolated_velocity_z0[4][4];
        double interpolated_velocity_z1[4];
        //z軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                dk_0  =  (all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z-1)])/2.0;
                dk_1  =  (all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+2)]
                         -all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)])/2.0;
                delta_k = all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                double c0 = advected_z-advected_index_z;
                //3つのスロープの符号が異なる場合
                if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                    interpolated_velocity_z0[i][j]=(1.0-c0)*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)]
                                                  +      c0*all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)];
    //                interpolated_velocity_x[i]=all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)];
                }

                else{
                    //論文の式は間違っている。正しくはこっち
                    interpolated_velocity_z0[i][j]=(dk_0+dk_1-2*delta_k)*(c0*c0*c0)+(3*delta_k-2*dk_0-dk_1)*(c0*c0)+dk_0*c0
                                              +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                }
            }
        }

        //y軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            dk_0  =  (interpolated_velocity_z0[i][2]-interpolated_velocity_z0[i][0])/2.0;
            dk_1  =  (interpolated_velocity_z0[i][3]-interpolated_velocity_z0[i][1])/2.0;
            delta_k = interpolated_velocity_z0[i][2]-interpolated_velocity_z0[i][1];
            double b0 = advected_y-advected_index_y;
            //3つのスロープの符号が異なる場合
            if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                interpolated_velocity_z1[i]=(1.0-b0)*interpolated_velocity_z0[i][1]
                                          +b0*interpolated_velocity_z0[i][2];
            }

            else{
                //論文の式は間違っている。正しくはこっち
                interpolated_velocity_z1[i]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                          +interpolated_velocity_z0[i][1];
            }

        }
        //x軸方向に関するinterpolation
        dk_0=(interpolated_velocity_z1[2]-interpolated_velocity_z1[0])/2.0;
        dk_1=(interpolated_velocity_z1[3]-interpolated_velocity_z1[1])/2.0;
        delta_k=(interpolated_velocity_z1[2]-interpolated_velocity_z1[1]);
        double a0=advected_x-advected_index_x;
        //3つのスロープの符号が異なる場合
        if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
            return (1.0-a0)*interpolated_velocity_z1[1]+a0*interpolated_velocity_z1[2];
        }
        else{
            //論文の式は間違っている。正しくはこっち
            return (dk_0+dk_1-2*delta_k)*(a0*a0*a0)+(3*delta_k-2*dk_0-dk_1)*(a0*a0)+dk_0*a0
                    +interpolated_velocity_z1[1];
        }
    }
}


//advect項の計算
//速度場を時間 -dt だけバックトレースしてadvect項を計算する
void advect_fluid(Grid& all_grid){
    double velocity_after_advect[smoke_simulation::physical_const::kGrid_num_x+1][smoke_simulation::physical_const::kGrid_num_y+1][smoke_simulation::physical_const::kGrid_num_z+1][3];
    //velocityのx成分を計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x+1;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                double velocity_y ,velocity_z;
                //考えてるx面での速度のy成分の計算
                if(ix<=0){
                    velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)])/2.0;
                }
                else if(ix>=smoke_simulation::physical_const::kGrid_num_x){
                    velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy,iz)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy+1,iz)])/2.0;
                }
                else{
                    velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy,iz)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix-1,iy+1,iz)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)])/4.0;
                }
                //考えてるx面での速度z成分の計算
                if(ix<=0){
                    velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)])/2.0;
                }
                else if(ix>=smoke_simulation::physical_const::kGrid_num_x){
                    velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz+1)])/2.0;
                }
                else{
                    velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix-1,iy,iz+1)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)])/4.0;
                }
                //バックトレース先の位置
                double advected_x=(double)ix-((all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_y=(double)iy-((velocity_y*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_z=(double)iz-((velocity_z*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                //バックトレース先の速度を補間(linear interpolationの場合)
//                velocity_after_advect[ix][iy][iz][0]=linear_interpolation_x_3D(advected_x, advected_y, advected_z, all_grid);
                //バックトレース先の速度を補間(monotonic cubic interpolationの場合)
                velocity_after_advect[ix][iy][iz][0]=monotonic_cubic_interpolation_x_3D(advected_x, advected_y, advected_z, all_grid);
            }
        }
    }

    //velocityのy成分を計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y+1;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                double velocity_x, velocity_z;
                //考えてるy面での速度のx成分の計算
                if(iy<=0){
                    velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)])/2.0;
                }
                else if(iy>=smoke_simulation::physical_const::kGrid_num_y){
                    velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy-1,iz)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy-1,iz)])/2.0;
                }
                else{
                    velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy-1,iz)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy-1,iz)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)])/4.0;
                }
                //考えてるy面での速度のz成分の計算
                if(iy<=0){
                    velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)])/2.0;
                }
                else if(iy>=smoke_simulation::physical_const::kGrid_num_y){
                    velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz+1)])/2.0;
                }
                else{
                    velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy-1,iz+1)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)])/4.0;
                }
                //バックトレース先の位置
                double advected_x=(double)ix-((velocity_x*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_y=(double)iy-((all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_z=(double)iz-((velocity_z*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                //バックトレース先の速度を補間(linear interpolationの場合)
//                velocity_after_advect[ix][iy][iz][1]=linear_interpolation_y_3D(advected_x, advected_y, advected_z, all_grid);
                //バックトレース先の速度を補間(monotonic cubic interpolationの場合)
                velocity_after_advect[ix][iy][iz][1]=monotonic_cubic_interpolation_y_3D(advected_x, advected_y, advected_z, all_grid);
            }
        }
    }

    //velocityのz成分を計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z+1;iz++){
                double velocity_x, velocity_y;
                //考えてるz面での速度のx成分の計算
                if(iz<=0){
                    velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)])/2.0;
                }
                else if(iz>=smoke_simulation::physical_const::kGrid_num_z){
                    velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz-1)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz-1)])/2.0;
                }
                else{
                    velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz-1)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz-1)]
                               +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)])/4.0;
                }
                //考えてるz面での速度のy成分の計算
                if(iz<=0){
                    velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)])/2.0;
                }
                else if(iz>=smoke_simulation::physical_const::kGrid_num_z){
                    velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz-1)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz-1)])/2.0;
                }
                else{
                    velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz-1)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz-1)]
                               +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)])/4.0;
                }
                //バックトレース先の位置
                double advected_x=(double)ix-((velocity_x*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_y=(double)iy-((velocity_y*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_z=(double)iz-((all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                //バックトレース先の速度を補間(linear interpolationの場合)
//                velocity_after_advect[ix][iy][iz][2]=linear_interpolation_z_3D(advected_x, advected_y, advected_z, all_grid);
                //バックトレース先の速度を補間(monotonic cubic interpolationの場合)
                velocity_after_advect[ix][iy][iz][2]=monotonic_cubic_interpolation_z_3D(advected_x, advected_y, advected_z, all_grid);
            }
        }
    }

    //計算結果をコピー
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x+1;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]=velocity_after_advect[ix][iy][iz][0];
            }
        }
    }
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y+1;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]=velocity_after_advect[ix][iy][iz][1];
            }
        }
    }
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z+1;iz++){
                all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]=velocity_after_advect[ix][iy][iz][2];
            }
        }
    }
}

//Poisson eq. を解くことによる圧力の計算
void calc_pressure(Grid& all_grid){
    //グリッドの総数
    const int N=smoke_simulation::physical_const::kGrid_num_x
               *smoke_simulation::physical_const::kGrid_num_y
               *smoke_simulation::physical_const::kGrid_num_z;
    //係数行列の計算
    linear_algebra::sparse_matrix A(N,N);
//    linear_algebra::sparse_matrix_with_diagonal_element A(N,N);
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                if(ix-1>=0){
                    A.input_element(get_voxel_center_index_3D(ix, iy, iz),get_voxel_center_index_3D(ix-1, iy, iz),1);
                }
                if(iy-1>=0){
                    A.input_element(get_voxel_center_index_3D(ix, iy, iz),get_voxel_center_index_3D(ix, iy-1, iz),1);
                }
                if(iz-1>=0){
                    A.input_element(get_voxel_center_index_3D(ix, iy, iz),get_voxel_center_index_3D(ix, iy, iz-1),1);
                }
                A.input_element(get_voxel_center_index_3D(ix, iy, iz),get_voxel_center_index_3D(ix, iy, iz),-6);
                if(iz+1<=smoke_simulation::physical_const::kGrid_num_z-1){
                    A.input_element(get_voxel_center_index_3D(ix, iy, iz),get_voxel_center_index_3D(ix, iy, iz+1),1);
                }
                if(iy+1<=smoke_simulation::physical_const::kGrid_num_y-1){
                    A.input_element(get_voxel_center_index_3D(ix, iy, iz),get_voxel_center_index_3D(ix, iy+1, iz),1);
                }
                if(ix+1<=smoke_simulation::physical_const::kGrid_num_x-1){
                    A.input_element(get_voxel_center_index_3D(ix, iy, iz),get_voxel_center_index_3D(ix+1, iy, iz),1);
                }
            }
        }
    }
    //連立方程式の右辺のベクトルの計算
    std::vector<double> b(N);
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                 b[get_voxel_center_index_3D(ix, iy, iz)]
                    =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)]-all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]
                     +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)]-all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]
                     +all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)]-all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)])
                    *smoke_simulation::physical_const::kCell_length
                    /(smoke_simulation::physical_const::kDt);
            }
        }
    }
/*
    //時間計測用
    std::chrono::system_clock::time_point  start, end;
    start = std::chrono::system_clock::now(); // 時間計測開始
*/
    //CG法により圧力場を得る
//    linear_algebra::conjugate_gradient(A, b, all_grid.pressure, N, 10000, 0.0001);
    linear_algebra::incomplete_cholesky_conjugate_gradient(A, b, all_grid.pressure, N, 10000, 0.0001);
/*
    end = std::chrono::system_clock::now();  // 時間計測終了
    std::ofstream writing_file;
    writing_file.open("length_of_time_ICCG_64_3d.dat", std::ios::app);
    writing_file << (std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()) << std::endl;
    writing_file.close();
*/
    //gauss seidel法を使う場合(Aはsparse_matrix_with_diagonal_elementにする)
//    gauss_seidel(A,b,all_grid.pressure,N,200);
}


//pressure gradient termの計算
void calc_pressure_gradient_term(Grid& all_grid){
    //圧力の計算
    calc_pressure(all_grid);
    //圧力場に境界条件をセット
    set_boundary_pressure(all_grid);
    //pressure gradeint term によって速度場を更新
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
            for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
                all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]-=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz)]-all_grid.pressure[get_voxel_center_index_3D(ix-1, iy, iz)])
                *smoke_simulation::physical_const::kDt
                /(smoke_simulation::physical_const::kCell_length);
            }
        }
    }
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z-1;iz++){
                all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]-=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz)]-all_grid.pressure[get_voxel_center_index_3D(ix, iy-1, iz)])
                *smoke_simulation::physical_const::kDt
                /(smoke_simulation::physical_const::kCell_length);
            }
        }
    }
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x-1;ix++){
        for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y-1;iy++){
            for(int iz=1;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]-=(all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz)]-all_grid.pressure[get_voxel_center_index_3D(ix, iy, iz-1)])
                *smoke_simulation::physical_const::kDt
                /(smoke_simulation::physical_const::kCell_length);
            }
        }
    }
}

//流体の 1 time step
void update_fluid_velocity(Grid& all_grid){
//固定境界条件の場合
//    add_source(all_grid);
    add_force_fluid(all_grid);
    set_boundary_velocity(all_grid);
    advect_fluid(all_grid);
    set_boundary_velocity(all_grid);
    calc_pressure_gradient_term(all_grid);
    set_boundary_velocity(all_grid);
}
}
