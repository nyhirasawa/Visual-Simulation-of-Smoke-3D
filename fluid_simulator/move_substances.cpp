#include "move_substances.h"

#include <iostream>

#include "linear_solver.h"
#include "physical_const.h"
#include "utils.h"

namespace smoke_simulation{
//外力項の計算
void add_force_substances(Grid& all_grid){

}
/*
//substance densityのlinear interpolation
double linear_interpolation_substances(double advected_x, double advected_y, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(advected_index_x);
//        advected_x=(advected_index_x+0.5);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(advected_index_x);
//        advected_x=(advected_index_x+0.5);
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(advected_index_y);
//        advected_y=(advected_index_y+0.5);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(advected_index_y);
//        advected_y=(advected_index_y+0.5);
    }
    double a0, a1, b0, b1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
    return a1*(b1*all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y)]
          +b0*all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y+1)])
          +a0*(b1*all_grid.substance_density[get_voxel_center_index(advected_index_x+1,advected_index_y)]
          +b0*all_grid.substance_density[get_voxel_center_index(advected_index_x+1,advected_index_y+1)]);
}
*/
//substance densityのlinear interpolation
double linear_interpolation_substances_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(advected_index_x);
//        advected_x=(advected_index_x+0.5);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(advected_index_x);
//        advected_x=(advected_index_x+0.5);
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(advected_index_y);
//        advected_y=(advected_index_y+0.5);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(advected_index_y);
//        advected_y=(advected_index_y+0.5);
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(advected_index_z);
//        advected_z=(advected_index_z+0.5);
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-1){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-2;
        advected_z=(advected_index_z);
//        advected_z=(advected_index_z+0.5);
    }
    double a0, a1, b0, b1, c0, c1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
    c0=advected_z-advected_index_z;
    c1=1.0-c0;
    return a1*b1*c1*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x,advected_index_y,advected_index_z)]
          +a1*b1*c0*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x,advected_index_y,advected_index_z+1)]
          +a1*b0*c1*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x,advected_index_y+1,advected_index_z)]
          +a1*b0*c0*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x,advected_index_y+1,advected_index_z+1)]
          +a0*b1*c1*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x+1,advected_index_y,advected_index_z)]
          +a0*b1*c0*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x+1,advected_index_y,advected_index_z+1)]
          +a0*b0*c1*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x+1,advected_index_y+1,advected_index_z)]
          +a0*b0*c0*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x+1,advected_index_y+1,advected_index_z+1)];
}
//substance densityのmonotonic cubic interpolation
double monotonic_cubic_substances_3D(double advected_x, double advected_y, double advected_z, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    int advected_index_z=(int)(advected_z);
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)(advected_index_x);
//        advected_x=(double)(advected_index_x+0.5);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(double)(advected_index_x);
//        advected_x=(double)(advected_index_x+0.5);
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)(advected_index_y);
//        advected_y=(double)(advected_index_y+0.5);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(double)(advected_index_y);
//        advected_y=(double)(advected_index_y+0.5);
    }
    if(advected_index_z<1){
        advected_index_z=1;
        advected_z=(double)(advected_index_z);
//        advected_z=(double)(advected_index_z+0.5);
    }
    if(advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-1){
        advected_index_z=smoke_simulation::physical_const::kGrid_num_z-2;
        advected_z=(double)(advected_index_z);
//        advected_z=(double)(advected_index_z+0.5);
    }
    //範囲外を参照しないようにする処置
    if(advected_index_x<=1||advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-2
     ||advected_index_y<=1||advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-2
     ||advected_index_z<=1||advected_index_z>=smoke_simulation::physical_const::kGrid_num_z-2){
         return linear_interpolation_substances_3D(advected_x, advected_y, advected_z, all_grid);
    }
    //monotonic cubic interpolation の処理
    else{
        double dk_0;
        double dk_1;
        double delta_k;
        double interpolated_substance_density0[4][4];
        //z軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++){
                dk_0  =  (all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z-1)])/2.0;
                dk_1  =  (all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+2)]
                         -all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)])/2.0;
                delta_k = all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)]
                         -all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                double b0 = advected_y-advected_index_y;
                if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                    interpolated_substance_density0[i][j]=(1.0-b0)*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)]
                                                         +b0      *all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z+1)];
                }
                else{
                    //論文の式は間違っている。正しくはこっち
                    interpolated_substance_density0[i][j]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                                         +all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1+j, advected_index_z)];
                }
            }
        }
        double interpolated_substance_density1[4];
        //y軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            dk_0  =  (interpolated_substance_density0[i][2]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y+1)]*/
                     -interpolated_substance_density0[i][0]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y-1)]*/)/2.0;
            dk_1  =  (interpolated_substance_density0[i][3]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y+2)]*/
                     -interpolated_substance_density0[i][1]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y)]*/)/2.0;
            delta_k = interpolated_substance_density0[i][2]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y+1)]*/
                     -interpolated_substance_density0[i][1]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y)]*/;
            double b0 = advected_y-advected_index_y;
            if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                interpolated_substance_density1[i]=(1.0-b0)*interpolated_substance_density0[i][1]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i,advected_index_y)]*/
                                                  +b0      *interpolated_substance_density0[i][2]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i,advected_index_y+1)]*/;
            }
            else{
                //論文の式は間違っている。正しくはこっち
                interpolated_substance_density1[i]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                                 +interpolated_substance_density0[i][1]/*all_grid.substance_density[get_voxel_center_index_3D(advected_index_x-1+i, advected_index_y)]*/;
            }
        }
        //x軸方向に関するinterpolation
        dk_0=(interpolated_substance_density1[2]-interpolated_substance_density1[0])/2.0;
        dk_1=(interpolated_substance_density1[3]-interpolated_substance_density1[1])/2.0;
        delta_k=(interpolated_substance_density1[2]-interpolated_substance_density1[1]);
        double a0=advected_x-advected_index_x;
        //3つのスロープの符号が異なる場合
        if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
            return (1.0-a0)*interpolated_substance_density1[1]
                   +a0     *interpolated_substance_density1[2];
        }
        else{
            //論文の式は間違っている。正しくはこっち
            return (dk_0+dk_1-2*delta_k)*(a0*a0*a0)+(3*delta_k-2*dk_0-dk_1)*(a0*a0)+dk_0*a0
                    +interpolated_substance_density1[1];
        }
    }
}



//流体の速度場によってsubstanceが運ばれる項(流体のadvect項に相当)
//速度場を時間 -dt だけバックトレースして計算する
void transport_substances(Grid& all_grid){
    double substance_density_after_advect[smoke_simulation::physical_const::kGrid_num_x][smoke_simulation::physical_const::kGrid_num_y][smoke_simulation::physical_const::kGrid_num_z];
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                double velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix,iy,iz)]+all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x_3D(ix+1,iy,iz)])/2.0;
                double velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy,iz)]+all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y_3D(ix,iy+1,iz)])/2.0;
                double velocity_z=(all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz)]+all_grid.velocity_in_voxel_face_z[get_voxel_face_index_z_3D(ix,iy,iz+1)])/2.0;
                //バックトレース先の座標
                double advected_x=ix-((velocity_x*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_y=iy-((velocity_y*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                double advected_z=iz-((velocity_z*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
                //バックトレース先の座標のindex
                //バックトレース先の速度を補間する
//                substance_density_after_advect[ix][iy][iz]=linear_interpolation_substances_3D(advected_x, advected_y, advected_z, all_grid);
                substance_density_after_advect[ix][iy][iz]=monotonic_cubic_substances_3D(advected_x, advected_y, advected_z, all_grid);
            }
        }
    }
    //計算結果をコピー
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int iz=0;iz<smoke_simulation::physical_const::kGrid_num_z;iz++){
                all_grid.substance_density[get_voxel_center_index_3D(ix,iy,iz)]=substance_density_after_advect[ix][iy][iz];
            }
        }
    }
}

//上の4ステップをまとめただけの関数(substance densityの1時間ステップ分の更新に相当)
void move_substances(Grid& all_grid){
//    add_force_substances(all_grid);
//    diffuse_substances(all_grid);
    transport_substances(all_grid);
//    dissipate_substances(all_grid);
}
}
