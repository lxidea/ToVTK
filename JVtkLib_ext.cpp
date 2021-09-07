#include "JVtkLib.h"
#include "Functions.h"
#include <fstream>

template<typename T>
void SwapEnd(T& var)
{
  char* varArray = reinterpret_cast<char*>(&var);
  for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
    std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

void JVtkLib::SaveVtkGridData(std::string fname,
    const JDataArrays &arrays,std::string posfield,
    const JDataArrays &arrays2,bool createpath){
    if(createpath) {
        // create path;
    }

    const unsigned npok = arrays.GetDataCount();
    const unsigned clok = arrays2.GetDataCount();
    //printf("DataCount: %u, CellCount: %u", npok,clok);

    std::ofstream file;
    file.open(fname, std::ios::out|std::ios::binary);

    file << "# vtk DataFile Version 2.0" << std::endl
        << "Comment if needed" << std::endl;
    
    file << "BINARY"<< std::endl << std::endl;

    bool drop=false;
    if (!arrays.ExistsName(posfield))
        Run_ExceptioonSta(fun::PrintStr("Count not find data named :%s",posfield));
    const tfloat3* pos = arrays.GetArrayFloat3(posfield);
    const uint* idp = arrays.GetArrayUint("Idp");
    const tfloat3* vel = arrays.GetArrayFloat3("Vel");
    const float* mass = arrays.GetArrayFloat("Mass");
    const float* rhop = arrays.GetArrayFloat("Rhop");
    const byte * type = arrays.GetArrayByte("Type");
    const float* pres = arrays.GetArrayFloat("Press");
    const byte*    mk = arrays.GetArrayByte("Mk");
    const float * dropr = NULL;
    const float * grav = NULL;
    if (arrays.ExistsName("DropletR")) {
        drop = true;
        dropr = arrays.GetArrayFloat("DropletR");
        grav = arrays.GetArrayFloat("GravityCoeff");
    }
    const tuint8* conn = arrays2.GetArrayUint8("Connectivity");

    // start writing
    file << "DATASET UNSTRUCTURED_GRID" << std::endl
    << "POINTS " << npok << " float" << std::endl;

    for(unsigned i=0;i<npok;i++) {
        float rx = pos[i].x;
        float ry = pos[i].y;
        float rz = pos[i].z;
        SwapEnd(rx);
        SwapEnd(ry);
        SwapEnd(rz);
        file.write(reinterpret_cast<char*>(&rx), sizeof(float));
        file.write(reinterpret_cast<char*>(&ry), sizeof(float));
        file.write(reinterpret_cast<char*>(&rz), sizeof(float));
    }
    file << std::endl;// end position writing

    file << "CELLS " << clok << " " << clok*9 << std::endl;
    for(unsigned i=0;i<clok;i++){
        const ttuint8 rec = ReorderCellTUint8(conn[i]);
        int num = 8;
        SwapEnd(num);
        file.write(reinterpret_cast<char*>(&num), sizeof(int));
        for(unsigned j=0;j<8;j++) {
            int id = rec.id[j];
            SwapEnd(id);
            file.write(reinterpret_cast<char*>(&id), sizeof(int));
        }
    }
    file << std::endl << std::endl;
    file << "CELL_TYPES " << clok << std::endl;
    int t = 11;
    SwapEnd(t);
    for(unsigned i=0;i<clok;i++)
        file.write(reinterpret_cast<char*>(&t), sizeof(int));
    
    file << std::endl;

    // writing Idp
    file << "POINT_DATA " << npok << std::endl
    << "SCALARS Idp unsigned_int" << std::endl
    << "LOOKUP_TABLE default" << std::endl;

    for(uint i=0;i<npok;i++) {
        unsigned ui = idp[i];
        SwapEnd(ui);
        file.write(reinterpret_cast<char*>(&ui), sizeof(unsigned));
    }
    file << std::endl;

    // writing vel
    file << "VECTORS Vel float" << std::endl;

    for(uint i=0;i<npok;i++) {
        float vx = vel[i].x;
        float vy = vel[i].y;
        float vz = vel[i].z;

        SwapEnd(vx);
        SwapEnd(vy);
        SwapEnd(vz);

        file.write(reinterpret_cast<char*>(&vx), sizeof(float));
        file.write(reinterpret_cast<char*>(&vy), sizeof(float));
        file.write(reinterpret_cast<char*>(&vz), sizeof(float));
    }
    file << std::endl;

    // writing mass
    file << "FIELD FieldData " << (drop?6:4) << std::endl
    << "Mass 1 " << npok << " float" << std::endl;

    for(uint i=0;i<npok;i++) {
        float m = mass[i];
        SwapEnd(m);
        file.write(reinterpret_cast<char*>(&m), sizeof(float));
    }
    file << std::endl << std::endl;

    // writing rhop
    file << "Rhop 1 " << npok << " float" << std::endl;

    for(uint i=0;i<npok;i++) {
        float r = rhop[i];
        SwapEnd(r);
        file.write(reinterpret_cast<char*>(&r), sizeof(float));
    }
    file << std::endl;

    //unsigned_char
    // writing type
    // file << "POINT_DATA " << npok << std::endl
    // << "SCALARS Type unsigned_char " << npok << std::endl
    // << "LOOKUP_TABLE default" << std::endl;

    // for(uint i=0;i<npok;i++) {
    //     byte t = type[i];
    //     SwapEnd(t);
    //     file.write(reinterpret_cast<char*>(&t), sizeof(byte));
    // }
    // file << std::endl;

    // writing press
    file << "Press 1 " << npok << " float" << std::endl;

    for(uint i=0;i<npok;i++) {
        float p = pres[i];
        SwapEnd(p);
        file.write(reinterpret_cast<char*>(&p), sizeof(float));
    }
    file << std::endl;

    // writing mk
    file << "Mk 1 " << npok << " unsigned_char" << std::endl;

    for(uint i=0;i<npok;i++) {
        byte m = mk[i];
        SwapEnd(m);
        file.write(reinterpret_cast<char*>(&m), sizeof(byte));
    }
    file << std::endl;

    if(drop) {
        // writing dropletr
        file << "DropletR 1 " << npok << " float" << std::endl;

        for(uint i=0;i<npok;i++) {
            float d = dropr[i];
            SwapEnd(d);
            file.write(reinterpret_cast<char*>(&d), sizeof(float));
        }
        file << std::endl;

        // writing gravitycoeff
        file << "GravityCoeff 1 " << npok << " float" << std::endl;

        for(uint i=0;i<npok;i++) {
            float g = grav[i];
            SwapEnd(g);
            file.write(reinterpret_cast<char*>(&g), sizeof(float));
        }
        file << std::endl;
    }

    file.close();
}