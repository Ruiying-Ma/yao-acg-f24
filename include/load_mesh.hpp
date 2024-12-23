#ifndef LOAD_MESH_HPP_
#define LOAD_MESH_HPP_


#include <vector>
#include <memory>
#include <sys/stat.h>
#include <fstream>
#include <cstdint>
#include <unistd.h> // For getcwd
#include <limits.h> // For PATH_MAX
#include "../external/Vertex.hpp"
#include "../external/Triangle.hpp"

typedef struct stat64 NativeStatStruct;
typedef std::shared_ptr<std::istream> InputStreamHandle;

template<typename T>
inline void streamRead(std::istream &in, T &dst)
{
    in.read(reinterpret_cast<char *>(&dst), sizeof(T));
}

template<typename T>
inline void streamRead(std::istream &in, std::vector<T> &dst)
{
    in.read(reinterpret_cast<char *>(&dst[0]), dst.size()*sizeof(T));
}

bool TungstenloadWo3(const std::string &relative_path, std::vector<Tungsten::Vertex> &verts, std::vector<Tungsten::TriangleI> &tris) {
    
    char buffer[PATH_MAX];

    // Get the current working directory
    if (getcwd(buffer, sizeof(buffer)) == nullptr) {
        return false;
    }
    
    NativeStatStruct info;
    std::string absolute_dir(buffer);
    std::string absolute_path = "";

    if (!absolute_dir.empty() && absolute_dir[absolute_dir.length() - 1] == '/') {
        absolute_path = absolute_dir + relative_path;
    } else {
        absolute_path = absolute_dir + "/" + relative_path;
    }

    if (stat64(absolute_path.c_str(), &info) == 0) {
        std::shared_ptr<std::istream> in(new std::ifstream(absolute_path, std::ios_base::in | std::ios_base::binary));
        if (!in->good()) {return false;}
        std::uint64_t numVerts, numTris;
        streamRead(*in, numVerts);
        verts.resize(size_t(numVerts));
        streamRead(*in, verts);
        streamRead(*in, numTris);
        tris.resize(size_t(numTris));
        streamRead(*in, tris);
        return true;
    }
    return false;
}

#endif
