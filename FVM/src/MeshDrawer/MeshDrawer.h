#pragma once

#include "Mesh/MeshBase.h"


class MeshDrawer
{
public:

    MeshDrawer(MeshBase const& mesh);

    void show();

private:

    MeshBase const& m_mesh;
};
