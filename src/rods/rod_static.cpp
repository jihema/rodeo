/*
 * rod_static.cpp
 */

#include "rod_static.h"

namespace rodeo
{

void StaticRodState<0>::compute_material_frames()
{
    for (size_t i = 0; i < m_material_frames.size(); ++i)
    {
        m_material_frames[i].col(0) = m_reference_frames[i].col(0); // Tangents are equal.
        m_material_frames[i].rightCols<2>() =
                m_reference_frames[i].rightCols<2>()
                        * Eigen::Rotation2Dd(get_torsion(i)).toRotationMatrix();
    }
}

}

