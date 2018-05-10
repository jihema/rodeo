/*
 * rod_state.cpp
 *
 *  Created on: 9 May 2018
 *      Author: jau
 */

#include "rod_state.h"
#include <assert.h>

namespace rodeo
{

void StaticRodState<0>::set_dofs(VecXd const& dofs, bool const init_frames)
{
    m_dirty = true;

    if (init_frames) // Initialising.
    {
        m_dofs = dofs;
        init_reference_frames();
    }
    else // Evolving.
    {
        transport_reference_frames(dofs);
        m_dofs = dofs;
    }

    compute_material_frames();
}

void StaticRodState<0>::init_reference_frames()
{
    Vec3d tangent_0 = (get_vertex(1) - get_vertex(0)).normalized();
    m_reference_frames[0] = build_orthonormal_basis(tangent_0);

    for (size_t i = 1; i < m_reference_frames.size(); ++i)
    {
        Vec3d tangent_1 = (get_vertex(i + 1) - get_vertex(i)).normalized();
        m_reference_frames[i] = parallel_transport_matrix(tangent_0, tangent_1) * m_reference_frames[i - 1]; // NB col(0) should coincide with tangent_1.
        tangent_0 = tangent_1;
    }
}

void StaticRodState<0>::transport_reference_frames(VecXd const& future_dofs)
{
    for (size_t i = 0; i < m_reference_frames.size(); ++i)
    {
        auto const& tangent_0 = m_reference_frames[i].col(0);

        Vec3d const tangent_1 = (future_dofs.segment<3>(4 * i + 4) - future_dofs.segment<3>(4 * i)).normalized();

        m_reference_frames[i] = parallel_transport_matrix(tangent_0, tangent_1) * m_reference_frames[i]; // NB col(0) should coincide with tangent_1.
    }
}

void StaticRodState<0>::compute_material_frames()
{
    for (size_t i = 0; i < m_material_frames.size(); ++i)
    {
        m_material_frames[i].col(0) = m_reference_frames[i].col(0); // Tangents are equal.
        m_material_frames[i].rightCols<2>() = m_reference_frames[i].rightCols<2>()
                * Eigen::Rotation2Dd(get_torsion(i)).toRotationMatrix();
    }
}

void StaticRodState<0>::compute()
{
    assert(m_rest_dofs != nullptr);
    assert(m_dofs.size() == m_rest_dofs->size());

    if (!m_dirty)
    {
        return;
    }

    // Compute energy.

    m_dirty = false;
}

void StaticRodState<1>::compute()
{
    assert(m_rest_dofs != nullptr);
    assert(m_dofs.size() == m_rest_dofs->size());

    if (!m_dirty)
    {
        return;
    }

    // Compute energy and force.

    m_dirty = false;
}

void StaticRodState<2>::compute()
{
    assert(m_rest_dofs != nullptr);
    assert(m_dofs.size() == m_rest_dofs->size());

    if (!m_dirty)
    {
        return;
    }

    // Compute energy, force and Jacobian.

    m_dirty = false;
}

}
