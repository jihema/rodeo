/*
 * rod_state.cpp
 */

#include "rod_state.h"

#include "../forces/force_base.h"

namespace rodeo
{

void RodStateBase::set_dofs(std::vector<Mat3d> const* const reference_frames,
        VecXd const& dofs)
{
    m_dirty = true;

    if (!reference_frames) // Initialising.
    {
        m_dofs = dofs;
        init_reference_frames();
    }
    else // Evolving.
    {
        transport_reference_frames(*reference_frames, dofs);
        m_dofs = dofs;
    }
}

void RodStateBase::init_reference_frames()
{
    Vec3d tangent_0 = (get_vertex(1) - get_vertex(0)).normalized();
    m_reference_frames[0] = build_orthonormal_basis(tangent_0);

    for (size_t i = 1; i < m_reference_frames.size(); ++i)
    {
        Vec3d tangent_1 = (get_vertex(i + 1) - get_vertex(i)).normalized();
        m_reference_frames[i] = parallel_transport_matrix(tangent_0, tangent_1)
                * m_reference_frames[i - 1]; // NB col(0) should coincide with tangent_1.
        tangent_0 = tangent_1;
    }
}

void RodStateBase::transport_reference_frames(
        std::vector<Mat3d> const& reference_frames, VecXd const& future_dofs)
{
    assert(reference_frames.size() * 4 + 3 == (size_t ) future_dofs.size());

    for (size_t i = 0; i < reference_frames.size(); ++i)
    {
        auto const& tangent_0 = reference_frames[i].col(0);

        Vec3d const tangent_1 = (future_dofs.segment<3>(4 * i + 4)
                - future_dofs.segment<3>(4 * i)).normalized();

        m_reference_frames[i] = parallel_transport_matrix(tangent_0, tangent_1)
                * reference_frames[i]; // NB col(0) should coincide with tangent_1.
    }
}



}
