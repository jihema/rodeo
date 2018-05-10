/*
 * rod_state.h
 *
 *  Created on: 9 May 2018
 *      Author: jau
 */

#pragma once

#include "rodmath.h"

#include <stddef.h>
#include <assert.h>
#include <stdexcept>
#include <vector>

namespace rodeo
{

class DirtyAccess: public std::runtime_error
{
public:
    DirtyAccess(std::string const& what) :
            std::runtime_error(what)
    {
    }
};

class RodState
{
public:

    RodState(VecXd const& rest_dofs) :
            m_rest_dofs(&rest_dofs), //
            m_dirty(true), //
            m_dofs(rest_dofs.size()), //
            m_reference_frames((rest_dofs.size() - 3) / 4), //
            m_material_frames((rest_dofs.size() - 3) / 4), //
            m_force(rest_dofs.size()), //
            m_jacobian(rest_dofs.size(), rest_dofs.size())
    {
    }

    void set_dofs_to_rest()
    {
        set_dofs(*m_rest_dofs, true);
    }

    /**
     * \brief Sets dofs to a new value, updates reference and material frames.
     *
     * The torsion dofs (indices 4 * i + 3) are always defined relative to some reference frames.
     *
     * If init_frames is true, that means that the torsion dofs are defined with respect to
     * reference frames that are parallel-transported from an initial, arbitrary reference frame.
     * That one is determined from the initial tangent by the build_orthonormal_basis() algorithm.
     *
     * If init_frames is false, that means that the torsion dofs are defined with respect to
     * time-parallel transported reference frames (from their previous position,
     * following the move of the tangent vector).
     *
     * In both cases, reference and material frames are updated when setting the dofs.
     */
    void set_dofs(VecXd const& dofs, bool const init_frames);

    inline VecXd const& get_dofs() const
    {
        return m_dofs;
    }

    inline size_t num_vertices() const
    {
        return (m_dofs.size() + 1) / 4;
    }

    inline Vec3d get_vertex(size_t i) const
    {
        assert(i < num_vertices());

        return m_dofs.segment<3>(4 * i);
    }

    inline double get_torsion(size_t i)
    {
        return m_dofs(4 * i + 3);
    }

    inline bool is_dirty() const
    {
        return m_dirty;
    }

    void compute();

    inline double get_energy() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_energy;
    }

    inline VecXd const& get_force() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_force;
    }

    inline BandLimitedMatXd const& get_jacobian() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_jacobian;
    }

private:

    void init_reference_frames();
    void transport_reference_frames(VecXd const& future_dofs);
    void compute_material_frames();

    VecXd const* m_rest_dofs = nullptr;
    bool m_dirty = true;
    VecXd m_dofs;

    // Auxiliary data, can be purged if we're running short on memory.
    std::vector<Mat3d> m_reference_frames;
    std::vector<Mat3d> m_material_frames;

    double m_energy = 0.;
    VecXd m_force;
    BandLimitedMatXd m_jacobian;
};

}

