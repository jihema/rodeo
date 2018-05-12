/*
 * rod_state.h
 */

#pragma once

#include "../rodmath.h"

#include <vector>

namespace rodeo
{

class ExplicitEulerSolver;

namespace force
{
class ForceBase;
}

class DirtyAccess: public std::runtime_error
{
public:
    DirtyAccess(std::string const& what) :
            std::runtime_error(what)
    {
    }
};

class RodBase
{
public:

    RodBase(VecXd const& rest_dofs) :
            m_rest_dofs(&rest_dofs)
    {
    }

    inline size_t num_dofs() const
    {
        return m_rest_dofs->size();
    }

    inline size_t num_vertices() const
    {
        return (m_rest_dofs->size() + 1) / 4;
    }

    inline DiagonalMatXd const& get_inverse_mass() const
    {
        return m_inverse_mass;
    }

    inline DiagonalMatXd const& get_mass() const
    {
        return m_mass;
    }

    inline double get_vertex_mass(size_t i) const
    {
        assert(i < num_vertices());

        return 1.; // FIXME
    }

    const VecXd* get_rest_dofs() const
    {
        return m_rest_dofs;
    }

    double get_rest_length(size_t i) const
    {
        assert(i < num_vertices() - 1);

        return (m_rest_dofs->segment<3>(4 * i + 4)
                - m_rest_dofs->segment<3>(4 * i)).norm();
    }

protected:

    // Fixed quantities.
    VecXd const* m_rest_dofs = nullptr;
    DiagonalMatXd m_mass = DiagonalMatXd(VecXd::Ones(m_rest_dofs->size())); // FIXME. In fact those should be stored with the original rest_dofs.
    DiagonalMatXd m_inverse_mass = DiagonalMatXd(
            VecXd::Ones(m_rest_dofs->size())); // FIXME
};

class RodStateBase: public RodBase
{
public:
    RodStateBase(VecXd const& rest_dofs) :
            RodBase(rest_dofs), //
            m_dofs(rest_dofs.size()), //
            m_reference_frames((rest_dofs.size() - 3) / 4)
    {
    }

    void set_dofs_to_rest()
    {
        set_dofs(nullptr, *m_rest_dofs);
    }

    /**
     * \brief Sets dofs to a new value, updates reference and material frames.
     *
     * The torsion dofs (indices 4 * i + 3) are always defined relative to some reference frames.
     *
     * If reference_frames is null, that means that the torsion dofs are defined with respect to
     * reference frames that are parallel-transported from an initial, arbitrary reference frame.
     * That one is determined from the initial tangent by the build_orthonormal_basis() algorithm.
     *
     * If reference_frames is not null, that means that the torsion dofs are defined with respect to
     * time-parallel transported reference frames (from their previous position,
     * following the move of the tangent vector).
     */
    void set_dofs(std::vector<Mat3d> const* const reference_frames,
            VecXd const& dofs);

    inline VecXd const& get_dofs() const
    {
        return m_dofs;
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

    const std::vector<Mat3d>& get_reference_frames() const
    {
        return m_reference_frames;
    }

protected:

    void init_reference_frames();
    void transport_reference_frames(std::vector<Mat3d> const& reference_frames,
            VecXd const& future_dofs);

    // State quantities.
    VecXd m_dofs;
    std::vector<Mat3d> m_reference_frames;

    std::vector<force::ForceBase const*> m_forces;

    bool m_dirty = true;
};

template<unsigned int n> class StaticRodState;

template<> class StaticRodState<0> : public RodStateBase
{
public:

    StaticRodState(VecXd const& rest_dofs) :
            RodStateBase(rest_dofs), //
            m_material_frames((rest_dofs.size() - 3) / 4)
    {
    }

    virtual ~StaticRodState()
    {
    }

    virtual void compute();

    inline double get_energy() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_energy;
    }

    void accumulate_energy(force::ForceBase & force);

protected:

    void compute_material_frames();

    // Auxiliary quantities, can be purged if we're running short on memory.
    std::vector<Mat3d> m_material_frames;

    // Computed quantities.
    double m_energy = 0.;
};

template<> class StaticRodState<1> : public StaticRodState<0>
{
public:

    StaticRodState(VecXd const& rest_dofs) :
            StaticRodState<0>(rest_dofs), //
            m_force(rest_dofs.size())
    {
    }

    void compute() override;

    inline VecXd const& get_force() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_force;
    }

    void accumulate_force(force::ForceBase & force);

protected:

    VecXd m_force;
};

template<> class StaticRodState<2> : public StaticRodState<1>
{
public:

    StaticRodState(VecXd const& rest_dofs) :
            StaticRodState<1>(rest_dofs), m_jacobian(rest_dofs.size(),
                    rest_dofs.size())
    {
    }

    void compute() override;

    inline BandLimitedMatXd const& get_jacobian() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_jacobian;
    }

    void accumulate_jacobian(force::ForceBase & force);

protected:
    BandLimitedMatXd m_jacobian;

    friend class ForceBase;
};

}

