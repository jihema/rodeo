/*
 * rod_static.h
 */

#pragma once

#include "rod_state.h"

namespace rodeo
{

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

    /**
     * Resets computed quantities (energy,...) to zero and computes auxiliaries (material frames...).
     */
    virtual void pre_compute()
    {
        compute_material_frames();
        m_potential_energy = 0;
    }

    inline double get_potential_energy() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_potential_energy;
    }

protected:

    void print(std::ostream& os) const
     {
         os << "potential energy = " << get_potential_energy();
     }

    void compute_material_frames();

    // Auxiliary quantities, can be purged if we're running short on memory.
    std::vector<Mat3d> m_material_frames;

    // Computed quantities.
    double m_potential_energy = 0.;

    friend class force::ForceBase;
};

template<> class StaticRodState<1> : public StaticRodState<0>
{
public:

    StaticRodState(VecXd const& rest_dofs) :
            StaticRodState<0>(rest_dofs), //
            m_force_vector(rest_dofs.size())
    {
    }

    virtual void pre_compute()
    {
        StaticRodState<0>::pre_compute();
        m_force_vector.setZero();
    }

    inline VecXd const& get_force() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_force_vector;
    }

protected:

    // Computed quantities.
    VecXd m_force_vector;

    friend class force::ForceBase;
};

template<> class StaticRodState<2> : public StaticRodState<1>
{
public:

    StaticRodState(VecXd const& rest_dofs) :
            StaticRodState<1>(rest_dofs), m_jacobian(rest_dofs.size(),
                    rest_dofs.size())
    {
    }

    virtual void pre_compute()
    {
        StaticRodState<1>::pre_compute();
        m_jacobian.setZero();
    }

    inline BandLimitedMatXd const& get_jacobian() const
    {
        if (m_dirty)
        {
            throw DirtyAccess(__PRETTY_FUNCTION__);
        }

        return m_jacobian;
    }

protected:

    // Computed quantities.
    BandLimitedMatXd m_jacobian;

    friend class force::ForceBase;
};

}
