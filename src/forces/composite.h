/*
 * composite.h
 */

#pragma once

#include "force_base.h"
#include <vector>

namespace rodeo
{

namespace force
{

class Composite: public ForceBase
{
public:

    Composite(std::vector<ForceBase*> const& components =
            std::vector<ForceBase*>()) :
            m_components(components)
    {
        if (!m_components.empty())
        {
            m_rod = m_components.front()->m_rod;
        }
        for (ForceBase const* component : m_components)
        {
            assert(component->m_rod == m_rod);
        }
    }

    void add_component(ForceBase* component)
    {
        if(m_components.empty())
        {
            m_rod = component->m_rod;
            m_components.push_back(component);
        }
        else
        {
            component->m_rod = m_rod;
            m_components.push_back(component);
        }
    }

    void set_rod(RodState const* rod) override
    {
        m_rod = rod;
        for (ForceBase* component : m_components)
        {
            component->set_rod(rod);
        }
    }

    double get_potential_energy() const override
    {
        double energy = 0;
        for (ForceBase const* component : m_components)
        {
            energy += component->get_potential_energy();
        }

        return energy;
    }

    VecXd get_force_vector() const override
    {
        VecXd force(m_rod->get_rest_dofs()->size());

        for (ForceBase const* component : m_components)
        {
            force += component->get_force_vector();
        }

        return force;
    }

    BandLimitedMatXd get_jacobian() const override
    {
        BandLimitedMatXd jacobian(m_rod->get_rest_dofs()->size(),
                m_rod->get_rest_dofs()->size());

        for (ForceBase const* component : m_components)
        {
            jacobian += component->get_jacobian();
        }

        return jacobian;
    }

private:
    std::vector<ForceBase*> m_components;
};

}

}
