/*
 * stretching.h
 */

#pragma once

#include "force_base.h"

namespace rodeo
{

namespace force
{

class Stretching: public ForceBase
{
public:
    Stretching(RodState const& rod, double k) :
            ForceBase(rod), m_k(k)
    {
    }

    virtual ~Stretching()
    {
    }

    double get_energy() const
    override
    {
        double energy = 0;
        for (size_t i = 0; i < m_rod.num_vertices() - 1; ++i)
        {
            Vec3d const edge = m_rod.get_vertex(i + 1) - m_rod.get_vertex(i);
            energy += 0.5 * m_k
                    * sqr(edge.norm() / m_rod.get_rest_length(i) - 1.);
        }

        return energy;
    }

    VecXd get_force() const override
    {
        VecXd force(m_rod.num_dofs());

        for (size_t i = 0; i < m_rod.num_vertices() - 1; ++i)
        {
            Vec3d const edge = m_rod.get_vertex(i + 1) - m_rod.get_vertex(i);
            force.segment < 3 > (i) += m_k
                    * (edge.norm() / m_rod.get_rest_length(i) - 1.)
                    * edge.normalized();
        }

        return force;
    }

    BandLimitedMatXd get_jacobian() const override
    {
        BandLimitedMatXd jacobian(m_rod.get_rest_dofs()->size(),
                m_rod.get_rest_dofs()->size());

        return jacobian; // FIXME
    }

private:

    double m_k;

};

}

}
