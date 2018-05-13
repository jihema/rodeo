/**
 * gravitation.h
 */

#include "force_base.h"

namespace rodeo
{

namespace force
{

class Gravitation: public ForceBase
{
public:

    Gravitation(Vec3d const& g = Vec3d(0., 0., -9.81)) :
            m_g(g)
    {
    }

    virtual ~Gravitation()
    {
    }

    double get_potential_energy() const override
    {
        double energy = 0;
        for (size_t i = 0; i < m_rod->num_vertices(); ++i)
        {
            energy -= m_rod->get_vertex_mass(i) * m_rod->get_vertex(i).dot(m_g);
        }

        return energy;
    }

    VecXd get_force_vector() const override
    {
        VecXd force(m_rod->get_rest_dofs()->size());

        for (size_t i = 0; i < m_rod->num_vertices(); ++i)
        {
            force.segment<3>(4 * i) = m_rod->get_vertex_mass(i) * m_g;
        }

        return force;
    }

    BandLimitedMatXd get_jacobian() const override
    {
        BandLimitedMatXd jacobian(m_rod->get_rest_dofs()->size(),
                m_rod->get_rest_dofs()->size());

        return jacobian; // Zero.
    }

private:

    const Vec3d m_g;
};

}

}
