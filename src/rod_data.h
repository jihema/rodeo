/*
 * rod_data.h
 *
 *  Created on: 9 May 2018
 */

#pragma once

#include "rod_state.h"

#include <boost/circular_buffer.hpp>

namespace rodeo
{

class RodData
{
public:

    static constexpr size_t epochs = 3; // 0 = past; 1 = present; 2 = future.

    RodData(VecXd const& rest_dofs) :
            m_rest_dofs(rest_dofs), //
            m_data(epochs)
    {
        for (size_t i = 0; i < epochs; ++i)
        {
            m_data.push_back(RodState(m_rest_dofs));
            m_data.back().set_dofs_to_rest();
            m_data.back().compute();
        }
    }

    void time_step()
    {
        assert(m_data.full());

        VecXd dofs = m_data.back().get_dofs(); // TODO: forward evolution come here.

        m_data.push_back(m_data.front()); // Reusing allocation.
        m_data.back().set_dofs(dofs, false); // Copy in evolved dofs.

        m_data.back().compute(); // On demand.
    }

    void back_step()
    {
        assert(m_data.full());

        VecXd dofs = m_data.front().get_dofs(); // TODO: backwards evolution come here.

        m_data.push_front(m_data.back());
        m_data.front().set_dofs(dofs, false);

        m_data.front().compute();
    }

    void print(std::ostream& os) const
    {
        os << "Past: " << m_data[0].get_energy() << "; ";
        os << "Present: " << m_data[1].get_energy() << "; ";
        os << "Future: " << m_data[2].get_energy() << ".";
    }

    const VecXd& get_rest_dofs() const
    {
        return m_rest_dofs;
    }

    void set_rest_dofs(const VecXd& restDofs)
    {
        m_rest_dofs = restDofs;
    }

private:

    VecXd m_rest_dofs;
    boost::circular_buffer<RodState> m_data;
};

inline std::ostream& operator<<(std::ostream& os, RodData const& rod_data)
{
    rod_data.print(os);
    return os;
}

}

