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

void RodState::compute()
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
