#ifndef CLAYER_H
#define CLAYER_H

#include "QGates.h"

class CLayer__
{
public:
    CLayer__()
    {

    }

    CLayer__(const std::shared_ptr<const CLayer__>& oo)
    {
        gates_.resize(oo->gates_.size());
        for(uint64_t id_gate = 0; id_gate < oo->gates_.size(); id_gate++)
        {
            gates_[id_gate] = oo->gates_[id_gate]->copy_gate();
        }
    }

    void add_gate(std::shared_ptr<Gate__>& gate)
    {
        gates_.push_back(gate);
    }


protected:
    // gates in the circuit;
    std::vector<YSG> gates_; 

};



#endif