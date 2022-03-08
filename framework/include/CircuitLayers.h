#ifndef CIRCUITLAYERS_H
#define CIRCUITLAYERS_H

#include "CLayer.h"

class CircuitLayers__
{
public:
    CircuitLayers__(YCU nq) : nq_(nq)
    {
        noc_layers_.resize(nq_, 0);
        layers_.push_back(std::make_shared<CLayer__>());
    }

    CircuitLayers__(const std::shared_ptr<const CircuitLayers__>& oo)
    {
        nq_ = oo->nq_;
        noc_layers_ = std::vector<uint64_t>(oo->noc_layers_);

        layers_.resize(oo->layers_.size());
        for(uint32_t id_layer = 0; id_layer < oo->layers_.size(); id_layer++)
        {
            layers_[id_layer] = std::make_shared<CLayer__>(oo->layers_[id_layer]);
        }
    }

    void add_gate(std::shared_ptr<Gate__>& gate)
    {
        YVIv ids_qubits_of_gate;

        // get all qubits, at which the gate acts (target, control etc.)
        gate->get_gubits_act_on(ids_qubits_of_gate);

        // find the first layer, which has enough free qubits to place the gate: 
        uint64_t id_min_noc_layer = std::numeric_limits<uint64_t>::max();
        int id_gate;
        for(uint32_t ii = 0; ii < ids_qubits_of_gate.size(); ii++)
        {
            id_gate = ids_qubits_of_gate[ii];

            if(id_min_noc_layer > noc_layers_[id_gate])
            {
                id_min_noc_layer = noc_layers_[id_gate];
            }
        }

        // insert the gate into the first non-occupied layer:
        gate->set_layer(id_min_noc_layer);
        layers_[id_min_noc_layer]->add_gate(gate);

        // shift the ids of non-occupied layers for the gate qubits:
        for(uint32_t ii = 0; ii < ids_qubits_of_gate.size(); ii++)
        {
            id_gate = ids_qubits_of_gate[ii];
            noc_layers_[id_gate] = id_min_noc_layer;
        }
    }





protected:
    // number of qubits in the circuit:
    uint32_t nq_;

    // noc_layers_[id_qubit] is id of the first non-occupied layer 
    //  at the qubit with id = id_qubit:
    std::vector<uint64_t> noc_layers_;

    // layers of the quantum circuit:
    std::vector<std::shared_ptr<CLayer__>> layers_;
};


#endif