#ifndef CIRCUITLAYERS_H
#define CIRCUITLAYERS_H

#include "CLayer.h"

class CircuitLayers__
{
public:
    CircuitLayers__(YCU nq) : nq_(nq)
    {
        noc_layers_.resize(nq_, 0);

        // create the zero-th layer:
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
        uint64_t id_first_noc_layer = 0;
        int id_qubit;
        for(auto& id_qubit: ids_qubits_of_gate)
            if(id_first_noc_layer < noc_layers_[id_qubit])
                id_first_noc_layer = noc_layers_[id_qubit];

        // insert the gate into the first non-occupied layer:
        gate->set_layer(id_first_noc_layer);
        layers_[id_first_noc_layer]->add_gate(gate);

        // shift the ids of non-occupied layers for the gate qubits:
        uint64_t id_new_noc_layer = id_first_noc_layer + 1;
        for(auto& id_qubit: ids_qubits_of_gate)
        {
            noc_layers_[id_qubit] = id_new_noc_layer;
        }

        // create the next non-occupied layer if needed:
        if(id_new_noc_layer >= layers_.size())
            layers_.push_back(std::make_shared<CLayer__>());
    }

    /**
     * @brief Return the number of non-empty layers
     */
    uint64_t get_n_layers(){ 
        if(layers_.back()->is_empty())
            return layers_.size() - 1; 
        else
            return layers_.size();
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