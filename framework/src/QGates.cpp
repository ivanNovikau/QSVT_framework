#include "../include/QGates.h"

using namespace std;

const string Gate__::name_shared_ = "gate";
const string X__::name_shared_ = "X";
const string Y__::name_shared_ = "Y";
const string Z__::name_shared_ = "Z";
const string H__::name_shared_ = "H";
const string Rx__::name_shared_ = "Rx";
const string Ry__::name_shared_ = "Ry";
const string Rz__::name_shared_ = "Rz";
const string Rc__::name_shared_ = "Rc";
const string Phase__::name_shared_ = "Phase";

void Gate__::correct_qubits(YCVI regs)
{
    int old_q = 0;
    auto ts = YVIv(ts_);
    auto cs = YVIv(cs_);
    auto conds = YVIv(conds_);
    for(auto& xx: regs)
    {
        change_element(ts,     ts_,     old_q, xx);
        change_element(cs,     cs_,     old_q, xx);
        change_element(conds,  conds_,  old_q, xx);
        ++old_q;
    } 
    ts_ = YVIv(ts);
    cs_ = YVIv(cs);
    conds_ = YVIv(conds);
}


void Gate__::write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq
){
    std::string l_nq_gate, l_name;

    // gate the most-signficant target qubit:
    auto id_top_q = get_most_signif_target_qubit();
    
    // width of the gate (take into account the number of target qubits):
    l_nq_gate = tex_gate_width(tex_lines, id_layer, nq, id_top_q);

    // write down the gate parameters:
    l_name = tex_get_gate_name(tex_lines, id_layer, nq);      

    // combine the gate information:
    tex_lines[nq - id_top_q - 1][id_layer] = "&\\gate" + l_nq_gate + "{" + l_name + "}";

    // add the control qubits 
    tex_add_control(tex_lines, id_layer, nq, id_top_q);
}

