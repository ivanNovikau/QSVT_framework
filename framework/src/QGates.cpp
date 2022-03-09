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

