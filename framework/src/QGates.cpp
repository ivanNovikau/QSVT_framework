#include "../include/QGates.h"

using namespace std;

const string Gate__::name_shared_ = "gate";
const string X__::name_shared_ = "X";
const string Z__::name_shared_ = "Z";
const string H__::name_shared_ = "H";
const string Ry__::name_shared_ = "Ry";
const string Rz__::name_shared_ = "Rz";
const string Rc__::name_shared_ = "Rc";
const string Phase__::name_shared_ = "Phase";
const string condR_sep__::name_shared_ = "cond_R";


void Gate__::correct_qubits(YCVI regs)
{
    int old_q = 0;
    auto ts = YVIv(ts_);
    auto cs = YVIv(cs_);
    auto cs_add = YVIv(cs_add_);
    auto conds = YVIv(conds_);
    for(auto& xx: regs)
    {
        change_element(ts,     ts_,     old_q, xx);
        change_element(cs,     cs_,     old_q, xx);
        change_element(cs_add, cs_add_, old_q, xx);
        change_element(conds,  conds_,  old_q, xx);
        ++old_q;
    } 
    ts_ = YVIv(ts);
    cs_ = YVIv(cs);
    cs_add_ = YVIv(cs_add);
    conds_ = YVIv(conds);
}

condR_sep__::condR_sep__(
    YCS sp_name,
    YVIv reg_conds, YCI t, 
    YCVQ as
) : Gate__(condR_sep__::name_shared_)
{   
    sp_name_ = sp_name;
    type_ = "qm";
    conds_ = YVIv(reg_conds);
    ts_.push_back(t);
    pars_ = YVQv(as);

    // --- save matrix ---
    unsigned nvv = reg_conds.size();
    unsigned Nv = pow(2, nvv);
    unsigned N = Nv * 2;
    int nq = nvv + 1;

    un_a_ = std::make_shared<YMATH::YMatrix>(N, N);
    un_b_ = std::make_shared<YMATH::YMatrix>(N, N);

    // create identity matrix
    for(unsigned i = 0; i < N; i++)
        (*un_a_)(i,i) = 1.0;

    unsigned sh;
    qreal a1;
    for(unsigned i = 0; i < Nv; i++)
    {
        a1 = 0.5*pars_[i];
        sh = i*2; // double since we fill 2x2 matrices on a diagonal
        
        // this is a matrix of Ry(pars_[i]):
        (*un_a_)(sh,sh)   = cos(a1); (*un_a_)(sh,sh+1)   = -sin(a1);
        (*un_a_)(sh+1,sh) = sin(a1); (*un_a_)(sh+1,sh+1) = cos(a1);
    }
}

void condR_sep__::conjugate_transpose()
{
    unsigned nv = conds_.size();
    unsigned Nv = pow(2, nv);
    int nq = nv + 1;

    unsigned sh;
    qreal temp;
    for(unsigned i = 0; i < Nv; i++)
    {
        sh = i*2;

        temp = (*un_a_)(sh,sh+1);
        (*un_a_)(sh,sh+1) = (*un_a_)(sh+1,sh);
        (*un_a_)(sh+1,sh) = temp; 

        temp = -(*un_b_)(sh,sh+1);
        (*un_b_)(sh,sh)   = -(*un_b_)(sh,sh);   (*un_b_)(sh,sh+1)   = -(*un_b_)(sh+1,sh);
        (*un_b_)(sh+1,sh) = temp;               (*un_b_)(sh+1,sh+1) = -(*un_b_)(sh+1,sh+1);
    }

    if(name_ == condR_sep__::name_shared_)
    {
        name_ = condR_sep__::name_shared_ + "*";
        sp_name_ = sp_name_ + "*";
    }
    else
    {
        name_ = condR_sep__::name_shared_;
        sp_name_ = sp_name_.substr(0, sp_name_.size()-1);
    }
}

void condR_sep__::generate(Qureg& oc) 
{ 
    unsigned nv = conds_.size();
    int nq = nv + 1;

    // --- resulting unitary matrix ---
    ComplexMatrixN U;
    U.numQubits = nq;
    U.real = un_a_->get_pointer(); U.imag = un_b_->get_pointer();

    // According to the definition of the second argument of the function multiQubitUnitary, 
    //      the qubits must be placed from the least to the most significant ones.
    // Moreover, due to the construction of the matrix above,
    //      we must assume here locally that the target qubit is the least significant.
    YVIv targs;
    targs.push_back(ts_[0]);
    copy(conds_.begin(), conds_.end(), back_inserter(targs));
    if(cs_.size() == 0)
        multiQubitUnitary(oc, &targs[0], nq, U);
    else
        multiControlledMultiQubitUnitary(oc, &cs_[0], cs_.size(), &targs[0], nq, U); 
}