#ifndef QGATES_H
#define QGATES_H

#include "QLib.h"


// -----------------------------------------------------------------
// --- Gates ---
// -----------------------------------------------------------------

class Gate__
{
    public:
        Gate__(YCS name) : name_(name)
        {
            id_layer_ = -1;
            flag_conj_ = false;
            tex_name_ = name_;
        }
        ~Gate__(){}

        /**
         * @brief Copy Construct.
         */
        Gate__(const Gate__& oo)
        {
            // copy to this object:

            name_ = oo.name_;
            type_ = oo.type_;
            tex_name_ = oo.tex_name_;
            
            ts_ = YVIv(oo.ts_);
            cs_ = YVIv(oo.cs_);
            conds_ = YVIv(oo.conds_);

            copy_matrix2(oo.u2_, this->u2_);
            if(oo.un_a_ != nullptr)
                un_a_ = std::make_shared<YMATH::YMatrix>(oo.un_a_);
            if(oo.un_b_ != nullptr)
                un_b_ = std::make_shared<YMATH::YMatrix>(oo.un_b_);

            pars_ = YVQv(oo.pars_);

            add_inf_ = std::string(oo.add_inf_);

            id_layer_ = oo.get_layer();

            flag_conj_ = oo.flag_conj_;
            flag_start_ = oo.flag_start_;
        }

        // copy this gate to a new gate:
        virtual YSG copy_gate() const { return std::make_shared<Gate__>(*this); }
        virtual void conjugate_transpose(){ flag_conj_ = !flag_conj_; };
        virtual void generate(Qureg& oc){};
        virtual void write_to_file(YMIX::File& cf){ write_to_file_base(cf); };

        /**
         * @brief Add control qubits to the gate.
         * @param[in] reg_control control qubits to add.
         */
        inline void add_control_qubits(YCVI reg_control)
        {
            copy(reg_control.begin(), reg_control.end(), back_inserter(cs_));
        }
        inline void add_control_qubits(YCI c)
        {
            cs_.push_back(c);
        }

        inline std::string get_name(){ return name_; }

        /**
         * @brief E.g. the gate has target qubit [1] and control qubits [2,0].
         * If new_pos = [3,2,1], then the target qubit is shifted as 1 -> 2,
         * the control qubits are shifted as 2 -> 1, 0 -> 3.
         * @param[in] new_pos new positions of target/control/conditional qubits. 
         */
        void correct_qubits(YCVI new_pos);

        inline 
        std::string get_type(){ return type_; }

        void get_gubits_act_on(YVI ids_qubit_act_on)
        {
            ids_qubit_act_on = YVIv(ts_);
            if(!cs_.empty())
                ids_qubit_act_on.insert(ids_qubit_act_on.end(), cs_.begin(), cs_.end());
            if(!conds_.empty())
                ids_qubit_act_on.insert(ids_qubit_act_on.end(), conds_.begin(), conds_.end());
        }

        void set_layer(const uint64_t& id_layer){ id_layer_ = id_layer; }

        uint64_t get_layer() const { return id_layer_; }

        inline void set_flag_start(YCB flag_start){ flag_start_ = flag_start; }
        inline bool get_flag_start(){ return flag_start_; }
        inline bool get_flag_conj(){ return flag_conj_; }

        void get_target_qubits(YVI ids_t){ ids_t = YVIv(ts_); }
        void get_control_qubits(YVI ids_c){ ids_c = YVIv(cs_); }

        virtual void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq
        );

    protected:
        inline
        void change_element(YVI v_new, YCVI v, YCI old_x, YCI new_x)
        {
            int index;
            auto it = find(v.begin(), v.end(), old_x);
            if (it != v.end())
            {
                index = it - v.begin();
                v_new[index] = new_x;        
            }    
        }

        inline
        void mc_st_u(Qureg& oc, YVI cs, YCI t, const ComplexMatrix2& u)
        {
            multiControlledUnitary(oc, &cs[0], cs.size(), t, u);
        }

        inline
        void copy_matrix2(const ComplexMatrix2& copy_from, ComplexMatrix2& copy_to)
        {
            copy_to.real[0][0] = copy_from.real[0][0]; 
            copy_to.real[0][1] = copy_from.real[0][1];
            copy_to.real[1][0] = copy_from.real[1][0];
            copy_to.real[1][1] = copy_from.real[1][1];

            copy_to.imag[0][0] = copy_from.imag[0][0]; 
            copy_to.imag[0][1] = copy_from.imag[0][1];
            copy_to.imag[1][0] = copy_from.imag[1][0];
            copy_to.imag[1][1] = copy_from.imag[1][1];
        }

        inline
        void write_to_file_base(YMIX::File& cf, const bool& flag_new_line=true)
        {
            cf << name_ << " ";
            
            if(flag_conj_)
                cf << "conj ";
            else
                cf << "orig "; 

            cf << "layer " << id_layer_ << " ";
            
            cf << "targets " << ts_.size() << " ";
            for(auto& x: ts_)
                cf << x << " ";

            cf << "controls " << cs_.size() << " ";
            for(auto& x: cs_)
                cf << x << " ";

            cf << "pars " << pars_.size() << " ";
            for(auto& x: pars_)
                cf << x << " ";
            
            if(flag_new_line)
                cf << "\n";
        }

        inline
        std::string tex_gate_width(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq,
            YCI id_top_q
        ){
            std::string l_nq_gate;
            if(ts_.size() == 1)
                l_nq_gate = "";
            else
                l_nq_gate = "[" + std::to_string(ts_.size()) + "]";
            return l_nq_gate;
        }

        inline
        std::string tex_get_gate_name(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq,
            YCS l_brackets = "()",
            YCB flag_inv_par = false
        ){
            std::string l_name, l_par;
            l_name = tex_name_;
            if(!flag_inv_par && flag_conj_)
                l_name += std::string("^\\dagger");

            if(pars_.size())
            {
                l_name.push_back(l_brackets[0]);
                for(auto id_p = 0; id_p < pars_.size(); id_p++)
                {
                    std::stringstream sstr;
                    if(flag_inv_par && flag_conj_)
                        sstr << -pars_[id_p];
                    else
                        sstr << pars_[id_p];
                    sstr >> l_par;
                    if(id_p > 0)
                        l_name += std::string(", ");
                    l_name += l_par;
                }
                l_name.push_back(l_brackets[1]);
            }
            return l_name;
        }

        inline
        void tex_add_control(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq,
            YCI id_top_q
        ){
            std::string l_c_dir;
            for(auto const& id_cq: cs_)
            {
                l_c_dir = std::to_string(id_cq - id_top_q);
                tex_lines[nq - id_cq - 1][id_layer] = "&\\ctrl{" + l_c_dir + "}";
            }
        }

        inline
        int get_most_signif_target_qubit(){ return *(std::max_element(ts_.begin(), ts_.end())); }

    public:
        const static std::string name_shared_;
        
    protected:
        std::string name_;// current name of the gate;
        std::string type_;// type of the gate;
        std::string tex_name_;// gate name as it is shown in the. tex file;

        YVIv ts_; // target qubits;
        YVIv cs_; // control qubits;
        YVIv conds_; // condition qubits;

        ComplexMatrix2 u2_; // matrix for a single-target gate;
        YVQv pars_; // parameters of the gate (e.g. angles);

        std::string add_inf_; // string line with additional information;

        YSM un_a_ = nullptr;
        YSM un_b_ = nullptr;

        uint64_t id_layer_; // id of the circuit layer, where the gate sits on;

        bool flag_conj_; // whether the gate is conjugated or not;
        bool flag_start_ = true; // is it the left side of the box? 
};

class GStop__ : public Gate__
{
public:
    GStop__(YCS name) : Gate__(name){ type_ = "stop"; }
    YSG copy_gate() const { return std::make_shared<GStop__>(*this); }
    void write_to_file(YMIX::File& cf){}
    void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq
    ){}
};

class Box__ : public Gate__
{
    public:
        Box__(YCS name, YCVI ts, YCVI cs, YCS tex_name="") : Gate__(name)
        {
            type_ = "box";
            ts_ = ts;
            cs_ = cs;

            if(!tex_name.empty())
                tex_name_ = tex_name;
        }

        YSG copy_gate() const {return std::make_shared<Box__>(*this);}
        YSB copy_box() const {return std::make_shared<Box__>(*this);}

        void write_to_file(YMIX::File& cf) override
        {
            cf << "Box ";
            if(flag_start_) 
                cf << "start ";
            else 
                cf << "end ";
            write_to_file_base(cf);
        }
};

class SQGate__ : public Gate__
{
    public:
        SQGate__(YCS name, YCI t) : Gate__(name)
        {   
            type_ = "q1";
            ts_ = YVIv(1);
            ts_[0] = t;
        }

        YSG copy_gate() const {return std::make_shared<SQGate__>(*this);};

        void conjugate_transpose(){
            flag_conj_ = !flag_conj_;
            u2_ = YMATH::inv_matrix2(u2_);
        }

        void generate(Qureg& oc){}
};

class X__ : public SQGate__
{
    public:
        X__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mX; }

        YSG copy_gate() const { return std::make_shared<X__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                pauliX(oc, ts_[0]);
            else
                mc_st_u(oc, cs_, ts_[0], u2_); 
        }

        void conjugate_transpose(){}

        void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq
        ){
            std::string l_nq_gate;

            // gate the most-signficant target qubit:
            auto id_top_q = get_most_signif_target_qubit();

            // combine the gate information:

            if(cs_.size() > 0)
                tex_lines[nq - id_top_q - 1][id_layer] = "&\\targ{}";
            else
                tex_lines[nq - id_top_q - 1][id_layer] = "&\\gate{X}";

            // add the control qubits 
            tex_add_control(tex_lines, id_layer, nq, id_top_q);
        }

    public:
        const static std::string name_shared_;
};

class Y__ : public SQGate__
{
    public:
        Y__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mY; }

        YSG copy_gate() const { return std::make_shared<Y__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                pauliY(oc, ts_[0]);
            else
                mc_st_u(oc, cs_, ts_[0], u2_); 
        }

        void conjugate_transpose(){}

    public:
        const static std::string name_shared_;
};

class Z__ : public SQGate__
{
    public:
        Z__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mZ; }

        YSG copy_gate() const { return std::make_shared<Z__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                pauliZ(oc, ts_[0]);
            else
                mc_st_u(oc, cs_, ts_[0], u2_); 
        }

        void conjugate_transpose(){}
    
    public:
        const static std::string name_shared_;
};

class H__ : public SQGate__
{
    public:
        H__(YCI t) : SQGate__(name_shared_, t){ u2_ = YGV::mH; }

        YSG copy_gate() const { return std::make_shared<H__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                hadamard(oc, ts_[0]); 
            else 
                mc_st_u(oc, cs_, ts_[0], u2_);
        }

        void conjugate_transpose(){}

    public:
        const static std::string name_shared_;
};

class sR__ : public SQGate__
{
    public:
        sR__(YCS name, YCI t, YCQR a) : SQGate__(name, t){ pars_.push_back(a); }
        YSG copy_gate() const { return std::make_shared<sR__>(*this); };
        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                unitary(oc, ts_[0], u2_); // take a general unitary function since the u2_ can be inversed;
            else 
                mc_st_u(oc, cs_, ts_[0], u2_);
        }
};

class Rx__ : public sR__
{
    public:
        Rx__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mRx(a); tex_name_ = "R_x"; }
        YSG copy_gate() const { return std::make_shared<Rx__>(*this); };

    public:
        const static std::string name_shared_;
};

class Ry__ : public sR__
{
    public:
        Ry__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mRy(a); tex_name_ = "R_y"; }
        YSG copy_gate() const { return std::make_shared<Ry__>(*this); };

    public:
        const static std::string name_shared_;
};

class Rz__ : public sR__
{
    public:
        Rz__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mRz(a); tex_name_ = "R_z";}
        YSG copy_gate() const { return std::make_shared<Rz__>(*this); };

    public:
        const static std::string name_shared_;
};

class sR2__ : public sR__
{
    public:
        sR2__(YCS name, YCI t, YCQR a1, YCQR a2) : sR__(name, t, a1)
        {   
            pars_.push_back(a2);
        }

        YSG copy_gate() const { return std::make_shared<sR2__>(*this); };
};

/** Ry(angle_ry) * Rz(angle_rz) */
class Rc__ : public sR2__
{
    public:
        Rc__(YCI t, YCQR angle_rz, YCQR angle_ry) : sR2__(name_shared_, t, angle_rz, angle_ry)
        { 
            u2_ = YGV::mRc(angle_rz, angle_ry); 
            tex_name_ = "R_c";
        }
        YSG copy_gate() const { return std::make_shared<Rc__>(*this); };

    public:
        const static std::string name_shared_;
};

class Phase__ : public sR__
{
    public:
        Phase__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = YGV::mPhase(a); tex_name_ = "\\phase"; }
        YSG copy_gate() const { return std::make_shared<Phase__>(*this); };

        void write_tex(
            std::vector<std::vector<std::string>>& tex_lines, 
            const uint64_t& id_layer,
            YCU nq
        ){
            std::string l_nq_gate, l_name;

            // gate the most-signficant target qubit:
            auto id_top_q = get_most_signif_target_qubit();

            // write down the gate parameters:
            l_name = tex_get_gate_name(tex_lines, id_layer, nq, "{}", true);

            // combine the gate information:
            tex_lines[nq - id_top_q - 1][id_layer] = "&" + l_name;

            // add the control qubits 
            tex_add_control(tex_lines, id_layer, nq, id_top_q);
        }

    public:
        const static std::string name_shared_;
};



#endif