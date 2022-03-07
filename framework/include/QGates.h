#ifndef QGATES_H
#define QGATES_H

#include "QLib.h"


// -----------------------------------------------------------------
// --- Gates ---
// -----------------------------------------------------------------

class Gate__
{
    public:
        Gate__(YCS name) : name_orig_(name){name_ = name_orig_;}
        ~Gate__(){}

        /**
         * @brief Copy Construct.
         */
        Gate__(const Gate__& oo)
        {
            // copy to this object:

            name_ = oo.name_;
            type_ = oo.type_;
            name_orig_ = oo.name_orig_;
            ts_ = YVIv(oo.ts_);
            cs_ = YVIv(oo.cs_);
            conds_ = YVIv(oo.conds_);
            cs_add_ = YVIv(oo.cs_add_);

            copy_matrix2(oo.u2_, this->u2_);
            if(oo.un_a_ != nullptr)
                un_a_ = std::make_shared<YMATH::YMatrix>(oo.un_a_);
            if(oo.un_b_ != nullptr)
                un_b_ = std::make_shared<YMATH::YMatrix>(oo.un_b_);

            pars_ = YVQv(oo.pars_);

            add_inf_ = std::string(oo.add_inf_);
        }

        // copy this gate to a new gate:
        virtual YSG copy_gate() const { return std::make_shared<Gate__>(*this); }
        virtual void conjugate_transpose(){};
        virtual void generate(Qureg& oc){};
        virtual void write_to_file(YMIX::File& cf){}

        /**
         * @brief Modiffy a circuit matrix accordingly to the gate action.
         * @param Re real part of a circuit matrix.
         * @param Im imaginary part of a circuit matrix.
         */
        virtual void modify_circuit_matrix(YSM Re, YSM Im)
        { 
            // !!! to do a general case !!!
        };

        /**
         * @brief Add control qubits to the gate.
         * @param[in] reg_control control qubits to add.
         */
        inline
        void add_control_qubits(YCVI reg_control)
        {
            copy(reg_control.begin(), reg_control.end(), back_inserter(cs_));
        }

        inline
        void add_control_qubits(YCI c)
        {
            cs_.push_back(c);
        }

        inline
        std::string get_name(){ return name_; }

        /**
         * @brief E.g. the gate has target qubit [1] and control qubits [2,0].
         * If new_pos = [3,2,1], then the target qubit is shifted as 1 -> 2,
         * the control qubits are shifted as 2 -> 1, 0 -> 3.
         * @param[in] new_pos new positions of target/control/conditional qubits. 
         */
        void correct_qubits(YCVI new_pos);

        inline 
        std::string get_type(){ return type_; }

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

    public:
        const static std::string name_shared_;

    protected:
        std::string name_;// current name of the gate;
        std::string name_orig_;// original name of the gate;

        std::string type_;// type of the gate;

        YVIv ts_; // target qubits;
        YVIv cs_; // control qubits;
        YVIv conds_; // condition qubits;
        YVIv cs_add_; // additional control qubits to build multi-qubit unitary;

        ComplexMatrix2 u2_; // matrix for a single-target gate;
        YVQv pars_; // parameters of the gate (e.g. angles);

        YGlobalVariables gv_; // global variables;
        std::string add_inf_; // string line with additional information;

        YSM un_a_ = nullptr;
        YSM un_b_ = nullptr;
};

class GStop__ : public Gate__
{
public:
    GStop__(YCS name) : Gate__(name)
    {
        name_ = name;
        type_ = "stop";
    }
    YSG copy_gate() const { return std::make_shared<GStop__>(*this); }
    void write_to_file(YMIX::File& cf) override{ cf << "GStop " << name_ << "\n"; }
    void modify_circuit_matrix(YSM Re, YSM Im){}
};

class Box__ : public Gate__
{
    public:
        Box__(YCS name, YCVI ts, YCVI cs) : Gate__(name)
        {
            name_ = name;
            type_ = "box";
            ts_ = ts;
            cs_ = cs;
        }

        YSG copy_gate() const {return std::make_shared<Box__>(*this);}
        YSB copy_box() const {return std::make_shared<Box__>(*this);}

        void write_to_file(YMIX::File& cf) override
        {
            cf << "Box ";
            if(flag_start_) cf << "start ";
            else cf << "end ";
            cf << name_ << " " << ts_.size() << " ";
            for(auto& x: ts_)
                cf << x << " ";
            cf << "controls " << cs_.size() << " ";
            for(auto& x: cs_)
                cf << x << " ";
            cf << "\n";
        }

        void modify_circuit_matrix(YSM Re, YSM Im){}

        bool flag_start_ = true; // is it the left side of the box?
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

        void conjugate_transpose(){}

        void generate(Qureg& oc){}

        void write_to_file(YMIX::File& cf) override
        { 
            if(cs_.empty())
            {
                cf << name_ << " " << ts_[0] << "\n";
            } 
            else
            {
                cf << "mc" << name_ << " ";
                for(auto& c1: cs_)
                    cf << c1 << " ";
                cf << ts_[0] << "\n";
            } 
        }

        void modify_circuit_matrix(YMATH::YMatrix& M){}
};

class X__ : public SQGate__
{
    public:
        X__(YCI t) : SQGate__(name_shared_, t){ u2_ = gv_.mX; }

        YSG copy_gate() const { return std::make_shared<X__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                pauliX(oc, ts_[0]);
            else
                mc_st_u(oc, cs_, ts_[0], u2_); 
        }

        void modify_circuit_matrix(YSM Re, YSM Im){
            unsigned idt = ts_[0];
            unsigned N = Re->get_nr();
            unsigned nq = log2(N);

            YSM Re_copy = std::make_shared<YMATH::YMatrix>(Re);

            unsigned Nb = 1 << (nq - idt - 1);
            unsigned bl = 1 << (idt + 1);
            unsigned bl_half = bl/2;
            if(cs_.empty())
            {
                // --- modify a circuit matrix ---
                unsigned idb, idr_loc, idr_loc_half, idr_new;
                for(unsigned idr = 0; idr < N; ++idr)
                {
                    idb = idr / bl;
                    idr_loc = idr % bl;
                    idr_loc_half = idr % bl_half;
                    idr_new = idb * bl;

                    // std::cout << "idb = " << idb << std::endl;
                    // std::cout << "idr_loc = " << idr_loc << std::endl;
                    // std::cout << "idr_loc_half = " << idr_loc_half << std::endl;
                    if(idr_loc < bl_half) 
                        idr_new += bl_half + idr_loc_half;
                    else 
                        idr_new += idr_loc_half;
                    // std::cout << "idr_new = " << idr_new << std::endl;

                    for(unsigned idc = 0; idc < N; ++idc)
                        (*Re)(idr, idc) = (*Re_copy)(idr_new, idc);
                }
            }else
            {

            }
        }

    public:
        const static std::string name_shared_;
};

class Z__ : public SQGate__
{
    public:
        Z__(YCI t) : SQGate__(name_shared_, t){ u2_ = gv_.mZ; }

        YSG copy_gate() const { return std::make_shared<Z__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                pauliZ(oc, ts_[0]);
            else
                mc_st_u(oc, cs_, ts_[0], u2_); 
        }
    
    public:
        const static std::string name_shared_;
};

class H__ : public SQGate__
{
    public:
        H__(YCI t) : SQGate__(name_shared_, t){ u2_ = gv_.mH; }

        YSG copy_gate() const { return std::make_shared<H__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                hadamard(oc, ts_[0]); 
            else 
                mc_st_u(oc, cs_, ts_[0], u2_);
        }

        void modify_circuit_matrix(YSM Re, YSM Im){
            unsigned idt = ts_[0];
            unsigned N = Re->get_nr();
            unsigned nq = log2(N);
            unsigned sbs = 1 << idt;           // size of a small block
            unsigned Nb = 1 << (nq - idt - 1); // number of big blocks
            qreal coef_sqrt_2 = 1./sqrt(2);

            YSM Re_copy = std::make_shared<YMATH::YMatrix>(Re);
            if(cs_.empty())
            {
                unsigned idr_0, idn;

                // walk through all blocks on a diagonal
                for(unsigned ib = 0; ib < Nb; ++ib)
                {
                    idr_0 = (2*sbs) * ib;

                    // walk within a block
                    for(unsigned idr = idr_0; idr < (idr_0 + sbs); ++idr)
                    {
                        idn  = idr + sbs; 

                        // consider all columns
                        for(unsigned idc = 0; idc < N; ++idc)
                        {
                            (*Re)(idr, idc) = coef_sqrt_2 * (
                                (*Re_copy)(idr, idc) + (*Re_copy)(idn, idc)
                            );
                            (*Re)(idn, idc) = coef_sqrt_2 * (
                                (*Re_copy)(idr, idc) - (*Re_copy)(idn, idc)
                            );
                        }  
                    }
                }
            }else
            {

            }
        }

    public:
        const static std::string name_shared_;
};

class sR__ : public SQGate__
{
    public:
        sR__(YCS name, YCI t, YCQR a) : SQGate__(name, t)
        {   
            pars_.push_back(a);
        }

        YSG copy_gate() const { return std::make_shared<sR__>(*this); };

        void generate(Qureg& oc) 
        { 
            if(cs_.empty())
                unitary(oc, ts_[0], u2_); // take a general unitary function since the u2_ can be inversed;
            else 
                mc_st_u(oc, cs_, ts_[0], u2_);
        }

        void conjugate_transpose(){
            u2_ = YMATH::inv_matrix2(u2_);
            if(name_ == name_orig_)
                name_ = name_orig_ + "*";
            else
                name_ = name_orig_;
        }

        void write_to_file(YMIX::File& cf) override
        { 
            if(cs_.empty())
                cf << name_ << " " << ts_[0] << " " << pars_[0] << "\n"; 
            else
            {
                cf << "mc" + name_ + " ";
                for(auto& c1: cs_)
                    cf << c1 << " ";
                cf << ts_[0] << " " << pars_[0] << "\n";
            }
        }
};

class Ry__ : public sR__
{
    public:
        Ry__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = gv_.mRy(a); }
        YSG copy_gate() const { return std::make_shared<Ry__>(*this); };

    public:
        const static std::string name_shared_;
};

class Rz__ : public sR__
{
    public:
        Rz__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = gv_.mRz(a); }
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

        void write_to_file(YMIX::File& cf) override
        { 
            if(cs_.empty())
                cf << name_ << " " << ts_[0] << " " << pars_[0] << " " << pars_[1] << "\n"; 
            else
            {
                cf << "mc" + name_ + " ";
                for(auto& c1: cs_)
                    cf << c1 << " ";
                cf << ts_[0] << " " << pars_[0] << " " << pars_[1] << "\n";
            }
        }
};

class Rc__ : public sR2__
{
    // Ry(angle_ry) * Rz(angle_rz)
    public:
        Rc__(YCI t, YCQR angle_rz, YCQR angle_ry) : sR2__(name_shared_, t, angle_rz, angle_ry)
        { 
            u2_ = gv_.mRc(angle_rz, angle_ry); 
        }
        YSG copy_gate() const { return std::make_shared<Rc__>(*this); };

    public:
        const static std::string name_shared_;
};

class Phase__ : public sR__
{
    public:
        Phase__(YCI t, YCQR a) : sR__(name_shared_, t, a){ u2_ = gv_.mPhase(a); }
        YSG copy_gate() const { return std::make_shared<Phase__>(*this); };

    public:
        const static std::string name_shared_;
};

class condR_sep__ : public Gate__
{
    public:
        condR_sep__(
            YCS sp_name,
            YVIv reg_conds, YCI t, 
            YCVQ as
        );

        condR_sep__(const condR_sep__& oo) : Gate__(oo)
        {
            // caseI_   = oo.caseI_;  
            sp_name_ = oo.sp_name_;  
        }

        YSG copy_gate() const { return std::make_shared<condR_sep__>(*this); };

        void conjugate_transpose();

        void generate(Qureg& oc);

        void write_to_file(YMIX::File& cf) override
        { 
            cf << name_ << " ";

            cf << sp_name_ << " ";

            cf << "reg_conds " << size(conds_) << " ";
            for(auto& x: conds_) cf << x << " ";

            cf << "reg_control " << size(cs_) << " ";
            for(auto& x: cs_) cf << x << " ";

            cf << "target " << ts_[0] << " ";

            cf.of << "\n";
        }

        public:
            const static std::string name_shared_;

        private:
            std::string sp_name_;


};


#endif