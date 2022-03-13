#ifndef QLIB_H
#define QLIB_H

#include "QuEST.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include <filesystem>

#include <math.h>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

#include <stdlib.h> 

#include <limits>
// #include <time.h>
#include <chrono>
#include <stdarg.h>
#include <memory>
#include <numeric>
#include <unistd.h>

#include <armadillo> // to remove

#include <mpi.h> // to remove

#include <cuda_runtime_api.h>

#include "H5Cpp.h"

#define YMPI  false
#define YCUDA true

// ------------------------------------------
// --- Type synonyms --- 
// ------------------------------------------
#define YCB const bool&

#define YCsh  const short&
#define YCVsh const std::vector<short>&
#define YVsh  std::vector<short>&
#define YVshv std::vector<short>

#define YCU  const uint32_t&
#define YCVU const std::vector<uint32_t>&
#define YVU  std::vector<uint32_t>&
#define YVUv std::vector<uint32_t>

#define YCI  const int&
#define YCVI const std::vector<int>&
#define YVI  std::vector<int>&
#define YVIv std::vector<int>

#define YCVT const std::vector<T>&

#define YCS  const std::string&
#define YS   std::string&
#define YCVS const std::vector<std::string>&
#define YVSv std::vector<std::string>

#define YCQR const qreal&
#define YCVQ const std::vector<qreal>&
#define YVQ  std::vector<qreal>&
#define YVQv std::vector<qreal>

#define YCVCo const std::vector<Complex>&
#define YVCo  std::vector<Complex>&
#define YVCov  std::vector<Complex>


#define YPQC YSQ 
#define YISS std::istringstream& 

#define YQCP QCircuit*
#define YSQ  std::shared_ptr<QCircuit>
#define YSG  std::shared_ptr<Gate__>
#define YSB  std::shared_ptr<Box__>
#define YSM  std::shared_ptr<YMATH::YMatrix>
#define YSCQ  std::shared_ptr<const QCircuit>
#define YCCQ const std::shared_ptr<const QCircuit>
#define YCCB const std::shared_ptr<const Box__>
#define YCCM const std::shared_ptr<const YMATH::YMatrix>

#define YCVVI const std::vector<std::vector<int>>&
#define YVVI std::vector<std::vector<int>>&
#define YVVIv std::vector<std::vector<int>>


#define YGV YGlobalVariables


#define LOG_INDENT "    "
#define WHITESPACE " \n\r\t\f\v"

#define YSIZE_CHAR_ARRAY 100000

#define TEX_BEGIN_CIRCUIT    "\\begin{quantikz}[row sep={0.5cm,between origins},column sep={0.3cm}]\n"
#define TEX_VERTICAL_GAP     "\\vspace{20pt}"
#define TEX_CIRCUIT_LENGTH 10

#define COMMENT_ORACLE    "//"
#define FORMAT_ORACLETOOL ".oracle"s 
#define FORMAT_PROFILE    ".condR_profile"s
#define FORMAT_CIRCUIT    ".circuit"s
#define FORMAT_TEX        ".tex"s
#define FORMAT_ANGLES     ".angles"s
#define FORMAT_LOG        ".log"s
#define FORMAT_QSP        ".qsp"s
#define FORMAT_INIT       ".init_state"s
#define ENDING_FORMAT_OUTPUT "_OUTPUT.hdf5"s
#define ENDING_FORMAT_RESTART "_RESTART.hdf5"s

// ------------------------------------------
// --- Structure with an initial state --- 
// ------------------------------------------
struct INIT_STATE__
{
    bool flag_defined = false;
    long long b_ampl; // position of the first element in a state vector to set
    long long n_ampls; // number of elements to set
    YVQv ampl_vec_real;
    YVQv ampl_vec_imag;
};

// ------------------------------------------
// --- Structure with Global variables --- 
// ------------------------------------------
struct YGlobalVariables{

    std::string reg_whole_circuit = "the_whole_circuit";

    static int tex_circuit_length; // length of the circuit in the .tex file;

    qreal inv_sqrt2 = 1./sqrt(2);

    const ComplexMatrix2 mX = 
    {
        .real = {
            {0., 1.}, 
            {1., 0.}
        },
        .imag = {{0., 0.}, {0., 0.}}
    };

    const ComplexMatrix2 mY = 
    {
        .real = {
            {0., 0.}, 
            {0., 0.}
        },
        .imag = {
            {0., -1.}, 
            {1.,  0.}
        }
    };

    const ComplexMatrix2 mZ = 
    {
        .real = {
            {1.,  0.}, 
            {0., -1.}
        },
        .imag = {{0., 0.}, {0., 0.}}
    };

    const ComplexMatrix2 mH = 
    {
        .real = {
            {inv_sqrt2,  inv_sqrt2}, 
            {inv_sqrt2, -inv_sqrt2}
        },
        .imag = {{0., 0.}, {0., 0.}}
    };

    ComplexMatrix2 mRx(YCQR a)
    {
        qreal a2 = a/2.;
        ComplexMatrix2 res = {
            .real = {
                {cos(a2),      0.},
                {     0., cos(a2)}
            },
            .imag = {
                {      0., -sin(a2)}, 
                {-sin(a2),       0.}
            }
        };
        return res;
    }

    ComplexMatrix2 mRy(YCQR a)
    {
        qreal a2 = a/2.;
        ComplexMatrix2 res = {
            .real = {
                {cos(a2), -sin(a2)},
                {sin(a2),  cos(a2)}
            },
            .imag = {{0., 0.}, {0., 0.}}
        };
        return res;
    }

    ComplexMatrix2 mRz(YCQR a)
    {
        qreal a2 = a/2.;
        ComplexMatrix2 res = {
            .real = {
                {cos(a2),      0.},
                {0.,      cos(a2)}
            },
            .imag = {
                {-sin(a2),      0.}, 
                {0.,       sin(a2)}
            }
        };
        return res;
    }

    ComplexMatrix2 mRc(YCQR az, YCQR ay)
    {
        qreal az2 = az/2.;
        qreal ay2 = ay/2.;

        // Ry(ay) * Rz(az)
        ComplexMatrix2 res = {
            .real = {
                {cos(az2)*cos(ay2), -cos(az2)*sin(ay2)},
                {cos(az2)*sin(ay2),  cos(az2)*cos(ay2)}
            },
            .imag = {
                {-sin(az2)*cos(ay2),  -sin(az2)*sin(ay2)}, 
                {-sin(az2)*sin(ay2),  sin(az2)*cos(ay2)}
            }
        };
        return res;
    }

    ComplexMatrix2 mPhase(YCQR a)
    {
        ComplexMatrix2 res = {
            .real = {
                {1., 0.},
                {0., cos(a)}
            },
            .imag = {
                {0., 0.}, 
                {0., sin(a)}
            }
        };
        return res;
    }


};

// ------------------------------------------
// --- Math functions --- 
// ------------------------------------------ 
namespace YMATH{
    /** Check if the variable \p x is zero.*/
    bool is_zero(YCQR x);

    /** Convert an integer \p x to an array of bits \p binaryNum :
    * 2 -> {1, 0}; 
    * 6 -> {1, 1, 0} etc.
    * @param[out] binaryNum resulting array of bits. Its size must have been initialized.
    */
    void intToBinary(int x, std::vector<short>& binaryNum);

    /**
     * @brief Convert an array of bits to a corresponding integer.
     * @param bb array of binary. bb[0] - the most significant bit.
     * @return integer = bb[0]*2**(n - 1) + ... + bb[n-2]*2**1 + bb[n-1]
     */
    long long int binaryToInt(const std::vector<short>& bb);

    /**
     * @brief Two-dimensional matrix of a qreal (which is usually double) type.
     * It is an object that provides some basic manipulations with a two-dimensional 
     * raw pointer.
     */
    class YMatrix{
        public:
            /**
             * @brief Create an empty matrix object whithout reserving any memory.
             */
            YMatrix();

            /**
             * @brief Create a zero matrix (\p Nrows, \p Ncols).
             * @param Nrows is a number of rows in a new matrix.
             * @param Ncols is a number of columns.
             * @param flag_zero if true, set elements to zero.
             */
            YMatrix(YCU Nrows, YCU Ncols, YCB flag_zero=true);

            /**
             * @brief Copy a matrix \p M. 
             */
            YMatrix(YCCM M);
            ~YMatrix();

            /**
             * @brief Create an empty matrix.
             */
            void create(YCU Nrows, YCU Ncols);

            /**
             * @brief Create a zero matrix.
             */
            void create_zero_matrix(YCU Nrows, YCU Ncols);

            /**
             * @brief Create an identity matrix.
             */
            void create_identity_matrix(YCU Nrows, YCU Ncols);

            void create_x10_matrix(YCU Nrows, YCU Ncols);

            /**
             * @brief Gives a raw pointer of the matrix.
             */
            qreal** get_pointer();

            /**
             * @brief Create and return 1-D pointer to the matrix
             */
            qreal* get_1d_pointer();

            /**
             * @brief Get a pointer to a transposed matrix
             */
            qreal* get_1d_transposed_pointer();

            inline int get_nr(){return nr_;}
            inline int get_nc(){return nc_;}

            void set_squared_from_transposed_1d_matrix(int N, qreal* M);

            inline
            qreal& operator()(YCU id_row, YCU id_col)
            {
                if (id_row >= nr_)
                {
                    std::cerr << "\nError: id-row = " << id_row << ", while n-rows = " << nr_ << std::endl;
                    exit(-1);
                }
                if (id_col >= nc_)
                {
                    std::cerr << "\nError: id-column = " << id_col << ", while n-columns = " << nc_ << std::endl;
                    exit(-1);
                }
                return a_[id_row][id_col];
            }

            inline
            qreal operator()(YCU id_row, YCU id_col) const
            {
                if (id_row >= nr_)
                {
                    std::cerr << "\nError: id-row = " << id_row << ", while n-rows = " << nr_ << std::endl;
                    exit(-1);
                }
                if (id_col >= nc_)
                {
                    std::cerr << "\nError: id-column = " << id_col << ", while n-columns = " << nc_ << std::endl;
                    exit(-1);
                }
                return a_[id_row][id_col];
            }

            /**
             * @brief Print the matrix with precision = \p prec.
             * @param prec precision of matrix elements.
             * @param flag_scientific if true, then print in the scientific notation.
             * @param wc width of every column.
             */
            void print(int prec=3, bool flag_scientific=false, int wc=2);

        protected:
            /**
             * @brief Reserve memory for a new matrix of a known size.
             * The function does not check whether the matrix has been already initialized.
             */
            void create_new();

            /**
             * @brief Free the memory occupied by the matrix.
             */
            void clear();

            /**
             * @brief Set elements to zeros.
             * The function does not check where the matrix has been already initialized.
             */
            void set_zeros();

        private:
            unsigned nr_, nc_;
            qreal** a_ = nullptr;
            std::shared_ptr<qreal[]> a_1d_ = nullptr;
            std::shared_ptr<qreal[]> a_1d_transposed_ = nullptr;

    };

    ComplexMatrix2 inv_matrix2(const ComplexMatrix2& a);

    /**
     * @brief Create a vector with integers in the interval [start, end).
     */
    YVIv get_range(YCI start, YCI end);

    /**
     * @brief Whether a string \p str can be converted to a number.
     * @param str is a string to check;
     */
    bool is_number(YCS str);
}

// ------------------------------------------
// --- Mix functions --- 
// ------------------------------------------
namespace YMIX{
    void print_log(
        const QuESTEnv& env, 
        YCS line, 
        YCI n_indent = 0, 
        const bool& flag_only_file=false, 
        const bool& flag_new_line=true
    );
    void print_log_flush(const QuESTEnv& env, YCS line, YCI n_indent=0);
    void print_log_err(const QuESTEnv& env, YCS line);

    /** Timer*/
    struct YTimer{
        public:           
            void Start(){
                start_ = std::chrono::steady_clock::now();
            }
            void Stop(){
                end_ = std::chrono::steady_clock::now();
            }
            void StartPrint(const QuESTEnv& env, YCS mess)
            {
                if(env.rank == 0) 
                {
                    Start();
                    print_log_flush(env, mess);
                }
            }
            void StopPrint(const QuESTEnv& env)
            {
                if(env.rank == 0)
                {
                    Stop();

                    std::ostringstream oss;
                    oss << std::scientific << std::setprecision(3) 
                        << get_dur_s() << " s\n";
                    print_log_flush(env, oss.str());
                }
            }
            double get_dur(){
                std::chrono::duration<double> dur_seconds = end_ - start_;
                return 1000.*dur_seconds.count(); // in ms
            }
            double get_dur_s(){
                std::chrono::duration<double> dur_seconds = end_ - start_;
                return dur_seconds.count(); // in seconds
            }
            std::string get_dur_str_ms(){
                std::ostringstream ostr;
                ostr << std::scientific << std::setprecision(3) << get_dur() << " ms";
                return ostr.str();
            }
            std::string get_dur_str_s(){
                std::ostringstream ostr;
                ostr << std::scientific << std::setprecision(3) << get_dur()/1000. << " s";
                return ostr.str();
            }
        protected:
            std::chrono::time_point<std::chrono::steady_clock> start_;
            std::chrono::time_point<std::chrono::steady_clock> end_;
    };

    void get_state_vector(const Qureg& qq, YCU nq, YVQ state_real, YVQ state_imag);

    /** Form a string with a circuit wavefunction
    * @param[out] str_wv output string with a full state vector;
    * @param[in]  ampls  input amplitudes of different \p states;
    * @param[in]  states input states;
    * @param[in]  organize_state input array that indicates how to output every state
    *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
    *        amplitude  |i3 i2> |i1> |i0>;
    * @param[in] ampl_prec precision of amplitudes to output in \p str_wv;
    */
    void getStrWavefunction(
        YS str_wv, 
        YCVCo ampls, 
        const std::list<std::vector<short>>& states, 
        YCVU organize_state, 
        YCU ampl_prec = 3 
    );

    /** Compute all states of a given circuit and their amplitudes 
    * @param[in] qq circuit;
    * @param[in] organize_state input array that indicates how to output every state
    *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
    *        amplitude  |i3 i2> |i1> |i0>;
    * @param[in] ampl_prec input precision of amplitude to output in \p str_wv;
    * @param[out] str_wv: output string with a full state vector;
    * @param[out] states: output resulting states;
    * @param[out] ampls:  output amplitudes of different \p states;
    */
    void Wavefunction(
        const Qureg& qq, 
        YS str_wv, 
        std::list<std::vector<short>>& states, 
        YVCo ampls, 
        YCVU organize_state,
        YCU ampl_prec = 3
    );

    /** Compute first \p n_states states starting from the low-priority qubits.
     * Returns only states with non-zero probability. 
    * @param[in] qq circuit;
    * @param[in] n_states a number of first low-priority to calculate.
    * @param[in] organize_state input array that indicates how to output every state
    *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
    *        amplitude  |i3 i2> |i1> |i0>;
    * @param[out] str_wv output string with a full state vector;
    * @param[out] states output resulting states;
    * @param[out] ampls  output amplitudes of different \p states;
    * @param[in] state_to_choose array with bits that a state must have to be chosen:
    *   Elements of this vector must be 0 or 1. 
    *   -1 means that it can be either 0 or 1.
    *   The first elements in the array correspond to the high-priority qubits.
    *   The size of the array must take into account the ancillae.
    * @param[in] prec input precision of amplitude to output in \p str_wv;
    */
    void Wavefunction_NonzeroProbability(
        const Qureg& qq, 
        YCU n_states,
        YCVU organize_state,
        YS str_wv, 
        std::list<std::vector<short>>& states, 
        YVCo ampls, 
        YCVsh state_to_choose = YVshv{},
        YCU prec = 3
    );

    /** Get only states with non-zero amplitudes 
    * @param[in] states_init initial circuit states;
    * @param[in] ampls_init  initial state amplitudes;
    * @param[in] organize_state input array that indicates how to output every state
    *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
    *        amplitude  |i3 i2> |i1> |i0>;
    * @param[in] ampl_prec input precision of amplitude to output in \p str_wv;
    * @param[out] str_wv output string with states with only non-zero ampltides;
    * @param[out] states output states with only non-zero ampltides;
    * @param[out] ampls  output non-zero amplitudes of different \p states;
    */
    void getNonzeroWavefunction(
        const std::list<std::vector<short>>& states_init, 
        const std::vector<Complex>& ampls_init, 
        const std::vector<unsigned>& organize_state, 
        std::string& str_wv, 
        std::list<std::vector<short>>& states, 
        std::vector<Complex>& ampls, 
        const unsigned ampl_prec = 3 
    );

    /** Get special states (but only with non-zero amplitudes).
    * @param[in] states_to_choose array with bits that \p states_init must have to be chosen:
    *   E.g. states_init = |x2 x1 x0>. 
    *   If states_to_choose = {0},   choose only states with x2 = 0.
    *   If states_to_choose = {0,1}, choose only states with x2 = 0, x1 = 1.
    *   If states_to_choose = {0,-1, 1}, choose only states with x2 = 0, x0 = 1 and any x1.
    * @param[in] state_init initial circuit states;
    * @param[in] ampls_init  initial state amplitudes;
    * @param[in] organize_state input array that indicates how to output every state
    *        e.g. if n_qubits = 4, and organize_state = [2,1,1], then
    *        amplitude  |i3 i2> |i1> |i0>;
    * @param[in] ampl_prec input precision of amplitude to output in \p str_wv;
    * @param[out] str_wv output string with chosen states;
    * @param[out] states chosen states;
    * @param[out] ampls  amplitudes of the chosen states \p states;
    */
    void getSpecialStates(
        const std::vector<short>& state_to_choose,
        const std::list<std::vector<short>>& states_init, 
        const std::vector<Complex>& ampls_init, 
        const std::vector<unsigned>& organize_state, 
        std::string& str_wv, 
        std::list<std::vector<short>>& states, 
        std::vector<Complex>& ampls, 
        const unsigned ampl_prec = 3 
    );


    struct File
    {
        std::ofstream of; // stream connected to a file;
        
        /** Open a file.
         * @param[in] fileName name of the file;
         * @param[in] flagNew  is it a new file (if new, old data will be deleted);
         * */
        File(YCS fileName, bool flagNew=false)
        {
            if(fileName.empty())
            {
                std::cerr << "Error: File name is not set." << std::endl;
                exit(-1);
            }

            if(flagNew)
                of.open(fileName);
            else
                of.open(fileName, std::ios::app);
            if(!of)
            {
                std::cerr << "Error: It's impossible to open the file:\n" << fileName << std::endl;
                exit(-1);
            }
        }

        ~File()
        {
            if(of.is_open())
            {
                of.clear();
                of.close();
            }
        }

        template<class T>
        std::ostream& operator<<(const T& x)
        {
            of << x;
            return of;
        }
    };


    struct LogFile : File
    {
        static std::string name_global_;
        LogFile(bool flagNew=false) : File(name_global_, flagNew){}
    };

    /**
     * @brief Remove a comment from a line.
     * @param line line from where to remove a comment;
     * @return a new line;
     */
    std::string remove_comment(YCS line);

    // true if line1 == line2
    bool compare_strings(YCS line1, YCS line2);

    // true if line1 == line2 and in lines
    bool compare_strings(YCS line1, YCS line2, YCVS lines);

    // true if line1 in lines
    bool compare_strings(YCS line1, YCVS lines);

    template<class T>
    std::vector<T> conc_vectors(YCVT v1, YCVT v2)
    {
        std::vector<T> vv(v1);
        copy(v2.begin(), v2.end(), back_inserter(vv));

        return vv;
    }
    template<class T>
    std::vector<T> conc_vectors(YCVT v1, YCVT v2, YCVT v3)
    {
        std::vector<T> vv(v1);
        copy(v2.begin(), v2.end(), back_inserter(vv));
        copy(v3.begin(), v3.end(), back_inserter(vv));
        return vv;
    }
    template<class T>
    std::vector<T> conc_vectors(YCVT v1, YCVT v2, YCVT v3, YCVT v4)
    {
        std::vector<T> vv(v1);
        copy(v2.begin(), v2.end(), back_inserter(vv));
        copy(v3.begin(), v3.end(), back_inserter(vv));
        copy(v4.begin(), v4.end(), back_inserter(vv));
        return vv;
    }

    void print(const ComplexMatrix2& a, YCI prec = 3);

    template<class T>
    void print(const std::vector<T>& a, YCI prec = 3, const bool& flag_scientific = false)
    {
        for(auto& x: a)
            if(flag_scientific)
                std::cout << std::scientific << std::setprecision(prec) <<  x << " ";
            else
                std::cout << std::setprecision(prec) << x << " ";
        std::cout << std::endl;
    }

    template<class T>
    void print(
        std::ofstream& of, const std::vector<T>& a, 
        YCI prec = 3, const bool& flag_scientific = false
    ){
        for(auto& x: a)
            if(flag_scientific)
                of << std::scientific << std::setprecision(prec) <<  x << " ";
            else
                of << std::setprecision(prec) <<  x << " ";
        of << "\n";
    }

    void print(std::vector<int> a);
    void print(short* a, const unsigned long& nr, const unsigned long& nc);

    std::string get_line(std::vector<int> a);

    std::string ltrim(YCS s);
    std::string rtrim(YCS s);
    std::string trim(YCS s);

    void insert_indent(YS line_original, YCS line_indent);

    bool is_present(YCVS v, YCS e);
    bool is_present(YCVI v, YCI e);

    void get_array_from_list(
        const std::list<std::vector<short>>& v, 
        short* array_1d, 
        const unsigned long& nr, 
        const unsigned long& nc
    );

    void get_current_date_time(YS line_date_time);
    struct H5File
    {
        /**
         * @brief Create and open an .hdf5 file with a name \p fname. 
         */
        void create(YCS fname);
        void close();

        inline void set_name(YCS fname){ name_ = fname; }

        /**
         * @brief Open an .hdf5 file with a name \p fname only to read it. 
         */
        void open_r();

        /**
         * @brief Open an .hdf5 file with a name \p fname to write-read it. 
         */
        void open_w();

        /**
         * @brief Add a group (folder) with a name \p gname to an already opened file.
         */
        void add_group(YCS gname);

        /**
         * @brief Add a dataset with a name \p dname, where a scalar \p v is to be written.
         * The dataset is put to a group \p gname.
         */
        template<class T>
        void add_scalar(const T& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            // add_group(gname);

            H5::Group grp(f_->openGroup(gname));
            write(v, dname, grp);
        }

        template<class T>
        void add_vector(const std::vector<T>& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            // add_group(gname);

            H5::Group grp(f_->openGroup(gname));
            write(v, dname, grp);
        }

        template<class T>
        void add_matrix(const std::list<std::vector<T>>& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            // add_group(gname);

            T* array_1d;
            unsigned long nr, nc;

            nr = v.size();
            nc = v.back().size();
            array_1d = new T[nr*nc];

            YMIX::get_array_from_list(v, array_1d, nr, nc);

            H5::Group grp(f_->openGroup(gname));
            write(array_1d, nc, nr, dname, grp);

            delete [] array_1d;
        }

        template<class T>
        void read_scalar(T& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to read a dataset " + dname + " from a group " + gname;
            H5::Group grp(f_->openGroup(gname));
            read(v, dname, grp);
        }

        template<class T>
        void read_vector(std::vector<T>& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            H5::Group grp(f_->openGroup(gname));
            read(v, dname, grp);
        }


        protected:
            inline void write(YCS v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                H5::StrType dtype(H5::PredType::C_S1, v.size()+1);
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(v, dtype);
            }
            inline void write(YCI v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_INT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((int*) &v, dtype);
            }
            inline void write(YCU v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_UINT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((unsigned*) &v, dtype);
            }
            inline void write(const long unsigned int& v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_ULONG;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((long unsigned int*) &v, dtype);
            }
            inline void write(const double& v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_DOUBLE;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((int*) &v, dtype);
            }

            inline void write(YCVU v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);
                auto dtype = H5::PredType::NATIVE_UINT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }
            inline void write(const std::vector<double>& v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);
                auto dtype = H5::PredType::NATIVE_DOUBLE;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }
            inline void write(const std::vector<Complex>& v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);

                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(Complex));
                H5Tinsert (dtype, "real", HOFFSET(Complex,real), H5T_NATIVE_DOUBLE);
                H5Tinsert (dtype, "imag", HOFFSET(Complex,imag), H5T_NATIVE_DOUBLE);

                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }

            inline void write(short* v, const unsigned long& nr, const unsigned long& nc, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {nr, nc};
                H5::DataSpace dspace(2, dims);
                auto dtype = H5::PredType::NATIVE_SHORT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(v, dtype);
            }


            template<class T>
            inline void read(T& v, YCS dname, H5::Group& grp)
            {
                H5::DataSet dataset = grp.openDataSet(dname);
                H5::DataType dtype = dataset.getDataType();
                dataset.read(&v, dtype);
            }
            inline void read(YS v, YCS dname, H5::Group& grp)
            {
                H5::DataSet dataset = grp.openDataSet(dname);
                H5::DataType dtype = dataset.getDataType();
                v="";
                dataset.read(v, dtype);
            }
            template<class T>
            inline void read(std::vector<T>& v, YCS dname, H5::Group& grp)
            {
                H5::DataSet dataset = grp.openDataSet(dname);

                H5::DataSpace dataspace = dataset.getSpace();
                int rank = dataspace.getSimpleExtentNdims();
                hsize_t dims_out[rank];
                int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

                unsigned long long N = 1;
                for(unsigned i_dim = 0; i_dim < rank; i_dim++)
                    N *= dims_out[i_dim];
                v = std::vector<T>(N);

                H5::DataType dtype = dataset.getDataType();

                dataset.read(&v[0], dtype, dataspace, dataspace);
            }

        protected:
            std::string name_;
            bool flag_opened;
            H5::H5File* f_;
    };
}


#endif