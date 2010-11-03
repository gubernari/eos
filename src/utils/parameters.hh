/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_PARAMETERS_HH
#define EOS_GUARD_SRC_UTILS_PARAMETERS_HH 1

#include <src/utils/exception.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    struct UnknownParameterError :
        public Exception
    {
        UnknownParameterError(const std::string & variable) throw ();
    };

    class Parameter;

    class Parameters :
        public PrivateImplementationPattern<Parameters>
    {
        private:
            Parameters(Implementation<Parameters> *);

        public:
            ~Parameters();

            Parameters clone() const;

            Parameter operator[] (const std::string & name) const;

            Parameter declare(const std::string & name, const double & value = 0.0);

            void set(const std::string & name, const double & value);

            static Parameters Defaults();
    };

    class Parameter
    {
        private:
            std::tr1::shared_ptr<Implementation<Parameters>> _imp;

            unsigned _index;

            Parameter(const std::tr1::shared_ptr<Implementation<Parameters>> & imp, unsigned index);

        public:
            friend class Parameters;

            ~Parameter();

            operator double () const;

            double operator() () const;

            const Parameter & operator= (const double &);

            const double & central() const;

            const double & max() const;

            const double & min() const;

            const std::string & name() const;
    };
}

#endif
