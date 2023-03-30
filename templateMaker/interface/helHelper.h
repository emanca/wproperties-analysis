#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

template <typename HIST_SF>
class helHelper
{

public:
    helHelper(HIST_SF &&hcoeff) : hcoeff_(std::make_shared<const HIST_SF>(std::move(hcoeff)))
    {
    }

    auto operator()(float y, float qt) const
    {

        auto const y_idx = hcoeff_->template axis<0>().index(y);
        auto const qt_idx = hcoeff_->template axis<1>().index(qt);

        ROOT::VecOps::RVec<float> coeff_array;
        coeff_array.resize(9);

        for (unsigned int i = 0; i < 9; i++)
        {
            coeff_array[i]=(hcoeff_->at(y_idx, qt_idx,i));
        }

        return coeff_array;
    }

protected:
    std::shared_ptr<const HIST_SF> hcoeff_;
};