



#include <iostream>

#include <string>
#include <sstream>

#include <fstream>
#include <vector>

#include  <cmath>

#include <Eigen/Sparse>
#include <Eigen/Dense>


struct row_2_dat
{
    double x1;
    double x2;
};

std::vector <row_2_dat> read_f(std::string fname, char div)
{
    std::ifstream f(fname);
    std::vector <row_2_dat> data;  // vector containing the data
    std::string l; // line
    std::string item; //value inside the "cell"
    row_2_dat row;

    if (!f.is_open())
    {
        std::cout << "Could not open the file" << std::endl;
        return {};
    }
    else
    {
        std::getline(f,l); // Jumping the first line

        while (std::getline(f,l))
        {
            // checking if the line is empty
            if (l.empty()) continue;
            // end check

            std::stringstream ss(l); // line as string

            //std::cout << ss.str() << std::endl;  //   TESTS

            // column, for 2 data column
            std::getline(ss, item, div);
            row.x1 = std::stod(item);
            std::getline(ss, item, div);
            row.x2 = std::stod(item);

            data.push_back(row);
        }
    }
  
    return data;
};


std::vector <double> make_axis(double min_, double max_, int nnn_)
{
    // Returns the pF axis, in pF
    std::vector <double> axis(nnn_);
    double step = (max_ - min_) / (nnn_ - 1);
    double temp_;

    for (int i = 0 ; i < nnn_ ; i++)
    {
        temp_ = min_ + i*step;
        axis[i] = temp_ ;
    };

    return axis;
};


struct SamplePF
{
    double time;
    double pf;
};
struct BinnedSeries
{
    std::vector<double> y_grid; std::vector<double> w;
};
BinnedSeries bin_to_pf_axis_nearest_mean( std::vector<SamplePF> samples, std::vector<double> pf_axis)
{
    int n = pf_axis.size();
    double pf0 = pf_axis[0];
    double dpf = pf_axis[1] - pf_axis[0];

    BinnedSeries out;
    out.y_grid.assign(n, 0.0);
    out.w.assign(n, 0.0);

    std::vector<double> sum(n, 0.0);

    for (size_t k = 0; k < samples.size(); ++k) {
        SamplePF &s = samples[k];

        int i = (int)std::llround((s.pf - pf0) / dpf);
        if ((unsigned)i < (unsigned)n) {
            sum[i] += s.time;   // value being averaged
            out.w[i] += 1.0;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (out.w[i] > 0.0)
            out.y_grid[i] = sum[i] / out.w[i];
    }

    return out;
};



std::vector<double> whittaker_weighted_time(
    std::vector<double> y_grid,
    std::vector<double> w,
    double lambda)
{
    int n = (int)y_grid.size();

    // A = W + lambda * (D^T D), with 2nd difference penalty
    Eigen::SparseMatrix<double> A(n, n);
    std::vector<Eigen::Triplet<double>> T;
    T.reserve(9 * n);

    // Diagonal W
    for (int i = 0; i < n; ++i) {
        if (w[i] != 0.0) T.push_back(Eigen::Triplet<double>(i, i, w[i]));
    }

    // Add lambda * D^T D  (5-diagonal structure, built by accumulating row outer products)
    for (int k = 0; k < n - 2; ++k) {
        int i0 = k;
        int i1 = k + 1;
        int i2 = k + 2;

        double c0 =  1.0;
        double c1 = -2.0;
        double c2 =  1.0;

        // diagonal
        T.push_back(Eigen::Triplet<double>(i0, i0, lambda * c0 * c0));
        T.push_back(Eigen::Triplet<double>(i1, i1, lambda * c1 * c1));
        T.push_back(Eigen::Triplet<double>(i2, i2, lambda * c2 * c2));

        // off-diagonal symmetric
        T.push_back(Eigen::Triplet<double>(i0, i1, lambda * c0 * c1));
        T.push_back(Eigen::Triplet<double>(i1, i0, lambda * c1 * c0));

        T.push_back(Eigen::Triplet<double>(i0, i2, lambda * c0 * c2));
        T.push_back(Eigen::Triplet<double>(i2, i0, lambda * c2 * c0));

        T.push_back(Eigen::Triplet<double>(i1, i2, lambda * c1 * c2));
        T.push_back(Eigen::Triplet<double>(i2, i1, lambda * c2 * c1));
    }

    A.setFromTriplets(T.begin(), T.end());
    A.makeCompressed();

    // b = W*y
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) b[i] = w[i] * y_grid[i];

    // Solve
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd z = solver.solve(b);

    return std::vector<double>(z.data(), z.data() + n);
};

std::vector<row_2_dat> build_time_pf_dataset(
    std::vector<double> time_smoothed,
    std::vector<double> pf_axis)
{
    int n = (int)pf_axis.size();
    std::vector<row_2_dat> out;
    out.resize(n);

    for (int i = 0; i < n; ++i) {
        out[i].x1 = time_smoothed[i];
        out[i].x2 = pf_axis[i];
    }
    return out;
};

std::vector<std::vector<row_2_dat>> smooth_two_time_pf_datasets(
    std::vector<BinnedSeries> binned_all,
    std::vector<double> pf_axis,
    double lambda)
{
    int m = (int)binned_all.size();
    std::vector<std::vector<row_2_dat>> out;
    out.resize(m);

    for (int i = 0; i < m; ++i) {
        std::vector<double> z = whittaker_weighted_time(
            binned_all[i].y_grid,
            binned_all[i].w,
            lambda);

        out[i] = build_time_pf_dataset(z, pf_axis);
    }

    return out;
};

double average_time(double t1, double t2 , double w_aver)
{
    if ((w_aver < 0.000001) && (-0.000001 < w_aver) )
    {
        double aver_t = std::sqrt(t1 * t2);
        return aver_t;
    }
    double aver_t = std::pow(((std::pow(t1,w_aver) + std::pow(t2,w_aver))/2),(1/w_aver))  ;
    return aver_t;
};

void write_csv_pf_t1_t2_tavg(
    std::string out_name,
    std::vector<double> pf_axis,
    std::vector<row_2_dat> sm1,
    std::vector<row_2_dat> sm2,
    std::vector<row_2_dat> avg,
    char div)
{
    std::ofstream f(out_name);

    f << "pf" << div << "t1_smoothed" << div << "t2_smoothed" << div << "t_avg" << "\n";

    for (int i = 0; i < (int)pf_axis.size(); ++i)
    {
        f << pf_axis[i] << div
          << sm1[i].x1 << div
          << sm2[i].x1 << div
          << avg[i].x1 << "\n";
    };
};


// ##################### MAIN ##################### //


int main()
{

    // Minimum and maximum limits for h
    double min_h = 27;    // cm (positive!)
    // double max_h = 100000.0; // if not 6.8
    double max_h = std::pow(10.0, 5.8);  // pressure at water coontent 0

    int nnn_ = 100;    // quantity of dots in axis (divisions - 1) in pF

    double lambda = 10.0; // Lamda for whittaker smooth

    char div = ',';
    double height = 5.0; // cm
    double dist = 1.25;  // cm
    double delta = 0.2;  // cm

    double w_aver = 2.0; // w average 


    // Creating the h axis data
    double min_pf = std::log10(min_h);
    double max_pf = std::log10(max_h);
    std::vector <double> pf_axis = make_axis(min_pf, max_pf, nnn_);
    
    //std::vector <double> test = whitt(); // here Whittaken smoother
    
    for(int i=0; i<pf_axis.size(); i++)
    {
        //std::cout << pf_axis[i] << std::endl;
    //    std::cout << pf_axis[i] << std::endl;
    };

    

    std::vector <std::string> fnames_ = {"input_dat_t1.csv", "input_dat_t2.csv"};
    std::vector <std::vector <row_2_dat>> vec_dat_all(fnames_.size());
    std::vector <std::vector <row_2_dat>> vec_dat_log_all(fnames_.size());
    std::vector <double> positions = {};

    // Initially, for 2 tensiometers only
  

    // FIX HERE WHIT - TODO
    double aaa ;
    //aaa = whitt();
    //std::cout << aaa << std::endl;
    // Vector containing the data


    // READ DATA AND PUT INTO VECTOR h vector and pf vector
    for(int i=0; i < fnames_.size();i++)
    {
        std::vector <row_2_dat> vec_ = read_f(fnames_[i], div);
        vec_dat_all[i] = vec_ ;
        vec_dat_log_all[i] = vec_;
        for(int ii=0; ii < vec_.size(); ii++ )
        {
            vec_dat_log_all[i][ii].x2 = std::log10(vec_dat_all[i][ii].x2);
        };
    };


    std::vector<BinnedSeries> binned_all;
    binned_all.resize(fnames_.size());
    // converting 
    for (int i = 0; i < (int)fnames_.size(); ++i)
    {
        std::vector<SamplePF> samples;
        samples.resize(vec_dat_log_all[i].size());

        for (int k = 0; k < (int)vec_dat_log_all[i].size(); ++k)
        {
            samples[k].time = vec_dat_log_all[i][k].x1; // y = time
            samples[k].pf   = vec_dat_log_all[i][k].x2; // x = pf
        };

        binned_all[i] = bin_to_pf_axis_nearest_mean(samples, pf_axis);
    };

    
    std::vector<std::vector<row_2_dat>> smoothed = smooth_two_time_pf_datasets(binned_all, pf_axis, lambda);

    std::vector<row_2_dat> averaged;
    averaged.resize(smoothed[0].size());

    for (int i = 0; i < (int)smoothed[0].size(); ++i)
    {
        double pf = smoothed[0][i].x2;          // same as smoothed[1][i].x2
        double t1 = smoothed[0][i].x1;
        double t2 = smoothed[1][i].x1;

        averaged[i].x2 = pf;
        averaged[i].x1 = average_time(t1, t2, w_aver);
    };


    write_csv_pf_t1_t2_tavg("out_pf_time.csv", pf_axis, smoothed[0], smoothed[1], averaged, ',');
    

    return 0;
};






