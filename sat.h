#ifndef SAT_H
#define SAT_H

#include <dlib/matrix.h>

using namespace std;
typedef dlib::matrix<double> Matrix;

class SAT {

    public:

    SAT() {
    }

        SAT(const int height, const int width) {
            this->SAT_X = Matrix(height, width);
            this->SAT_Y = Matrix(height, width);
            this->SAT_Z = Matrix(height, width);
            this->SAT_XX = Matrix(height, width);
            this->SAT_XY = Matrix(height, width);
            this->SAT_XZ = Matrix(height, width);
            this->SAT_YY = Matrix(height, width);
            this->SAT_YZ = Matrix(height, width);
            this->SAT_ZZ = Matrix(height, width);
            this->SAT_SAMPLES = Matrix(height, width);
        }

        Matrix SAT_X;
        Matrix SAT_Y;
        Matrix SAT_Z;
        Matrix SAT_XX;
        Matrix SAT_XY;
        Matrix SAT_XZ;
        Matrix SAT_YY;
        Matrix SAT_YZ;
        Matrix SAT_ZZ;
        Matrix SAT_SAMPLES;

        double getValueAt(const dlib::matrix<double> &M, const int &j, const int &i) {
            if (j < 0 || i < 0 || j >= M.nr() || i >= M.nc()) {
                return 0;
            }
            return M(j, i);
        }
};

#endif // SAT_H
