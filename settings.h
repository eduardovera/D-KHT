#ifndef SETTINGS_H
#define SETTINGS_H


#include <limits>
#include <sat.h>
#include <dlib/geometry/vector.h>
#include <QDir>
#define NONZERO 0.00001

const double PI = acos(-1);
const double PI2 = 2.0 * PI;
static const double root22pi32 = 2.0 * sqrt(2.0) * pow(PI, 1.5);

typedef dlib::vector<double, 3> dvec3;
typedef dlib::vector<double, 3> dvec2;
typedef dlib::vector<int, 2> ivec2;


class hough_settings {
    public:

        /// COPY ROOM
        /**
        hough_settings() {
            max_point_distance = 0.0;
            input_filename = QDir::currentPath().toStdString() + "/input/copy_room_4161.png";
            output_filename = QDir::currentPath().toStdString() + "/output/output.png";
            max_distance2plane = 110;
            s_t = .2;
            n_phi = 153;
            n_rho = 130;
            s_ms = 200;

            inv_camera_fx = 1.0 / 526.37;
            inv_camera_fy = 1.0 / 526.37;
            camera_cx = 313.68;
            camera_cy = 259.02;
        }
        /**/

        /// CUBE
        /**/
        hough_settings() {

            max_point_distance = 0.0;
            input_filename = QDir::currentPath().toStdString() + "/input/cube.png";
            output_filename = QDir::currentPath().toStdString() + "/output/output.png";
            max_distance2plane = 150;
            s_t = .36;
            n_phi = 30;
            n_rho = 80;
            s_ms = 200;

            inv_camera_fx = 1.0 / 526.37;
            inv_camera_fy = 1.0 / 526.37;
            camera_cx = 313.68;
            camera_cy = 259.02;


        }
        /**/


        int n_phi;
        int n_rho;
        double s_t;
        int s_ms;


        double max_point_distance;
        double max_distance2plane;

        double inv_camera_fx;
        double inv_camera_fy;
        double camera_cx;
        double camera_cy;

        std::string input_filename;
        std::string output_filename;
};

#endif
