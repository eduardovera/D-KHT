#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <ctime>

class Logger {

    private:
        std::string taskname;
        std::chrono::high_resolution_clock::time_point start_time;

    public:
        Logger();
        Logger(const std::string taskname);
        Logger start_counting();
        void end_counting();
    };

#endif // LOGGER_H
