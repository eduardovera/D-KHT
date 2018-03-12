#include "logger.h"

Logger::Logger() {
    this->taskname = "";
}

Logger::Logger(const std::string taskname)
{
    this->taskname = taskname;
}

Logger Logger::start_counting() {
    this->start_time = std::chrono::high_resolution_clock::now();
    return *this;
}

void Logger::end_counting() {
    auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time).count() * .001;
    std::cout << taskname << ": " << elapsed_time << " ms / " << 1.0 / (elapsed_time * .001) << " fps" << std::endl;
}
