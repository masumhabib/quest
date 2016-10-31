/** @page qlog Logging facility API documentation
 *
 * @section qlogintro Introduction
 *
 * This is a simple loggin facility for QUEST package. It has the following
 * features:
 * - Lightweight and fast: 
 *      - All TRACE level logs are ignored at compile time in Release build. 
 *      - Logging is either:
 *          - Asynchronous: uses dedicated thread for speedup.
 *          - Synchronous: if desired.
 *      - No file I/O is done by default: logs are written to std::log (stderr). 
 *        To save log in file, use OS features such as pipeling, ridirection 
 *        etc. 
 * - Simple API: 
 *      - Supports insertion operator <<.
 *      - Also supports formatted output style syntax like printf ().
 * - Output is configurable: Any object that has << operator defined can be 
 *   used as output. Default is std::log.
 * - In multiprocess (MPI/OpenMP) setup, it only prints logging from the 
 *   master by default.
 * - Five predefined log levels: TRACE, DEBUG, INFO, WARN and ERROR.
 * - Truncate repeating messages after printing predefined number of times.
 *
 * @section qlogexamples Examples
 *
 * 1. Complete simple example:
 * @code
 * #include <qlog.hpp>
 *
 * void test_qlog () {
 *     //qlog::init();
 *     //using qlog::log;
 *
 *     //log::info () << "Hey this will appear as information" << log::end;
 *
 *     //log::info ("Hey this is formatted info #%d", 0);
 * }
 *
 * @endcode
 *
 * 2. For other usage, please see the unit tests: 
 *
 */

#ifndef QLOG_HPP
#define QLOG_HPP

namespace qlog {

enum Verbosity {
    Trace = -2,
    Debug = -1, 
    Info  =  0,
    Warn  =  1,
    Error =  2,
};

//class qlog {
//
//public:
//    static Logger& get_logger ();
//    static Logger& create (bool threaded = true, std::size_t mpi_master = 0);
//private:
//    static std::unique_ptr<Logger> logger;
//};
//
//void init (bool threaded = true, std::size_t mpi_master = 0) {
//    Logger::create (threaded, mpi_master);
//}
//
//
//namespace log {
//
//Logger& info () {
//
//}
//
//}

}

#endif




