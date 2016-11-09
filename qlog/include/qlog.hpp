/**@page qlog Logging facility API documentation
 *
 * @section qlogintro Introduction
 *
 * This is a simple loggin facility for QUEST package. It has the following
 * features:
 * - Lightweight and fast: 
 *     - All TRACE level logs are ignored at compile time in Release build. 
 *     - Early checking of severity -- no processing is done if the severity
 *       level is lower than the set value.
 *     - Logging is either:
 *         - Asynchronous: uses dedicated thread for speedup.
 *         - Synchronous: if desired.
 *     - No file I/O is done by default: logs are written to std::log (stderr). 
 * - Simple API: 
 *      - Supports insertion operator <<.
 *      - Also supports formatted output style syntax like printf ().
 * - Output is configurable: Any object that has << operator defined can be 
 *   used as output. Default is std::log.
 * - In multiprocess (MPI/OpenMP) setup, it only prints logging from the 
 *   master by default.
 * - Five predefined log levels: TRACE, DEBUG, INFO, WARN and ERROR.
 * - Truncate repeating messages after printing predefined number of times.
 * - Logging can be turned on/off for a particular class/file/module. 
 *
 * @section qlogexamples Examples
 *
 * 1. Complete simple example:
 *
 * @code
 * #include <qlog.hpp>
 *
 * // if logger object is available in the scope, it will be used
 * // by the LOG macros. Otherwise, the default logger object will be used.
 * auto logger = qlog::get_logger (_FILE_);
 *
 * void test_qlog () {
 *     LOG_INFO ("Hey this is formatted info #%d", 0);
 *     LOG_INFO () << "Hey this is formatted info #" << 0;
 *
 *     LOG_DEBUG ("Hey this is formatted info #%d", 0);
 *     LOG_DEBUG () << "Hey this is formatted info #" << 0;
 * }
 *
 * @endcode
 *
 * Output:
 * @code
 *   [II]: Hey this is a formatted info #0 
 *   [II]: Hey this is a formatted info #0 
 *   [DD]: [@qlog.hpp:38] Hey this is a formatted info #0 
 *   [DD]: [@qlog.hpp:39] Hey this is a formatted info #0 
 * @endcode
*
 * 2. For other usage, please see the unit tests.
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

enum class Severity {
    Trace = -2,
    Debug = -1, 
    Info  =  0,
    Warn  =  1,
    Error =  2,
};



}

#endif




