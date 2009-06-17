#ifndef _TIMER_HPP_
#define _TIMER_HPP_

//
// Copyright (c) 2008 Simen Kvaal
//
// This file is part of OpenFCI.
//
// OpenFCI is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// OpenFCI is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with OpenFCI. If not, see <http://www.gnu.org/licenses/>.
//

#include <ctime>
#include <iostream>
#include <iomanip>

/**
 * \file timer.hpp
 * \author Kenneth Wilder & Simen Kvaal
 * \date 9-16-08
 *
 * \brief Simple class to time actions. 
 *
 * Taken from http://oldmill.uchicago.edu/~wilder/Code/timer/ and documented using Doxygen by me.
 * I also added operator() to get the elapsed time... This was a hack in order to prevent
 * overloaded operator<<() to touch the stream's settings. So, 'cout << t()' puts a double to the
 * stream instead of the timer itself.
 *
 **/


/// \brief Class definition of timer.
class timer
{
    /// Support for streaming
    friend std::ostream& operator<<(std::ostream& os, timer& t);

  private:
    bool running;           ///< Is it running or not? Note that an instance must explicitely be started with timer::start() or timer::restart()
    clock_t start_clock;    ///< Clock when timer::start() is called.
    time_t start_time;      ///< Time when timer::start() is called.
    double acc_time;        ///< Accumulated time

    /// \brief Return elapsed time since timer::start() was called.
    double elapsed_time();

 public:
    ///\ brief Default constructor.
    timer() : running(false), start_clock(0), start_time(0), acc_time(0) { }

    void start(const char* msg = 0);
    void restart(const char* msg = 0);
    void stop(const char* msg = 0);
    void check(const char* msg = 0);
    /// \brief Get the time elapsed in seconds
    double operator()()
    {
      return acc_time + (running ? elapsed_time() : 0);
    }

}; // class timer

/// \brief Time elapsed in seconds
///
/// Return the total time that the timer has been in the "running"
/// state since it was first "started" or last "restarted".  For
/// "short" time periods (less than an hour), the actual cpu time
/// used is reported instead of the elapsed time.
inline double timer::elapsed_time()
{
  time_t acc_sec = time(0) - start_time;
  if (acc_sec < 3600)
    return (clock() - start_clock) / (1.0 * CLOCKS_PER_SEC);
  else
    return (1.0 * acc_sec);

} // timer::elapsed_time

/// \brief Start timer
///
/// Start a timer.  If it is already running, let it continue running.
/// Print an optional message.
/// \param msg Message to std::cout
inline void timer::start(const char* msg)
{
  // Print an optional message, something like "Starting timer t";
  if (msg) std::cout << msg << std::endl;

  // Return immediately if the timer is already running
  if (running) return;

  // Set timer status to running and set the start time
  running = true;
  start_clock = clock();
  start_time = time(0);

} // timer::start

/// \brief Restart timer
///
/// Turn the timer off and start it again from 0.  Print an optional message.
/// \param msg Message to std::cout
inline void timer::restart(const char* msg)
{
  // Print an optional message, something like "Restarting timer t";
  if (msg) std::cout << msg << std::endl;

  // Set timer status to running, reset accumulated time, and set start time
  running = true;
  acc_time = 0;
  start_clock = clock();
  start_time = time(0);

} // timer::restart


/// \brief Stop the timer and print an optional message.
/// \param msg Message to std::cout
inline void timer::stop(const char* msg)
{
  // Print an optional message, something like "Stopping timer t";
  if (msg) std::cout << msg << std::endl;

  // Compute accumulated running time and set timer status to not running
  if (running) acc_time += elapsed_time();
  running = false;

} // timer::stop


/// \brief Print out an optional message followed by the current timer timing.
/// \param msg Message to std::cout
inline void timer::check(const char* msg)
{
  // Print an optional message, something like "Checking timer t";
  if (msg) std::cout << msg << " : ";

  std::cout << "Elapsed time [" << std::setiosflags(std::ios::fixed)
            << std::setprecision(2)
            << acc_time + (running ? elapsed_time() : 0) << "] seconds\n";

} // timer::check


/// \brief Support for output to a stream
///
/// Allow timers to be printed to ostreams using the syntax 'os << t'
/// for an ostream 'os' and a timer 't'.  For example, "cout << t" will
/// print out the total amount of time 't' has been "running".
/// If you don't want to touch the stream's precision setting etc, use 'os << t()' instead.
inline std::ostream& operator<<(std::ostream& os, timer& t)
{
  os << std::setprecision(2) << std::setiosflags(std::ios::fixed)
    << t.acc_time + (t.running ? t.elapsed_time() : 0);
  return os;
}


#endif // _TIMER_HPP_

