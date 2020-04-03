/* pfp-ds - prefix free parsing data structures
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file common.hpp
   \brief common.hpp contains common features.
   \author Massimiliano Rossi
   \date 12/03/2020
*/

#ifndef _COMMON_HH
#define _COMMON_HH

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <assert.h>

#include <sys/time.h>

#include <sys/mman.h> // for mmap
#include <unistd.h>
#include <sys/stat.h>
 #include <fcntl.h>

#include <sstream>      // std::stringstream

#include <vector>      // std::vector

#include <chrono>       // high_resolution_clock


std::string NowTime();
void _internal_messageInfo(const std::string message);
void _internal_messageWarning( const std::string file, const unsigned int line, const std::string message);
void _internal_messageError( const std::string file, const unsigned int line,const std::string message);


std::string NowTime()
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    char buffer[100];
    tm r;
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&tv.tv_sec, &r));
    char result[100];
    snprintf(result, 100, "%s"/*.%06ld"*/, buffer/*, (long)tv.tv_usec*/);
    return result;
}


template<typename T>
inline void _internal_message_helper(std::stringstream &ss, T const &first) { ss << first; }
template<typename T, typename... Args>
inline void _internal_message_helper(std::stringstream &ss, T const &first, const Args&... args) { ss << first << " "; _internal_message_helper(ss,args...); }
template<typename T, typename... Args>
inline std::string _internal_message(T const &first, const Args&... args) { std::stringstream ss; _internal_message_helper(ss,first,args...); return ss.str(); }


void _internal_messageInfo(const std::string message)
{
  std::cout << "[INFO] " << NowTime() << " - " << "Message: " << message << std::endl;
}

void _internal_messageWarning( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cout << "[WARNING] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
}

void _internal_messageError( const std::string file, const unsigned int line,
  const std::string message)
{
  std::cerr << "[ERROR] " << NowTime() << " - "
  << "File: " << file << '\n'
  << "Line: " << line << '\n'
  << "Message: " << message << std::endl;
  assert( false );
  exit( 1 );
}



#define info( args... ) \
    _internal_messageInfo( _internal_message(args) )

#ifdef VERBOSE
  #define verbose( args... ) \
      _internal_messageInfo( _internal_message(args) )
#else
  #define verbose( args... )
#endif

#define warning( args... ) \
    _internal_messageWarning( __FILE__, __LINE__, _internal_message(args) )

#define error( args... ) \
    _internal_messageError( __FILE__, __LINE__, _internal_message(args) )




//*********************** File I/O *********************************************
template<typename T>
void map_file(const char *filename, T*& ptr, size_t& length){
    struct stat filestat;
    int fd;

    if ((fd = open(filename, O_RDONLY)) < 0)
        error("open() file " + std::string(filename) + " failed" );

    if (fstat(fd, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    length = filestat.st_size / sizeof(T);

    if ((ptr = mmap(NULL, filestat.st_size, PROT_READ, MAP_SHARED, fd, 0)) == MAP_FAILED)
        error("mmap() file " + std::string(filename) + " failed");
}

template<typename T>
void read_file(const char *filename, T*& ptr, size_t& length){
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    length = filestat.st_size / sizeof(T);
    ptr = new T[length];

    if ((fread(ptr, sizeof(T), length, fd)) != length)
        error("fread() file " + std::string(filename) + " failed");

    fclose(fd);
}

template<typename T>
void read_file(const char *filename, std::vector<T>& ptr){
    struct stat filestat;
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    int fn = fileno(fd);
    if (fstat(fn, &filestat) < 0)
        error("stat() file " + std::string(filename) + " failed" );

    if(filestat.st_size % sizeof(T) != 0)
        error("invilid file " + std::string(filename));

    size_t length = filestat.st_size / sizeof(T);
    ptr.resize(length);

    if ((fread(&ptr[0], sizeof(T), length, fd)) != length)
        error("fread() file " + std::string(filename) + " failed");

    fclose(fd);
}

template<typename T>
void read_fasta_file(const char *filename, std::vector<T>& v){
    FILE* fd;

    if ((fd = fopen(filename, "r")) == nullptr)
        error("open() file " + std::string(filename) + " failed" );

    v.clear();

    char c;
    while (fread( &c, sizeof(char), 1,fd) == 1) {
      if(c == '>'){
        while(fread( &c, sizeof(char), 1,fd) == 1 && c != '\n');
      }else{
        v.push_back(c);
        while(fread( &c, sizeof(char), 1,fd) == 1 && c!= '\n') v.push_back(c);
      }
  	}
  	fclose(fd);
}


//*********************** Time resources ***************************************

/*!
 * op the operation that we want to measure
 */
#define elapsed_time(op) \
{ \
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now(); \
  op; \
  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now(); \
  verbose("Elapsed time (s): ",std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start ).count()); \
}

#endif /* end of include guard: _COMMON_HH */
