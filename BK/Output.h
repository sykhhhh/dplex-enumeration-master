#pragma once
#ifndef OUTPUT_H
#define OUTPUT_H

#include <sstream>

// the maximum number of solutions that we can output, 10^10 default
#define MAX_NUMBER 1e10

#ifdef DUMP_PLEX
// dump to file
class Output
{
    std::string file_path;
    using ui = uint32_t;
    using ll = long long;
    ll max_output_number;

public:
    std::ofstream out;
    ll counter; // the number of dumped DPlexes
    Output(std::string &s) : file_path(s), counter(0), max_output_number(MAX_NUMBER)
    {
        out.open(s);
        if (!out.is_open())
        {
            std::cout << "failed to open " << s << std::endl;
            exit(1);
        }
    }
    ~Output()
    {
        out.flush();
        out.close();
    }
    void dump(std::vector<ui> &s, int size = -1)
    {
        if (size == -1)
            size = s.size();
        counter++;
        for (int i = 0; i < size; i++)
        {
            ui u = s[i];
            out << u << ' ';
        }
        out << '\n';
        if (counter >= max_output_number)
        {
            out.flush();
            out.close();
            std::cout << "Too many solutions !!" << std::endl;
            exit(1);
        }
    }
};
#else
// just count the number
class Output
{
    using ui = uint32_t;
    using ll = long long;

public:
    ll counter; // the number of dumped DPlexes
    Output(std::string &s) : counter(0)
    {
    }
    ~Output()
    {
    }
    void dump(std::vector<ui> &s, int size = -1)
    {
        counter++;
    }
};
#endif

#endif