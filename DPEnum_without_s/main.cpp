#include "Graph.h"
#include "Prepro.h"
#include "Branch.h"

string file_path;
string DPlex_path = "./DPlexes.out";
Graph g;
int paramK, paramL;

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        printf("3 params are required: $./enumdp graph_path k l !!! \n");
        exit(1);
    }
    file_path = string(argv[1]);
    paramK = atoi(argv[2]);
    paramL = atoi(argv[3]);
    g.readFromFile(file_path);

    Timer tot_time;

    // we just need to consider k<l
    if (paramK > paramL)
    {
        swap(paramK, paramL);
        g.flip();
    }

    Output out(DPlex_path);

#ifndef ONLY_LARGE
    puts("------------{preprocessing}------------");
    Prepro prepro(g, out, paramK, paramL);
    prepro.begin_preprocessing();
    printf("#DPlex: %lld , use-time: %.4lf s\n", out.counter, tot_time.get_time_seconds());
    printf("now: n= %u m= %u\n", g.n, g.m);
#endif

    puts("------------{Recursive Search}------------");
    Branch rec_search(g, out, paramK, paramL);
    rec_search.run();

    printf("Done: #DPlex= %lld use-time= %.4lf s dfs_cnt= %lld \n", out.counter, tot_time.get_time_seconds(), rec_search.dfs_cnt);

    return 0;
}