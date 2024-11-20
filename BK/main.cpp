#include "Graph.h"
#include "Branch.h"

string file_path;
string DPlex_path = "./DPlexes.out";
Graph g;
int paramK, paramL;
// int main()
// {
//     freopen("1.in","w",stdout);
//     g.readFromFile("Wiki-Vote.bin");
//     return 0;
// }
int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        printf("3 params are required: $./DPEM graph_path k l !!! \n");
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

    puts("------------{BK}------------");
    Branch bk(g, out, paramK, paramL);
    bk.run();

    printf("Done: #DPlex= %lld use-time= %.4lf s dfs_cnt= %lld \n", out.counter, tot_time.get_time_seconds(), bk.dfs_cnt);

    return 0;
}