#include <bits/stdc++.h>

using namespace std;
using ui = uint32_t;
using pii = pair<ui, ui>;

string get_file_name(string str)
{
    string ret = "";
    for (char ch : str)
    {
        ret += ch;
        if (ch == '\\' || ch == '/')
            ret = "";
    }
    return ret;
}

string get_file_name_without_suffix(string name)
{
    name = get_file_name(name);
    string ret = "";
    for (char ch : name)
    {
        if (ch == '.')
            break;
        ret += ch;
    }
    return ret;
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("input-graph-path are required !!! \n");
        exit(1);
    }
    string in_path = string(argv[1]);
    ifstream in(in_path);
    assert(in.is_open());
    ui n, m;
    in >> n >> m;
    vector<pii> edges;
    for (ui i = 0; i < m; i++)
    {
        ui a, b;
        in >> a >> b;
        if (a == b)
            continue;
        assert(a < n && b < n);
        edges.push_back({a, b});
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    m = edges.size();
    printf("input: %s , n= %u , m= %u\n", in_path.c_str(), n, m);
    ui *dout = new ui[n];
    memset(dout, 0, sizeof(ui) * n);
    ui *edge_to = new ui[m];
    for (ui i = 0; i < m; i++)
    {
        ui a = edges[i].first, b = edges[i].second;
        dout[a]++;
        edge_to[i] = b;
    }

    string out_file = get_file_name_without_suffix(in_path)+".bin";
    FILE *out = fopen(out_file.c_str(), "wb");
    assert(out != nullptr);
    ui size_of_int = sizeof(ui);
    fwrite(&size_of_int, sizeof(ui), 1, out);
    fwrite(&n, sizeof(ui), 1, out);
    fwrite(&m, sizeof(ui), 1, out);
    fwrite(dout, sizeof(ui), n, out);
    fwrite(edge_to, sizeof(ui), m, out);
    fflush(out);
    fclose(out);
    printf("output: %s\n", out_file.c_str());
    puts("");
    return 0;
}