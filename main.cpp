#include <boost/program_options.hpp>
#include <fstream>
#include <array>
#include "header.h"

using namespace std;
static unsigned g_kmer_len = 3;
static uint8_t g_code[256];
static Embedding g_e;
uint32_t g_mask;

void make_code(void)
{
    for(int i = 0; i < 256; i++)
        g_code[i] = 4;

    g_code['A'] = g_code['a'] = 0;
    g_code['N'] = g_code['n'] = 0;
    g_code['C'] = g_code['c'] = 1;
    g_code['G'] = g_code['g'] = 2;
    g_code['T'] = g_code['t'] = 3;
}

struct Read
{
    string name, seq, qual;
    map<uint32_t, vector<unsigned>> index;
};

struct Motif
{
    string name, seq, eseq;
    unsigned pos;
};

istream &operator>>(istream &in, Read &r)
{
    string tmp;
    if (!getline(in, r.name)) return in;
    if (!getline(in, r.seq)) return in;
    if (!getline(in, tmp)) return in;
    if (!getline(in, r.qual)) return in;

    return in;
}

istream &operator>>(istream &in, Motif &m)
{
    if (!getline(in, m.name)) return in;
    return getline(in, m.seq);
}

void index_read(Read &r)
{
    unsigned len = r.seq.size();
    uint32_t k = 0;
    for (unsigned i = 0; i < g_kmer_len - 1; i++) {
        k = (k << 2) + *(g_code + r.seq[i]);
    }

    for (unsigned i = g_kmer_len - 1; i < len; i++) {
        k = (k << 2) + *(g_code + r.seq[i]);
        unsigned key = k & g_mask;
        r.index[key].push_back(i - (g_kmer_len - 1));
    }
}

void seed(const string &motif, Read &r, vector<unsigned> &candidates,
        unsigned threshold)
{
    unsigned len = motif.size();
    vector<unsigned> all_candidates;
    for (unsigned i = 0; i + g_kmer_len <= len; i += g_kmer_len) {
        uint32_t k = 0;
        for (unsigned j = i; j < i + g_kmer_len; j++)
            k = (k << 2) + *(g_code + motif[j]);
        assert(k == (k & g_mask));
        for (auto c = r.index[k].begin(); c != r.index[k].end(); ++c)
            all_candidates.push_back(*c > i ? *c - i : 0);
    }

    sort(all_candidates.begin(), all_candidates.end());
    unsigned ncandidates = all_candidates.size();
    for (unsigned i = 0; i < ncandidates; ) {
        unsigned j = i + 1;
        while (j < ncandidates &&
                all_candidates[i] == all_candidates[j])
            ++j;
        if (j - i >= threshold)
            candidates.push_back(all_candidates[i]);
        i = j;
    }
}

unsigned get_best_match(const string &motif, const string &emotif, Read &r,
        vector<unsigned> candidates)
{
    unsigned min_pos = numeric_limits<unsigned>::min();
    unsigned min_edist = numeric_limits<unsigned>::max();
    for (unsigned pos : candidates) {
        string_view cview(r.seq.data() + pos, motif.size());
        unsigned edist = g_e.embed_compare(cview, emotif, min_edist);
        if (edist < min_edist) {
            min_pos = pos;
            min_edist = edist;
        }
    }

    return min_pos;
}

void map_motifs(vector<Motif> &motifs, Read &r)
{
    for (unsigned i = 0; i < motifs.size(); i++) {
       vector<unsigned> candidates;
       seed(motifs[i].seq, r, candidates, 2);
       if (candidates.empty())
           seed(motifs[i].seq, r, candidates, 1);

       if (!candidates.empty()) {
           motifs[i].pos = get_best_match(motifs[i].seq, motifs[i].eseq,
                   r, candidates);
       } else {
           motifs[i].pos = numeric_limits<unsigned>::max();
       }
    }
}

void print_motifs(vector<Motif> &motifs, Read &r, ostream &out)
{
    sort(motifs.begin(), motifs.end(),
            [](Motif &left, Motif &right) {return left.pos < right.pos;});

    out << "Read: " << r.name << endl;
    for (Motif &m : motifs) {
        out << m.name << ",";
        if (m.pos == (numeric_limits<unsigned>::max()))
            out << "*\n";
        else
            out << m.pos << endl;
    }
}

void process_reads(vector<Motif> &motifs, const string &rfname, ostream &out)
{
    ifstream in(rfname);
    vector<Read> reads;
    Read r{};
    while (in >> r) {
        index_read(r);
        map_motifs(motifs, r);
        print_motifs(motifs, r, out);
    }
}

int main(int argc, char *argv[])
{
    namespace po = boost::program_options;
    bool help{};
    po::options_description description{"motif-search [options]"};
    description.add_options()
        ("help,h", po::bool_switch(&help), "Display help")
        ("motifs,m", po::value<string>(), "File containing motifs")
        ("read,r", po::value<vector<string>>()->multitoken(), "File(s) containing reads")
        ("kmerlen,l", po::value<unsigned>(), "Length of kmer (5)")
        ("output,o", po::value<string>(), "Ouptut file(default stdout)");
    po::command_line_parser parser{argc, argv};
    parser.options(description);
    auto parsed_result = parser.run();
    po::variables_map vm;
    po::store(parsed_result, vm);
    po::notify(vm);

    if (help) {
        cerr << description << endl;
        return 0;
    }

    if (vm["motifs"].empty() || vm["read"].empty()) {
        cout << "Must specificy motif and atleast one read file as input" << endl;
        return 1;
    }

    if (!vm["kmerlen"].empty()) {
        g_kmer_len = vm["kmerlen"].as<unsigned>();
        cout << "Using kmer of size " << g_kmer_len;
    }

    assert(g_kmer_len < 16);
    g_mask = (1U << (g_kmer_len * 2)) - 1;

    const string &motif_name = vm["motifs"].as<string>();
    cout << "Reading motif file " << motif_name << endl;
    ifstream mfile(motif_name);
    vector<Motif> motifs;
    Motif m;
    while (mfile >> m)
        motifs.push_back(m);
    cout << "Found " << motifs.size() << " motifs." << endl;

    for (Motif &m : motifs)
        m.eseq = g_e.embed_string(m.seq);
    
    make_code();

    string ofile_name = "";
    ofstream ofile;
    if (!vm["output"].empty()) {
        ofile_name = vm["output"].as<string>();
        ofile.open(ofile_name);
    }

    auto start = chrono::system_clock::now();
    for (const string &rf : vm["read"].as<vector<string>>()) {
        cout << "Processing read file " << rf << endl;
        process_reads(motifs, rf, ofile.is_open() ? ofile : cout);
    }
    if (ofile.is_open())
        ofile.close();
    auto end = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cerr << "Completed. Wall clock time: " <<
        elapsed.count() << " millisecs.\n";

    return 0;
}
