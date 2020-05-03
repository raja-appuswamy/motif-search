#include "header.h"

using namespace std;

Embedding::Embedding()
{
    srand(time(nullptr));

    for (int i = 0; i < N_ECHAR; i++)
        for (int j = 0; j < MAX_ELEN; j++)
            rnd_str[j * N_ECHAR + i] = rand() % 2;
}

string Embedding::embed_string(const string_view &str)
{
    unsigned len = str.size();
    unsigned elen = len * EFACTOR;
    string estr(elen, 0);
    unsigned i = 0;
    for (unsigned j = 0; j < elen; j++) {
        uint8_t s = i < len ? str[i] : EPAD;
        estr[j] = s;
        i = i + rnd_str[j * N_ECHAR + s];
    }

    return estr;
}

unsigned Embedding::embed_compare(const string_view &to_embed, const string &ref,
        unsigned threshold)
{
    unsigned len = to_embed.size();
    unsigned elen = len * EFACTOR;
    assert(elen == ref.size());

    unsigned i = 0, nmismatch = 0;
    for (unsigned j = 0; j < elen; j++)
    {
        uint8_t s = i < len ? to_embed[i] : EPAD;
        nmismatch += (ref[j] == s ? 0 : 1);
        if (nmismatch > threshold) {
            break;
        }
        i = i + rnd_str[j * N_ECHAR + s];
    }

    return nmismatch;
}
