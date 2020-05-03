#pragma once
class Embedding
{
    public:
        Embedding();
        std::string embed_string(const std::string_view &str);
        unsigned embed_compare(const std::string_view &to_embed,
                const std::string &ref, unsigned threshold);

    private:
        constexpr static int N_ECHAR = 5;
        constexpr static int MAX_ELEN = 100;
        constexpr static int N_RNDBITS = MAX_ELEN * N_ECHAR;
        constexpr static int EFACTOR = 3;
        constexpr static int EPAD = 4;
        std::bitset<N_RNDBITS> rnd_str;
};
