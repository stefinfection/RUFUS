//
// Created by Stephanie Georges on 8/13/24.
//
#include <iostream>
#include <string>
#include <array>

class Chunk {
private:
    std::string* chroms;
    static std::array<std::string, 24> grch38_chroms = {
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X",
            "Y"
    };
    std::string* lengths;
    static std::array<int, 24> grch38_chrom_lengths = {
            248956422,
            242193529,
            198295559,
            190214555,
            181538259,
            170805979,
            159345973,
            145138636,
            138394717,
            133797422,
            135086622,
            133275309,
            114364328,
            107043718,
            101991189,
            90338345,
            83257441,
            80373285,
            58617616,
            64444167,
            46709983,
            50818468,
            156040895,
            57227415
    };
public:
    // Constructor
    Chunk(int chunkSize, string species, string build) {
        this->chunkSize = chunkSize;
        if (build == "GRCh38" && species == "human") {
            this->chroms = grch38_chroms;
            this->lengths = grch38_chrom_lengths;
        } else {
            std::cout << "Genome not found" << std::endl;
        }
    }

    // Returns chromosome and coordinates for a given chunk
    std::string getChunk(int chunkNum) {
        int chunkStart = 0;
        int chunkEnd = 0;
        int chunkSize = this->chunkSize;
        for (int i = 0; i < this->chroms.size(); i++) {
            if (chunkNum < this->lengths[i] / chunkSize) {
                chunkStart = chunkNum * chunkSize;
                chunkEnd = chunkStart + chunkSize;
                return this->chroms[i] + ":" + std::to_string(chunkStart) + "-" + std::to_string(chunkEnd);
            } else {
                chunkNum -= this->lengths[i] / chunkSize;
            }
        }
        return "Chunk not found";
        -R chr${curr_chr}:${start_coord}-${end_coord}
    }
};