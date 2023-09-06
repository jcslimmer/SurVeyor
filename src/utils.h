#ifndef UTILS_H_
#define UTILS_H_

struct config_t {
    
    int threads;
    int min_clip_len, min_stable_mapq;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        threads = std::stoi(config_params["threads"]);
        min_clip_len = std::stoi(config_params["min_clip_len"]);
        min_stable_mapq = std::stoi(config_params["min_stable_mapq"]);
    }
};

#endif
