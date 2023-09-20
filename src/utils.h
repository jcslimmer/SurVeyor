#ifndef UTILS_H_
#define UTILS_H_

struct config_t {
    
    int threads;
    int min_clip_len, min_stable_mapq, min_diff_hsr;;

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
        min_diff_hsr = std::stoi(config_params["min_diff_hsr"]);
    }
};

struct stats_t {

    int read_len;
    int max_is;

    void parse(std::string stats_file) {
        std::unordered_map<std::string, std::string> stats_params;
        std::ifstream fin(stats_file);
        std::string name, value;
        while (fin >> name >> value) {
            stats_params[name] = value;
        }
        fin.close();

        read_len = std::stoi(stats_params["read_len"]);
        max_is = std::stoi(stats_params["max_is"]);
    }

};

struct contig_map_t {

    std::unordered_map<std::string, size_t> name_to_id;
    std::vector<std::string> id_to_name;

    void parse(std::string workdir) {
    	std::ifstream fin(workdir + "/contig_map");
		std::string name;
		int id = 0;
		while (fin >> name) {
			name_to_id[name] = id;
			id_to_name.push_back(name);
			id++;
		}
    }

    size_t size() {return id_to_name.size();}
    std::string get_name(size_t id) {return id_to_name[id];};
    size_t get_id(std::string& name) {return name_to_id[name];};
};

template<typename T>
T mean(std::vector<T>& v) {
    return std::accumulate(v.begin(), v.end(), (T)0.0)/v.size();
}

#endif
