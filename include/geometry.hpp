#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>

#include "../config/config.hpp"

class Geometry
{
  private:
    Config c;
    int nrec, nsrc;

    void load();
    void create();

    struct Receivers {
      std::vector<int> x;
      std::vector<int> z;
    };

    struct Sources {
      std::vector<int> x;
      std::vector<int> z;
    };

    Receivers createReceivers();
    Sources createSources();

  public:
    Receivers rec;
    Sources src;

    void get();
    void save();
};

#endif
