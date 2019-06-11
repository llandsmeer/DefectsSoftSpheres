#ifndef BCC_OFFSETS_HPP
#define BCC_OFFSETS_HPP

#include "crystal.hpp"
#include "axis_offsets.hpp"

class bcc_offsets {
    crystal * crystalp;
    axis_offsets p1;
    axis_offsets p2;
    axis_offsets p3;
    axis_offsets p4;
    lattice_cell * mid;
    std::vector<axis_offsets*> sorted;

public:
    bcc_offsets(crystal * crystalp, lattice_cell * mid) :
        crystalp(crystalp),
        p1(crystalp, vec3(+1, +1, +1)),
        p2(crystalp, vec3(-1, +1, +1)),
        p3(crystalp, vec3(+1, -1, +1)),
        p4(crystalp, vec3(+1, +1, -1)),
        mid(mid) {
            p1.add_trace(mid);
            p2.add_trace(mid);
            p3.add_trace(mid);
            p4.add_trace(mid);
            p1.label = 1;
            p2.label = 2;
            p3.label = 3;
            p4.label = 4;
            sorted.push_back(&p1);
            sorted.push_back(&p2);
            sorted.push_back(&p3);
            sorted.push_back(&p4);
    }

    void measure() {
        p1.measure();
        p2.measure();
        p3.measure();
        p4.measure();
        std::sort(std::begin(sorted), std::end(sorted), [&](const auto & a, const auto & b) -> bool {
            return a->sum() < b->sum();
        });
    }

    void write(std::ostream & out) const {
        auto write_axis = [&](axis_offsets * ao) {
            int idx = 0;
            out << "A" << ao->label;
            for (const auto & ref : ao->particles) {
                out << " "
                    << ref.state
                    << " "
                    << ao->offsets[idx];
                idx += 1;
            }
        };
        write_axis(sorted[0]);
        out << " : ";
        write_axis(sorted[1]);
        out << " : ";
        write_axis(sorted[2]);
        out << " : ";
        write_axis(sorted[3]);
        out << std::endl;
    }
};

#endif
