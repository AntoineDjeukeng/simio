#include "middle_reservoir_monitor.hpp"
#include "reactor_setup.hpp"

#include "xtc.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void print_usage(const char* executable) {
    std::cerr << "Usage: " << executable
              << " <trajectory.xtc> <final.gro> <ion_insertion_report.json> <out_dir>"
                 " [max_frames=100] [frame_begin=0] [frame_end=-1]\n";
}

int parse_int(const char* value, const char* name) {
    try {
        return std::stoi(value);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid ") + name + ": " + value);
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5 || argc > 8) {
        print_usage(argv[0]);
        return 2;
    }

    try {
        const std::string xtc_path = argv[1];
        const std::string gro_path = argv[2];
        const std::string report_path = argv[3];
        const std::filesystem::path output_dir = argv[4];
        const int max_frames = argc > 5 ? parse_int(argv[5], "max_frames") : 100;
        const int frame_begin = argc > 6 ? parse_int(argv[6], "frame_begin") : 0;
        const int frame_end = argc > 7 ? parse_int(argv[7], "frame_end") : -1;
        if (max_frames == 0 || max_frames < -1 || frame_begin < 0 || (frame_end >= 0 && frame_end < frame_begin)) {
            throw std::runtime_error("Invalid frame range");
        }

        const simio::middle_reservoir::ReactorSetup setup =
            simio::middle_reservoir::load_reactor_setup(gro_path, report_path);
        std::filesystem::create_directories(output_dir);
        simio::middle_reservoir::write_setup_json(setup,
                                                   gro_path,
                                                   report_path,
                                                   (output_dir / "middle_reservoir_setup.json").string());

        std::cout << "[middle-reservoir] water=" << setup.nwater << " na=" << setup.nna
                  << " cl=" << setup.ncl << " natoms=" << setup.natoms << '\n'
                  << "[middle-reservoir] x_absolute=[" << setup.left_gate.x_nm << ", "
                  << setup.right_gate.x_nm << ") nm x_relative=[0, "
                  << setup.middle_length_nm() << ") nm\n"
                  << "[middle-reservoir] left_gate x=" << setup.left_gate.x_nm << " z=["
                  << setup.left_gate.zmin_nm << ',' << setup.left_gate.zmax_nm << ") type="
                  << setup.left_gate.wall_type << '\n'
                  << "[middle-reservoir] right_gate x=" << setup.right_gate.x_nm << " z=["
                  << setup.right_gate.zmin_nm << ',' << setup.right_gate.zmax_nm << ") type="
                  << setup.right_gate.wall_type << '\n';

        XtcTraj trajectory{};
        int rc = xtc_open(&trajectory, xtc_path.c_str(), 1);
        if (rc != exdrOK) throw std::runtime_error("Failed to open XTC, rc=" + std::to_string(rc));
        if (trajectory.natoms != setup.natoms) {
            xtc_close(&trajectory);
            throw std::runtime_error("XTC/GRO atom-count mismatch: XTC=" +
                                     std::to_string(trajectory.natoms) + " GRO=" +
                                     std::to_string(setup.natoms));
        }

        simio::Frame frame;
        frame.atoms.x.resize(static_cast<size_t>(setup.natoms));
        frame.atoms.y.resize(static_cast<size_t>(setup.natoms));
        frame.atoms.z.resize(static_cast<size_t>(setup.natoms));
        simio::middle_reservoir::MiddleReservoirMonitor monitor(setup);

        int global_frame = 0;
        int processed = 0;
        while ((rc = xtc_read_next(&trajectory)) == exdrOK) {
            const XtcFrame* source = xtc_tail(&trajectory);
            if (!source) throw std::runtime_error("XTC reader returned a null frame");
            const int current_frame = global_frame++;
            if (current_frame < frame_begin) continue;
            if (frame_end >= 0 && current_frame >= frame_end) break;
            if (max_frames >= 0 && processed >= max_frames) break;

            frame.step = source->step;
            frame.time_ps = source->time;
            frame.pbc = simio::Pbc3D(source->box[0][0], source->box[1][1], source->box[2][2]);
            for (int atom = 0; atom < setup.natoms; ++atom) {
                frame.atoms.x[static_cast<size_t>(atom)] = source->x[atom][0];
                frame.atoms.y[static_cast<size_t>(atom)] = source->x[atom][1];
                frame.atoms.z[static_cast<size_t>(atom)] = source->x[atom][2];
            }
            monitor.process_frame(frame, current_frame);
            ++processed;
        }
        xtc_close(&trajectory);
        if (rc != exdrOK && rc != exdrENDOFFILE) {
            throw std::runtime_error("XTC reader error, rc=" + std::to_string(rc));
        }
        if (processed == 0) throw std::runtime_error("No frames were processed");

        const std::filesystem::path csv_path = output_dir / "middle_reservoir.csv";
        monitor.write_csv(csv_path.string());
        std::cout << "Processed " << processed << " frame(s).\n"
                  << "  wrote: " << csv_path << '\n'
                  << "  wrote: " << (output_dir / "middle_reservoir_setup.json") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
