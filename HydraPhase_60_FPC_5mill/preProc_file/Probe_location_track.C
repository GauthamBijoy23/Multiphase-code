#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <iomanip>
#include <sstream>

struct Cell {
    int id;         // Cell ID
    double x, y, z; // Cell center coordinates
    int proc;       // Processor ID
    double dist;    // Distance to probe
};

double distance(double x1, double y1, double z1,
                double x2, double y2, double z2) {
    return std::sqrt((x1 - x2)*(x1 - x2) +
                     (y1 - y2)*(y1 - y2) +
                     (z1 - z2)*(z1 - z2));
}

bool readGridFile(const std::string& filename, std::vector<Cell>& cells, int procId) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    int numCells;
    file >> numCells;
    cells.resize(numCells);

    for (int i = 0; i < numCells; ++i) {
        file >> cells[i].id >> cells[i].x >> cells[i].y >> cells[i].z;
        cells[i].proc = procId;
    }
    return true;
}

int main() {
    // Input probe location
    double probeX, probeY, probeZ;
    std::cout << "Enter probe location (x y z): ";
    std::cin >> probeX >> probeY >> probeZ;

    // Store the two nearest cells
    Cell nearest1{ -1, 0, 0, 0, -1, std::numeric_limits<double>::max() };
    Cell nearest2{ -1, 0, 0, 0, -1, std::numeric_limits<double>::max() };

    // Loop over processor files
    for (int proc = 0; proc < 100; ++proc) {
        std::ostringstream fname;
        fname << "geometry_files/grid_proc"
              << std::setw(4) << std::setfill('0') << proc
              << ".dat";

        std::vector<Cell> cells;
        if (!readGridFile(fname.str(), cells, proc)) {
            continue; // skip if file not found
        }

        for (auto& cell : cells) {
            cell.dist = distance(probeX, probeY, probeZ, cell.x, cell.y, cell.z);

            // Update nearest1 and nearest2
            if (cell.dist < nearest1.dist) {
                nearest2 = nearest1; // shift down
                nearest1 = cell;
            } else if (cell.dist < nearest2.dist) {
                nearest2 = cell;
            }
        }
    }

    // Print results
    if (nearest1.id != -1) {
        std::cout << "\nNearest cell 1: "
                  << "Processor " << nearest1.proc
                  << " | Cell ID = " << nearest1.id
                  << " | Distance = " << nearest1.dist << "\n";
    }
    if (nearest2.id != -1) {
        std::cout << "Nearest cell 2: "
                  << "Processor " << nearest2.proc
                  << " | Cell ID = " << nearest2.id
                  << " | Distance = " << nearest2.dist << "\n";
    }

    return 0;
}

